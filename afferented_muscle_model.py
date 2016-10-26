import numpy as np
from scipy.signal import sawtooth,square,gaussian,welch
import matplotlib.pyplot as plt
import time
import ipdb
from amm_functions import *

def afferented_muscle_model(muscle_parameters,delay_parameters,gain_parameters,TargetTrajectory,CorticalInput,**kwargs):
	"""
	muscle_parameters must be a dictionary with "Pennation Angle", "Muscle Mass",
	"Optimal Length", "Tendon Length", "Initial Muscle Length", and "Initial Tendon Length"
	values.

	delay_parameters must be a dictionary with "Efferent", "Ia","II", "Ib", and 
	"Cortical" values.

	gain_parameters must be a dictionary with "Gamma Dynamic Gain", "Gamma Static Gain",
	"Ia", "II", and "Ib" values.

	~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~

	FeedbackOption is a string that must be either 'ff_only' (feedforward only)
	'servo_control' (feedforward, spindle, and GTO), 'fb_control' (feedback
	control - proprioceptive system + supraspinal loop) and 'cortical_fb_only'

	"""

	FeedbackOption = kwargs.get("FeedbackOption",'ff_only')
	test_input_values(muscle_parameters,delay_parameters,gain_parameters,FeedbackOption)

	import numpy as np 
	from scipy import signal
	import control
	import random
	import matplotlib.pyplot as plt

	# muscle architectural parameters
	PennationAngle = muscle_parameters['Pennation Angle']
	MuscleMass = muscle_parameters['Muscle Mass']
	OptimalLength = muscle_parameters['Optimal Length']
	TendonLength = muscle_parameters['Tendon Length']
	OptimalTendonLength = TendonLength*1.05
	InitialMuscleLength = muscle_parameters['Initial Muscle Length']
	InitialTendonLength = muscle_parameters['Initial Tendon Length']
	InitialMusculoTendonLength = InitialMuscleLength*np.cos(PennationAngle)+InitialTendonLength

	# calculate initial length based on balance between stiffness of muscle and
	# tendon
	TendonConstant = 27.8
	TendonRateConstant = 0.0047
	NormalizedRestingTendonLength = 0.964

	MuscleConstant = 23.0  #355, 67.1
	MuscleRateConstant = 0.046 #0.04, 0.056
	RestingMuscleLength = 1.17  #1.35, 1.41
	PassiveMuscleForce = float(MuscleConstant * MuscleRateConstant * np.log(np.exp((1 - RestingMuscleLength)/MuscleRateConstant)+1))
	NormalizedSeriesElasticLength = float(TendonRateConstant*np.log(np.exp(PassiveMuscleForce/TendonConstant/TendonRateConstant)-1)+NormalizedRestingTendonLength)
	SeriesElasticLength = float(TendonLength * NormalizedSeriesElasticLength)

	MaximumMusculoTendonLength = float(OptimalLength * np.cos(PennationAngle) + TendonLength + 0.5)
	NormalizedMaximumFascicleLength = float((MaximumMusculoTendonLength - SeriesElasticLength)/OptimalLength)
	MaximumContractileElementLength = float(NormalizedMaximumFascicleLength/np.cos(PennationAngle)) # Contractile Element = Muscle

	InitialLength = float((InitialMusculoTendonLength+OptimalTendonLength\
							*((TendonRateConstant/MuscleRateConstant)*RestingMuscleLength-NormalizedRestingTendonLength\
							- TendonRateConstant*np.log((MuscleConstant/TendonConstant)*(MuscleRateConstant/TendonRateConstant))))\
							/(100*(1+TendonRateConstant/MuscleRateConstant*OptimalTendonLength/MaximumContractileElementLength*(1/OptimalLength))*np.cos(PennationAngle)))
	ContractileElementLength = float(InitialLength/(OptimalLength/100)) # initializing CE Length
	SeriesElasticElementLength = float((InitialMusculoTendonLength - InitialLength*100)/OptimalTendonLength) # MTU - initial muscle = initial SEE length

	# calculate maximal force output of the muscle based on PCSA
	MuscleDensity = 1.06 #g/cm**3
	PCSA = (MuscleMass*1000)/MuscleDensity/OptimalLength #physionp.logical cross-sectional area
	SpecificTension = 31.4 # WHERE DOES THIS COME FROM AND WHAT DOES IT MEAN? N/
	MaximumContractileElementForce = PCSA * SpecificTension

	# visnp.cosity of muscle (ref: Elias et al. 2014)
	MuscleViscosity = 0.005 #0.001

	# assign delays
	EfferentDelay = delay_parameters['Efferent Delay']
	IaAfferentDelay = delay_parameters['Ia Delay']
	IIAfferentDelay = delay_parameters['II Delay']
	IbAfferentDelay = delay_parameters['Ib Delay']
	CorticalDelay = delay_parameters['Cortical Delay']

	SamplingFrequency = 10000
	FeedbackSamplingFrequency = 1000
	SamplingRatio = SamplingFrequency/FeedbackSamplingFrequency
	SamplingPeriod = 1/SamplingFrequency
	Time = np.arange(0,len(TargetTrajectory)*SamplingPeriod,SamplingPeriod) # 0:SamplingPeriod:(length(TargetTrajectory)-1)/SamplingFrequency

	# parameter initialization
	ContractileElementVelocity = 0
	ContractileElementAcceleration = 0
	MuscleAcceleration = [0]
	MuscleVelocity = [0]
	MuscleLength = [ContractileElementLength*OptimalLength/100]
	SeriesElasticElementForce = 0.105

	# filter parameters for Noise
	BButtersCoefficients,AButtersCoefficients = signal.butter(4,100/(SamplingFrequency/2),'low')

	# Activation dynamics parameters
	Y_dot = 0
	Y = 0
	Saf_dot = 0
	Saf = 0
	fint_dot = 0
	fint = 0
	feff_dot = 0
	feff = 0

	ActivationFrequency = 0
	EffectiveMuscleActivation = 0

	DynamicSpindleFrequency = 0
	StaticSpindleFrequency = 0
	Bag1TensionSecondDeriv = 0
	Bag1TensionFirstDeriv = 0
	Bag1Tension = 0.0
	Bag2TensionSecondDeriv = 0
	Bag2TensionFirstDeriv = 0
	Bag2Tension = 0
	ChainTensionSecondDeriv = 0
	ChainTensionFirstDeriv = 0
	ChainTension = 0

	## GTO model
	GTOConstant1 = 60
	GTOConstant2 = 4

	Num,Den = [1.7,2.58,0.4],[1,2.2,0.4]
	ContinuousTransferFunction = control.tf(Num,Den)
	DiscreteTransferFunction = control.matlab.c2d(ContinuousTransferFunction,1/FeedbackSamplingFrequency)
	Num,Den = control.matlab.tfdata(DiscreteTransferFunction)
	Num,Den = Num[0][0],Den[0][0]

	## Transcortical loop
	TransCorticalLoopConstant = 0.01

	## Gains
	GammaDynamicGain = gain_parameters['Gamma Dynamic Gain']
	GammaStaticGain = gain_parameters['Gamma Static Gain']
	IaSpindleGain = gain_parameters['Ia Gain']
	IISpindleGain = gain_parameters['II Gain']
	IbGTOGain = gain_parameters['Ib Gain']

	# Convert delays in ms to samples
	IaAfferentDelayTimeStep = int(IaAfferentDelay*SamplingRatio)   #Ia + II 30
	IIAfferentDelayTimeStep = int(IIAfferentDelay*SamplingRatio)
	IbAfferentDelayTimeStep = int(IbAfferentDelay*SamplingRatio)   #Ib 40
	CorticalDelayTimeStep = int(CorticalDelay*SamplingRatio)   #cortical 50
	EfferentDelayTimeStep = int(EfferentDelay*SamplingRatio)

	FeedbackInput = 0
	Count = 0

	# Convert force trajectory to unit of newton
	TargetForceTrajectory = TargetTrajectory*MaximumContractileElementForce

	def bag1_model(CE,Bag1):
		Length,LengthFirstDeriv,LengthSecondDeriv = CE['Length'][-1],CE['Velocity'][-1],CE['Acceleration'][-1]
		DynamicSpindleFrequency, GammaDynamicGain = Bag1['DynamicSpindleFrequency'][-1], Bag1['GammaDynamicGain']
		Bag1TensionFirstDeriv, Bag1Tension = Bag1['TensionFirstDeriv'][-1], Bag1['Tension'][-1]
		## Feedback system parameters
		## Spindle Model
		p = 2

		Bag1TimeConstant = 0.149
		Bag1Frequency = 60

		Bag1BetaNot = 0.0605
		Bag1Beta = 0.2592
		Bag1Gamma = 0.0289

		G = 20000  #7000

		C = 0.58*(LengthFirstDeriv >= 0) + 0.42 # when true C = 1 else, C = 0.42       
		
		DynamicSpindleFrequencyFirstDeriv = (GammaDynamicGain**p/(GammaDynamicGain**p+Bag1Frequency**p) - DynamicSpindleFrequency)\
										/ Bag1TimeConstant
		DynamicSpindleFrequency = (SamplingPeriod)*DynamicSpindleFrequencyFirstDeriv + DynamicSpindleFrequency
		
		Bag1Beta = Bag1BetaNot + Bag1Beta * DynamicSpindleFrequency
		Bag1Gamma = Bag1Gamma * DynamicSpindleFrequency

		R = 0.46 #length dependency of the force-velocity relationship
		a = 0.3
		K_SR = 10.4649
		K_PR = 0.15
		M = 0.0002 # intrafusal fiber mass

		LN_SR = 0.0423

		L0_SR = 0.04 #in units of L0
		L0_PR = 0.76 #polar region rest length

		Bag1TensionSecondDeriv = float((K_SR/M)*(	(C*Bag1Beta*np.sign(LengthFirstDeriv-Bag1TensionFirstDeriv/K_SR) \
											* ((abs(LengthFirstDeriv-Bag1TensionFirstDeriv/K_SR))**a)  \
											* (Length-L0_SR-Bag1Tension/K_SR-R)) \
											+ K_PR*(Length-L0_SR-Bag1Tension/K_SR-L0_PR) \
											+ M*LengthSecondDeriv \
											+ Bag1Gamma \
											- Bag1Tension 	))
		Bag1['TensionFirstDeriv'].append(float(Bag1TensionSecondDeriv*(SamplingPeriod) + Bag1TensionFirstDeriv))
		Bag1['Tension'].append(float(Bag1['TensionFirstDeriv'][-1]*(SamplingPeriod) + Bag1Tension))
		Bag1['DynamicSpindleFrequency'].append(DynamicSpindleFrequency)
		Bag1AfferentPotential = G*(Bag1['Tension'][-1]/K_SR-(LN_SR-L0_SR))

		return(Bag1AfferentPotential)
	def bag2_model(CE,Bag2):
		Length,LengthFirstDeriv,LengthSecondDeriv = CE['Length'][-1],CE['Velocity'][-1],CE['Acceleration'][-1]
		StaticSpindleFrequency, GammaStaticGain = Bag2['StaticSpindleFrequency'][-1], Bag2['GammaStaticGain']
		Bag2TensionFirstDeriv, Bag2Tension = Bag2['TensionFirstDeriv'][-1], Bag2['Tension'][-1]

		## Feedback system parameters
		## spindle_model Model
		p = 2

		Bag2TimeConstant = 0.205
		Bag2Frequency = 60

		Bag2BetaNot = 0.0822
		Bag2Beta = -0.046
		Bag2Gamma = 0.0636

		if LengthFirstDeriv >= 0:
			C = 1 #constant describing the np.experimentally observed asymmetric effect of velocity on force production during lengthening and shortening
		else:
			C = 0.42

		G = 10000 #7250 #3800

		StaticSpindleFrequencyFirstDeriv = (GammaStaticGain**p/(GammaStaticGain**p+Bag2Frequency**p)-StaticSpindleFrequency) \
										/ Bag2TimeConstant
		StaticSpindleFrequency = (SamplingPeriod)*StaticSpindleFrequencyFirstDeriv + StaticSpindleFrequency

		Bag2Beta = Bag2BetaNot + Bag2Beta * StaticSpindleFrequency
		Bag2Gamma = Bag2Gamma * StaticSpindleFrequency

		R = 0.46 #length dependency of the force-velocity relationship
		a = 0.3
		K_SR = 10.4649
		K_PR = 0.15
		M = 0.0002 # intrafusal fiber mass

		LN_SR = 0.0423
		LN_PR = 0.89

		L0_SR = 0.04 #in units of L0
		L0_PR = 0.76 #polar region rest length

		L_secondary = 0.04
		X = 0.7

		Bag2TensionSecondDeriv = (K_SR/M)*(	(C*Bag2Beta*np.sign(LengthFirstDeriv-Bag2TensionFirstDeriv/K_SR)  \
											*((abs(LengthFirstDeriv-Bag2TensionFirstDeriv/K_SR))**a)  \
											*(Length-L0_SR-Bag2Tension/K_SR-R)+K_PR  \
											*(Length-L0_SR-Bag2Tension/K_SR-L0_PR))  \
											+ M*LengthSecondDeriv \
											+ Bag2Gamma \
											- Bag2Tension )
		Bag2['TensionFirstDeriv'].append(Bag2TensionSecondDeriv*(SamplingPeriod) + Bag2TensionFirstDeriv)
		Bag2['Tension'].append(Bag2['TensionFirstDeriv'][-1]*(SamplingPeriod) + Bag2Tension)
		Bag2['StaticSpindleFrequency'].append(StaticSpindleFrequency)
		Bag2PrimaryAfferentPotential = G*(Bag2['Tension'][-1]/K_SR-(LN_SR-L0_SR))
		Bag2SecondaryAfferentPotential = G*(	X*(L_secondary/L0_SR)*(Bag2['Tension'][-1]/K_SR-(LN_SR-L0_SR)) \
							+(1-X)*(L_secondary/L0_PR)*(Length-Bag2['Tension'][-1]/K_SR-(L0_SR+LN_PR))	)

		return(Bag2PrimaryAfferentPotential,Bag2SecondaryAfferentPotential)
	def chain_model(CE,Chain,Bag2):
		Length,LengthFirstDeriv,LengthSecondDeriv = CE['Length'][-1],CE['Velocity'][-1],CE['Acceleration'][-1]
		ChainTensionFirstDeriv,ChainTension = Chain['TensionFirstDeriv'][-1],Chain['Tension'][-1]
		StaticSpindleFrequency,GammaStaticGain = Bag2['StaticSpindleFrequency'][-1],Bag2['GammaStaticGain']
		## Feedback system parameters
		## spindle_model Model
		p = 2

		ChainFrequency = 90

		ChainBetaNot = 0.0822
		ChainBeta = - 0.069
		ChainGamma = 0.0954

		if LengthFirstDeriv >= 0:
			C = 1 #constant describing the np.experimentally observed asymmetric effect of velocity on force production during lengthening and shortening
		else:
			C = 0.42

		G = 10000 #7250    #3000

		ChainStaticSpindleFrequency = GammaStaticGain**p/(GammaStaticGain**p+ChainFrequency**p)

		ChainBeta = ChainBetaNot + ChainBeta * ChainStaticSpindleFrequency
		ChainGamma = ChainGamma * StaticSpindleFrequency

		R = 0.46 #length dependency of the force-velocity relationship
		a = 0.3
		K_SR = 10.4649
		K_PR = 0.15
		M = 0.0002 # intrafusal fiber mass

		LN_SR = 0.0423
		LN_PR = 0.89

		L0_SR = 0.04 #in units of L0
		L0_PR = 0.76 #polar region rest length

		L_secondary = 0.04
		X = 0.7

		ChainTensionSecondDeriv = K_SR/M * (C * ChainBeta * np.sign(LengthFirstDeriv-ChainTensionFirstDeriv/K_SR)*((abs(LengthFirstDeriv-ChainTensionFirstDeriv/K_SR))**a)*(Length-L0_SR-ChainTension/K_SR-R)+K_PR*(Length-L0_SR-ChainTension/K_SR-L0_PR)+M*LengthSecondDeriv+ChainGamma-ChainTension)
		Chain['TensionFirstDeriv'].append(ChainTensionSecondDeriv*SamplingPeriod + ChainTensionFirstDeriv)
		Chain['Tension'].append(ChainTensionFirstDeriv*SamplingPeriod + ChainTension)

		ChainPrimaryAfferentPotential = G*(Chain['Tension'][-1]/K_SR-(LN_SR-L0_SR))
		ChainSecondaryAfferentPotential = G*(	X*(L_secondary/L0_SR)*(Chain['Tension'][-1]/K_SR-(LN_SR-L0_SR)) \
							+ (1-X)*(L_secondary/L0_PR)*(Length-Chain['Tension'][-1]/K_SR-(L0_SR+LN_PR))	)

		return(ChainPrimaryAfferentPotential,ChainSecondaryAfferentPotential)
	def spindle_model(CE,Bag1,Bag2,Chain):

		S = 0.156
		Bag1AfferentPotential = bag1_model(CE,Bag1)
		Bag2PrimaryAfferentPotential,Bag2SecondaryAfferentPotential = bag2_model(CE,Bag2)
		ChainPrimaryAfferentPotential,ChainSecondaryAfferentPotential = chain_model(CE,Chain,Bag2)

		if Bag1AfferentPotential < 0: Bag1AfferentPotential = 0 
		if Bag2PrimaryAfferentPotential < 0: Bag2PrimaryAfferentPotential = 0
		if ChainPrimaryAfferentPotential < 0: ChainPrimaryAfferentPotential = 0
		if Bag2SecondaryAfferentPotential < 0: Bag2SecondaryAfferentPotential = 0
		if ChainSecondaryAfferentPotential < 0: ChainSecondaryAfferentPotential = 0	        

		if Bag1AfferentPotential > (Bag2PrimaryAfferentPotential+ChainPrimaryAfferentPotential):
			Larger = Bag1AfferentPotential
			Smaller = Bag2PrimaryAfferentPotential+ChainPrimaryAfferentPotential
		else:
			Larger = Bag2PrimaryAfferentPotential+ChainPrimaryAfferentPotential
			Smaller = Bag1AfferentPotential

		PrimaryOutput = Larger + S * Smaller
		SecondaryOutput = Bag2SecondaryAfferentPotential + ChainSecondaryAfferentPotential

		if PrimaryOutput < 0:
			PrimaryOutput = 0
		elif PrimaryOutput > 100000:
			PrimaryOutput = 100000
		if SecondaryOutput < 0:
			SecondaryOutput = 0
		elif SecondaryOutput > 100000:
			SecondaryOutput = 100000

		return(PrimaryOutput,SecondaryOutput)
	def activation_frequency_slow(CE,SlowTwitch,Activation,SamplingFrequency):
		Length, LengthFirstDeriv = CE['Length'][-1], CE['Velocity'][-1]
		Y,fint,feff_dot,feff = SlowTwitch['Y'],SlowTwitch['fint'],SlowTwitch['feff_dot'],SlowTwitch['feff']
		ActivationFrequencySlow = SlowTwitch['ActivationFrequency']

		Uth = 0.001
		f_half = 8.5 #f_half = 34 #can be found in Table 2 of Brown and Loeb 2000
		fmin = 4 #fmin = 15
		fmax = 2*f_half #fmin = 0.5*f_half #can be found in Song et al. 2008
		af = 0.56
		nf0 = 2.11
		nf1 = 5
		cy = 0.35
		Vy = 0.1
		Ty = 0.2

		Tf1 = 0.0484
		Tf2 = 0.032
		Tf3 = 0.0664
		Tf4 = 0.0356

		Y_dot = (1 - cy*(1-np.exp(-abs(LengthFirstDeriv)/Vy))-Y)/Ty
		Y = Y_dot*SamplingPeriod + Y 

		#nan_test(Y_dot,'Y_dot')
		#nan_test(Y,'Y')

		fenv = (fmax-fmin)/(1-Uth)*Activation+fmin-(fmax-fmin)*Uth
		fenv = fenv/f_half

		#nan_test(fenv,'fenv')

		if feff_dot >= 0:
		    Tf = Tf1 * Length**2 + Tf2 * fenv
		else:
		    Tf = (Tf3 + Tf4*ActivationFrequencySlow)/Length

		#nan_test(Tf,'Tf')

		fint_dot = (fenv - fint)/Tf
		fint = fint_dot*SamplingPeriod + fint
		feff_dot = (fint - feff)/Tf
		feff = feff_dot*SamplingPeriod + feff

		#nan_test(fint_dot,'fint_dot')
		#nan_test(fint,'fint')
		#nan_test(feff_dot,'feff_dot')
		#nan_test(feff,'feff')

		nf = nf0 + nf1*(1/Length-1)
		SlowTwitch['ActivationFrequency'] = 1 - np.exp(-(Y*feff/(af*nf))**nf)
		SlowTwitch['Y'],SlowTwitch['fint'],SlowTwitch['feff_dot'],SlowTwitch['feff'] = Y,fint,feff_dot,feff 
		return(SlowTwitch)
	def activation_frequency_fast(CE,FastTwitch,Activation,SamplingFrequency):
		Length = CE['Length'][-1]
		Saf,fint,feff_dot,feff = FastTwitch['Saf'],FastTwitch['fint'],FastTwitch['feff_dot'],FastTwitch['feff']
		ActivationFrequencyFast = FastTwitch['ActivationFrequency']

		Uth = 0.001
		f_half = 34 #f_half = 34 #can be found in Table 2 of Brown and Loeb 2000
		fmin = 15 #fmin = 15
		fmax = 2*f_half #fmin = 0.5*f_half #can be found in Song et al. 2008
		af = 0.56
		nf0 = 2.11
		nf1 = 3.3
		as1 = 1.76
		as2 = 0.96
		Ts = 0.043

		Tf1 = 0.0206
		Tf2 = 0.0136
		Tf3 = 0.0282
		Tf4 = 0.0151

		fenv = (fmax-fmin)/(1-Uth)*Activation+fmin-(fmax-fmin)*Uth
		fenv = fenv/f_half

		#nan_test(fenv,'fenv')

		if feff_dot >= 0:
			Tf = Tf1 * Length**2 + Tf2 * fenv
		elif feff_dot < 0:
			Tf = (Tf3 + Tf4*ActivationFrequencyFast)/Length

		#nan_test(Tf,'Tf')

		if feff < 0.1:
			AS = as1
		elif feff >= 0.1:
			AS = as2

		#nan_test(AS,'AS')

		Saf_dot = (AS - Saf)/Ts
		Saf = Saf_dot*SamplingPeriod + Saf
		fint_dot = (fenv - fint)/Tf
		fint = fint_dot*SamplingPeriod + fint
		feff_dot = (fint - feff)/Tf
		feff = feff_dot*SamplingPeriod + feff
		nf = nf0 + nf1*(1/Length-1)
		FastTwitch['ActivationFrequency'] = 1 - np.exp(-(S_af*feff/(af*nf))**nf)
		FastTwitch['Saf'],FastTwitch['fint'],FastTwitch['feff_dot'],FastTwitch['feff'] = Saf,fint,feff_dot,feff
		return(FastTwitch)
	def force_length(CE):
		Length = CE['Length'][-1]
		beta = 2.3
		omega = 1.12
		rho = 1.62
		ForceLength = np.exp(-abs((Length**beta - 1)/omega)**rho)
		return(ForceLength)
	def concentric_force_velocity(CE):
		Length,V = CE['Length'][-1],CE['Velocity'][-1]
		Vmax = -7.88
		cv0 = 5.88
		cv1 = 0
		ConcentricForceVelocity = (Vmax - V)/(Vmax + (cv0 + cv1*Length)*V)
		return(ConcentricForceVelocity)
	def eccentric_force_velocity(CE):
		Length,V = CE['Length'][-1],CE['Velocity'][-1]
		av0 = -4.7
		av1 = 8.41
		av2 = -5.34
		bv = 0.35
		EccentricForceVelocity = (bv - (av0 + av1*Length + av2*Length**2)*V)/(bv+V)
		return(EccentricForceVelocity)
	def parallel_elastic_element_force_1(CE):
		Length, MaximumContractileElementLength = CE['Length'][-1], CE['Maximum Length']
		c1_pe1 = 23.0  #355, 67.1
		k1_pe1 = 0.046 #0.04, 0.056
		Lr1_pe1 = 1.17  #1.35, 1.41
		ParallelElasticElementForce1 = c1_pe1 * k1_pe1 * np.log(np.exp((Length/MaximumContractileElementLength - Lr1_pe1)/k1_pe1)+1)
		return(ParallelElasticElementForce1)
	def parallel_elastic_element_force_2(CE):
		Length = CE['Length'][-1]
		c2_pe2 = -0.02 #0.01  -0.1
		k2_pe2 = -21
		Lr2_pe2 = 0.70 #0.79 0.59
		ParallelElasticElementForce2 = c2_pe2*np.exp((k2_pe2*(Length-Lr2_pe2))-1)
		return(ParallelElasticElementForce2)
	def normalized_series_elastic_element_force(SEE):
		LT = SEE['Length'][-1]
		cT_se = 27.8
		kT_se = 0.0047
		LrT_se = 0.964
		NormalizedSeriesElasticElementForce = cT_se * kT_se * np.log(np.exp((LT - LrT_se)/kT_se)+1)
		return(NormalizedSeriesElasticElementForce)
	
	Input, IaInput, IIInput = [],[],[]
	IbInput, x, TemporaryIbInput = [],[],[]
	Noise, FilteredNoise, LongInput = [],[],[]
	OutputForceMuscle,OutputForceTendon,OutputForceLength = [], [], []
	OutputForceVelocity,OutputForcePassive1,OutputForcePassive2 = [], [], []
	OutputSeriesElasticElementLength,OutputContractileElementVelocity = [(InitialMusculoTendonLength - InitialLength*100)/OptimalTendonLength], []
	OutputContractileElementLength = [InitialLength/(OptimalLength/100)]
	OutputContractileElementAcceleration,OutputActivationFrequency,OutputEffectiveMuscleActivation = [], [], []
	
	CE = initialize_dictionary(['Length','Velocity','Acceleration','Maximum Length'],\
									[[InitialLength/(OptimalLength/100)],[0],[0],float(NormalizedMaximumFascicleLength/np.cos(PennationAngle))])
	Bag1 = initialize_dictionary(['GammaDynamicGain','DynamicSpindleFrequency','Tension','TensionFirstDeriv'],\
										[GammaDynamicGain,[DynamicSpindleFrequency],[Bag1Tension],[Bag1TensionFirstDeriv]])
	Bag2 = initialize_dictionary(['GammaStaticGain','StaticSpindleFrequency','Tension','TensionFirstDeriv'],\
										[GammaStaticGain,[StaticSpindleFrequency],[Bag2Tension],[Bag2TensionFirstDeriv]])
	Chain = initialize_dictionary(['Tension','TensionFirstDeriv'],\
										[[ChainTension],[ChainTensionFirstDeriv]])
	SlowTwitch = initialize_dictionary(['Y','fint','feff_dot','feff','ActivationFrequency'],\
											[Y,fint,feff_dot,feff,ActivationFrequency])
	FastTwitch = initialize_dictionary(['Saf','fint','feff_dot','feff','ActivationFrequency'],\
											[Saf,fint,feff_dot,feff,ActivationFrequency])
	SEE = initialize_dictionary(['Length','Force'],[[float((InitialMusculoTendonLength - InitialLength*100)/OptimalTendonLength)],[SeriesElasticElementForce]])

	StartTime = time.time()
	FeedforwardInput = TargetForceTrajectory/MaximumContractileElementForce

	for i in range(len(Time)): #= 1:length(Time)
		if FeedbackOption == 'ff_only': # Feedforward input onl
			Input.append(FeedforwardInput[i])
			IaInput.append(0)
			IIInput.append(0)
			IbInput.append(0)
		elif FeedbackOption == 'servo_control': # Servo control (feedforward + spindle and GTO)
			PrimaryOutput,SecondaryOutput = spindle_model(CE,Bag1,Bag2,Chain)
			IaInput.append(PrimaryOutput)
			IIInput.append(SecondaryOutput)
			if i == Count:
				# Get spindle primary and secondary afferent activity

				# Get Ib activity
				x.append(GTOConstant1*np.log(SEE['Force'][-1]/GTOConstant2+1))
				if i == 0:
					TemporaryIbInput.append(x[-1])
				elif i == 1*SamplingRatio:
					TemporaryIbInput.append(x[-1])
				elif i == 2*SamplingRatio:
					TemporaryIbInput.append((Num[0]*x[-1])/Den[0])
				elif i == 3*SamplingRatio:
					TemporaryIbInput.append((Num[1]*x[-2] + Num[0]*x[-1] \
						- Den[1]*TemporaryIbInput[-2])/Den[0])
				elif i >= 3*SamplingRatio:
					TemporaryIbInput.append((Num[2]*x[-3] + Num[1]*x[-2] + Num[0]*x[-1] \
						- Den[2]*TemporaryIbInput[-3] - Den[1]*TemporaryIbInput[-2])/Den[0])

				IbInput.append(TemporaryIbInput[-1]*(TemporaryIbInput[-1]>0))
			
				if i in range(IaAfferentDelayTimeStep-1,IbAfferentDelayTimeStep):
					Input.append(IaInput[i-IaAfferentDelayTimeStep+1]/IaSpindleGain \
								+FeedforwardInput[i]) #input to the muscle
				elif i in range(IbAfferentDelayTimeStep, IIAfferentDelayTimeStep):
					Input.append(IaInput[i-IaAfferentDelayTimeStep+1]/IaSpindleGain \
								-IbInput[i-IbAfferentDelayTimeStep]/IbGTOGain \
								+FeedforwardInput[i]) #input to the muscle
				elif i in range(IIAfferentDelayTimeStep, CorticalDelayTimeStep):
					Input.append(IaInput[i-IaAfferentDelayTimeStep+1]/IaSpindleGain \
								+IIInput[i-IIAfferentDelayTimeStep]/IISpindleGain \
								-IbInput[i-IbAfferentDelayTimeStep]/IbGTOGain \
								+FeedforwardInput[i]) #input to the muscle
				else:	
					Input.append(FeedforwardInput[i])
				
				Count = SamplingRatio+Count
			else:
				IbInput.append(IbInput[-1])
				Input.append(Input[-1])
		elif FeedbackOption == 'fb_control': # Feedback control (proprioceptive systems + supraspinal loop)
			PrimaryOutput,SecondaryOutput = spindle_model(CE,Bag1,Bag2,Chain)
			IaInput.append(PrimaryOutput)
			IIInput.append(SecondaryOutput)
			if i == Count:
				# Get spindle primary and secondary afferent activity
				
				# Get Ib activity
				x.append(GTOConstant1*np.log((SEE['Force'][-1]/GTOConstant2+1)))
				if i == 0:
					TemporaryIbInput.append(x[-1])
				elif i == 1*SamplingRatio:
					TemporaryIbInput.append(x[-1])
				elif i == 2*SamplingRatio:
					TemporaryIbInput.append((Num[0]*x[-1])/Den[0])
				elif i == 3*SamplingRatio:
					TemporaryIbInput.append((Num[1]*x[-2] + Num[0]*x[-1]- Den[1]*TemporaryIbInput[-1])/Den[0])
				elif i > 3*SamplingRatio:
					TemporaryIbInput.append((Num[2]*x[-3] + Num[1]*x[-2] + Num[0]*x[-1] \
														- Den[2]*TemporaryIbInput[-2] - Den[1]*TemporaryIbInput[-1])/Den[0])

				IbInput.append(TemporaryIbInput[-1]*(TemporaryIbInput[-1]>0))

				if i in range(IaAfferentDelayTimeStep-1,IbAfferentDelayTimeStep):
					Input.append(IaInput[i-IaAfferentDelayTimeStep+1]/IaSpindleGain\
								+ FeedforwardInput[i]	) #input to the muscle
				elif i in range(IbAfferentDelayTimeStep, IIAfferentDelayTimeStep):
					Input.append(IaInput[i-IaAfferentDelayTimeStep+1]/IaSpindleGain \
								- IbInput[i-IbAfferentDelayTimeStep]/IbGTOGain \
								+ FeedforwardInput[i] ) #input to the muscle
				elif i in range(IIAfferentDelayTimeStep, CorticalDelayTimeStep):
					Input.append(IaInput[i-IaAfferentDelayTimeStep+1]/IaSpindleGain \
								+ IIInput[i-IIAfferentDelayTimeStep]/IISpindleGain \
								- IbInput[i-IbAfferentDelayTimeStep]/IbGTOGain \
								+ FeedforwardInput[i]	) #input to the muscle
				elif i >= CorticalDelayTimeStep:
					FeedbackInput = TransCorticalLoopConstant*(TargetForceTrajectory[i]-SEE['Force'][i-CorticalDelayTimeStep+1])/MaximumContractileElementForce + FeedbackInput  # feedback input through cortical pathway
					Input.append(IaInput[i-IaAfferentDelayTimeStep+1]/IaSpindleGain \
								+IIInput[i-IIAfferentDelayTimeStep]/IISpindleGain \
								-IbInput[i-IbAfferentDelayTimeStep]/IbGTOGain \
								+FeedbackInput	) #input to the muscle
				else:	
					Input.append(FeedforwardInput[i])
				
				Count = SamplingRatio+Count
			else:
				IbInput.append(IbInput[-1])
				Input.append(Input[-1])
		elif FeedbackOption == 'cortical_fb_only':
			if i > CorticalDelayTimeStep-1:
				FeedbackInput = TransCorticalLoopConstant*(TargetForceTrajectory[i]-OutputForceTendon[i-CorticalDelayTimeStep-1])/MaximumContractileElementForce + FeedbackInput  # feedback input through cortical pathway
				Input.append(FeedbackInput) #input to the muscle
			else:
				Input.append(FeedforwardInput[i])

			IaInput.append(0)
			IIInput.append(0)
			IbInput.append(0)

		# Feedforward only
		# Input[i] = input[i] #TargetForceTrajectory[i]/MaximumContractileElementForce
		## Noise + additional input
		if Input[-1] < 0:
			Input[-1] = 0
		elif Input[-1] > 1:
			Input[-1] = 1

		random.seed(1)
		if i > 4:
			Noise.append(2*(random.random()-0.5)*(np.sqrt(0.01*Input[i])*np.sqrt(3)))
			FilteredNoise.append((BButtersCoefficients[4]*Noise[i-4] + BButtersCoefficients[3]*Noise[i-3] + BButtersCoefficients[2]*Noise[i-2] + BButtersCoefficients[1]*Noise[i-1] + BButtersCoefficients[0]*Noise[i] \
									- AButtersCoefficients[4]*FilteredNoise[i-4] - AButtersCoefficients[3]*FilteredNoise[i-3] - AButtersCoefficients[2]*FilteredNoise[i-2] - AButtersCoefficients[1]*FilteredNoise[i-1])/AButtersCoefficients[0])
		else:
			Noise.append(0)
			FilteredNoise.append(0)

		Input[-1] = Input[-1] + FilteredNoise[-1] + CorticalInput[-1]

		if Input[-1] < 0:
			Input[-1] = 0
		elif Input[-1] > 1:
			Input[-1] = 1

		# add delay along efferent pathway
		if i > EfferentDelayTimeStep:
			LongInput.append(Input[i-EfferentDelayTimeStep])
		else:
			LongInput.append(0)

		# add activation filter
		TU = 0.12*(LongInput[-1] < EffectiveMuscleActivation) + 0.03 # When true TU = 0.15, else TU = 0.03

		EffectiveMuscleActivationFirstDeriv = (LongInput[-1] - EffectiveMuscleActivation)/TU
		EffectiveMuscleActivation = EffectiveMuscleActivationFirstDeriv*(SamplingPeriod) + EffectiveMuscleActivation # effective neural drive

		# SlowTwitch = activation_frequency_slow(CE,SlowTwitch,EffectiveMuscleActivation,SamplingFrequency) # not used

		# force-velocity relationship
		ForceVelocity = (CE['Velocity'][-1] <= 0)*concentric_force_velocity(CE) \
							+ (CE['Velocity'][-1] > 0)*eccentric_force_velocity(CE) 
		# force-length relationship
		ForceLength = force_length(CE)
		ContractileElementForce = ForceLength*ForceVelocity
		# viscous property
		ForceViscocity = MuscleViscosity * CE['Velocity'][-1]
		# passive element 1
		ForcePassive1 = parallel_elastic_element_force_1(CE)
		# passive element 2
		ForcePassive2 = parallel_elastic_element_force_2(CE)
		ForcePassive2 = (ForcePassive2 <= 0)*ForcePassive2

		# total force from contractile element
		ForceTotal = (EffectiveMuscleActivation*(ContractileElementForce + ForcePassive2) + ForcePassive1 + ForceViscocity)*MaximumContractileElementForce
		ForceTotal = ForceTotal*(ForceTotal>=0.0)

		# force from series elastic element
		SEE['Force'].append(normalized_series_elastic_element_force(SEE) * MaximumContractileElementForce)

		# calculate muscle excursion acceleration based on the difference
		# between muscle force and tendon force
		if i < len(Time)-1:
			MuscleAcceleration.append((SEE['Force'][-1]*np.cos(PennationAngle) - ForceTotal*(np.cos(PennationAngle))**2)/(MuscleMass) \
				+ (MuscleVelocity[-1])**2*np.tan(PennationAngle)**2/(MuscleLength[-1]))
			# integrate acceleration to get velocity
			MuscleVelocity.append((MuscleAcceleration[-1]+ \
				MuscleAcceleration[-2])/2*(SamplingPeriod)+MuscleVelocity[-1])
			# integrate velocity to get length
			MuscleLength.append((MuscleVelocity[-1]+ \
				MuscleVelocity[-2])/2*(SamplingPeriod)+MuscleLength[-1])

			# normalize each variable to optimal muscle length or tendon legnth
			CE['Acceleration'].append(MuscleAcceleration[-1]/(OptimalLength/100))
			CE['Velocity'].append(MuscleVelocity[-1]/(OptimalLength/100))
			CE['Length'].append(MuscleLength[-1]/(OptimalLength/100))
			SEE['Length'].append((InitialMusculoTendonLength - CE['Length'][-1]*OptimalLength*np.cos(PennationAngle))/OptimalTendonLength)

			# store data
			OutputActivationFrequency.append(ActivationFrequency)
			OutputEffectiveMuscleActivation.append(EffectiveMuscleActivation)

		OutputForceMuscle.append(ForceTotal)
		OutputForceLength.append(ForceLength)
		OutputForceVelocity.append(ForceVelocity)
		OutputForcePassive1.append(ForcePassive1)
		OutputForcePassive2.append(ForcePassive2)
		statusbar(i,len(Time),StartTime=StartTime,Title = 'afferented_muscle_model')

	# plt.figure()
	# plt.plot(Time,SEE['Force'][1:])
	# plt.plot(Time,TargetForceTrajectory,'r')
	# plt.legend(['Output Force','Target Force'])
	# plt.xlabel('Time (sec)')
	# plt.ylabel('Force (N)')

	# save data as output in structure format
	output = {	'Force' : OutputForceMuscle, 										'ForceTendon' : SEE['Force'][1:], \
				'FL' : OutputForceLength, 											'FV' : OutputForceVelocity, \
				'PE Force 1' : OutputForcePassive1, 								'PE Force 2' : OutputForcePassive2, \
				'Target' : TargetForceTrajectory, 									'ContractileElementLength' : CE['Length'], \
				'ContractileElementVelocity' : CE['Velocity'][1:],				 	'ContractileElementAcceleration' : CE['Acceleration'][1:], \
				'SeriesElasticElementLength' : SEE['Length'],					 	'Activation Frequency' : OutputActivationFrequency, \
				'Input' : Input, 													'Noise' : Noise, \
				'FilteredNoise' : FilteredNoise, 									'U' : OutputEffectiveMuscleActivation, \
				'Ia Input' : IaInput,												'II Input' : IIInput, \
				'Ib Input' : IbInput 	} #, 'CorticalInput' : CorticalInput }
	"""
	NOTES:

	'ContractileElementVelocity','ContractileElementAcceleration', and 'ForceTendon' are from 1: because originally the lists were not initialized even though Vce and Ace are initialized to zero


	"""
	return(output)

muscle_parameters = {	"Pennation Angle":5*np.pi/180, 		"Muscle Mass":0.075,\
						"Optimal Length":10.1, 				"Tendon Length":23.5,\
					 	"Initial Muscle Length":10.1+0.41, 	"Initial Tendon Length":23.5+0.09}

# Define delays based on limb length and conduction velocity of each pathway
DistanceToTheSpinalCord = 0.8 # cm
EfferentConductionVelocity = 48.5 #m/s (ref. Elias et al. 2014) S-type = 44-51, FR = 51-52, FF 52-53
IaConductionVelocity = 64.5 #m/s (ref. Elias et al. 2014)
IIConductionVelocity = 32.5 #m/s (ref. Elias et al. 2014)
IbConductionVelocity = 59 #m/s (ref. Elias et al. 2014)
SynapticDelay = 2 #ms (ref. Kandel)

delay_parameters = {"Efferent Delay": round(DistanceToTheSpinalCord/EfferentConductionVelocity*1000), \
					"Ia Delay" : round(DistanceToTheSpinalCord/IaConductionVelocity*1000) + SynapticDelay, \
					"II Delay" : round(DistanceToTheSpinalCord/IIConductionVelocity*1000) + 2*SynapticDelay, \
					"Ib Delay" : round(DistanceToTheSpinalCord/IbConductionVelocity*1000) + 2*SynapticDelay, \
					"Cortical Delay" : 50 }

# Define gain parameters for each neural pathway
gain_parameters = {	"Gamma Dynamic Gain" : 70, \
					"Gamma Static Gain" : 70, \
					"Ia Gain" : 800, \
					"II Gain" : 3000, \
					"Ib Gain" : 1000 }

# Define force trajectory that model needs to track
# You can create any target trajectory as you want, but it has to be
# sampled at 10000 Hz to be consistent with the sampling frequency of
# muscle model
SamplingFrequency = 10000
Time = np.arange(0,15,1/SamplingFrequency) #0:SamplingPeriod:15
TrajectoryType = 'constant' #'triangle','sinwave','trapezoid' (Still need to do Trap)
TargetTrajectory = generate_target_force_trajectory(TrajectoryType,Time,0.3,0,0) # in unit of #MVC

# Define additional input to muscle (e.g., 20 Hz oscillatory input)
# CorticalInputAmplitude = (1/10)*TargetTrajectory # in unit of #MVC
# CorticalInputFrequency = 20 # in unit of Hz
# CorticalInput = CorticalInputAmplitude*np.sin(2*np.pi*CorticalInputFrequency*Time)
CorticalInput = np.zeros(len(TargetTrajectory))

FeedbackOption = 'fb_control'
if FeedbackOption == 'ff_only':
    FeedbackOptionString = 'FF'
    TargetTrajectory = 1.028*TargetTrajectory
elif FeedbackOption == 'servo_control':
    FeedbackOptionString = 'SC'
    TargetTrajectory = 0.9285*TargetTrajectory
elif FeedbackOption == 'fb_control':
    FeedbackOptionString = 'FB'
elif FeedbackOption == 'cortical_fb_only':
    FeedbackOptionString = 'CC'

# Output from the muscle model
NumberOfTrials = 1
Output = []
PXX = []
Range = []
for i in range(NumberOfTrials): #trialN = 1:10
	#rng('shufle')
	output = afferented_muscle_model(muscle_parameters,delay_parameters,gain_parameters,TargetTrajectory,CorticalInput,FeedbackOption = FeedbackOption)    
	compare_output(output)
	Output.append(output)
	f,pxx = welch(output['ForceTendon'][-int(10*SamplingFrequency)-1:]-np.average(output['ForceTendon'][int(-10*SamplingFrequency-1):]),window = gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5)),\
		noverlap = SamplingFrequency,nperseg = len(gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5))), fs = SamplingFrequency) 
	pxx = smooth(pxx,5)
	PXX.append(np.array(pxx,ndmin=2))
	Range.append(max(output['Ia Input'][-10*SamplingFrequency:-1])-min(output['Ia Input'][-10*SamplingFrequency-1:]))
	print("\n")

PXX = np.concatenate(PXX,axis=0)
# plt.figure()
# plt.plot(f[:131],np.average(PXX[:,:131],axis = 0))

#max(output.Ia(end-10*SamplingFrequency:end))-min(output.Ia(end-10*SamplingFrequency:end))
#Input to the MN model
# Input = output.Input
#
# pool_parameter.N = 120
# pool_parameter.gain = 1
# pool_parameter.ISI_CoV = 0.2
# forceOption = 0
#
# [time_MN,spike_train] = MotorUnitModel_071316(pool_parameter,Input,forceOption)
#
# spike_train_logical = logical(spike_train)
# plotRaster(spike_train_logical(1:120,:),time_MN)