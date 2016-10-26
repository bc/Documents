import numpy as np
from scipy.signal import sawtooth,square,gaussian,welch
import matplotlib.pyplot as plt
import time
import ipdb
from amm_functions import *

def afferented_2_muscles_model(muscle_1_parameters,muscle_2_parameters,\
							delay_1_parameters,delay_2_parameters,\
							gain_1_parameters,gain_2_parameters,
							TargetTrajectory1,TargetTrajectory2,\
							CorticalInput1,CorticalInput2,\
							**kwargs):
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
	test_input_values(muscle_1_parameters,delay_1_parameters,gain_1_parameters,FeedbackOption)
	test_input_values(muscle_2_parameters,delay_2_parameters,gain_2_parameters,FeedbackOption)

	import numpy as np 
	from scipy import signal
	import control
	import random
	import matplotlib.pyplot as plt

	# muscle architectural parameters
	PennationAngle1 = muscle_1_parameters['Pennation Angle']
	MuscleMass1 = muscle_1_parameters['Muscle Mass']
	OptimalLength1 = muscle_1_parameters['Optimal Length']
	TendonLength1 = muscle_1_parameters['Tendon Length']
	OptimalTendonLength1 = TendonLength1*1.05
	InitialMuscleLength1 = muscle_1_parameters['Initial Muscle Length']
	InitialTendonLength1 = muscle_1_parameters['Initial Tendon Length']
	InitialMusculoTendonLength1 = InitialMuscleLength1*np.cos(PennationAngle1)+InitialTendonLength1

	PennationAngle2 = muscle_2_parameters['Pennation Angle']
	MuscleMass2 = muscle_2_parameters['Muscle Mass']
	OptimalLength2 = muscle_2_parameters['Optimal Length']
	TendonLength2 = muscle_2_parameters['Tendon Length']
	OptimalTendonLength2 = TendonLength2*1.05
	InitialMuscleLength2 = muscle_2_parameters['Initial Muscle Length']
	InitialTendonLength2 = muscle_2_parameters['Initial Tendon Length']
	InitialMusculoTendonLength2 = InitialMuscleLength2*np.cos(PennationAngle2)+InitialTendonLength2

	# MaximumContractileElementLength1, \
	# 		ContractileElementLength1,SeriesElasticElementLength1, \
	# 		MaximumContractileElementForce1,ContractileElementVelocity1,ContractileElementAcceleration1, \
	# 		MuscleAcceleration1,MuscleVelocity1,MuscleLength1, \
	# 		SeriesElasticElementForce1,Y_dot1,Y1, \
	# 		Saf_dot1,Saf1,fint_dot1, \
	# 		fint1,feff_dot1,feff1, \
	# 		ActivationFrequency1,EffectiveMuscleActivation1,DynamicSpindleFrequency1, \
	# 		StaticSpindleFrequency1,Bag1TensionSecondDeriv1,Bag1TensionFirstDeriv1, \
	# 		Bag1Tension1,Bag2TensionSecondDeriv1,Bag2TensionFirstDeriv1, \
	# 		Bag2Tension1,ChainTensionSecondDeriv1,ChainTensionFirstDeriv1, \
	# 		ChainTension1
	# 			= return_initial_values(PennationAngle1,MuscleMass1,OptimalLength1,OptimalTendonLength1,\
	# 									TendonLength1,InitialMuscleLength1,InitialTendonLength1,InitialMusculoTendonLength1)
	# MaximumContractileElementLength2, \
	# 		ContractileElementLength2,SeriesElasticElementLength2, \
	# 		MaximumContractileElementForce2,ContractileElementVelocity2,ContractileElementAcceleration2, \
	# 		MuscleAcceleration2,MuscleVelocity2,MuscleLength2, \
	# 		SeriesElasticElementForce2,Y_dot2,Y2, \
	# 		Saf_dot2,Saf2,fint_dot2, \
	# 		fint2,feff_dot2,feff2, \
	# 		ActivationFrequency2,EffectiveMuscleActivation2,DynamicSpindleFrequency2, \
	# 		StaticSpindleFrequency2,Bag1TensionSecondDeriv2,Bag1TensionFirstDeriv2, \
	# 		Bag1Tension2,Bag2TensionSecondDeriv2,Bag2TensionFirstDeriv2, \
	# 		Bag2Tension2,ChainTensionSecondDeriv2,ChainTensionFirstDeriv2, \
	# 		ChainTension2
	# 			= return_initial_values(PennationAngle2,MuscleMass2,OptimalLength2,OptimalTendonLength2,\
	# 									TendonLength2,InitialMuscleLength2,InitialTendonLength2,InitialMusculoTendonLength2)

	# visnp.cosity of muscle (ref: Elias et al. 2014)
	MuscleViscosity = 0.005 #0.001

	# assign delays
	# EfferentDelay1 = delay_1_parameters['Efferent Delay']
	# IaAfferentDelay1 = delay_1_parameters['Ia Delay']
	# IIAfferentDelay1 = delay_1_parameters['II Delay']
	# IbAfferentDelay1 = delay_1_parameters['Ib Delay']
	# CorticalDelay1 = delay_1_parameters['Cortical Delay']

	# EfferentDelay2 = delay_2_parameters['Efferent Delay']
	# IaAfferentDelay2 = delay_2_parameters['Ia Delay']
	# IIAfferentDelay2 = delay_2_parameters['II Delay']
	# IbAfferentDelay2 = delay_2_parameters['Ib Delay']
	# CorticalDelay2 = delay_2_parameters['Cortical Delay']

	SamplingFrequency = 10000
	FeedbackSamplingFrequency = 1000
	SamplingRatio = SamplingFrequency/FeedbackSamplingFrequency
	SamplingPeriod = 1/SamplingFrequency
	Time = np.arange(0,len(TargetTrajectory1)*SamplingPeriod,SamplingPeriod) # 0:SamplingPeriod:(length(TargetTrajectory)-1)/SamplingFrequency

	# filter parameters for Noise
	BButtersCoefficients,AButtersCoefficients = signal.butter(4,100/(SamplingFrequency/2),'low')

	# compare(BButtersCoefficients,np.array([0.0898e-05, 0.3594e-05, 0.5391e-05, 0.3594e-05, 0.0898e-05]))
	# compare(AButtersCoefficients,np.array([1.0000, -3.8358,  5.5208, -3.5335,  0.8486]))

	Num,Den = [1.7,2.58,0.4],[1,2.2,0.4]
	ContinuousTransferFunction = control.tf(Num,Den)
	DiscreteTransferFunction = control.matlab.c2d(ContinuousTransferFunction,1/FeedbackSamplingFrequency)
	Num,Den = control.matlab.tfdata(DiscreteTransferFunction)
	Num,Den = Num[0][0],Den[0][0]

	# compare(Num,np.array([1.7000,-3.3974, 1.6974]))
	# compare(Den,np.array([1.0000,-1.9978,0.9978]))

	## Transcortical loop
	TransCorticalLoopConstant = 0.01
	
	## Gains
	# GammaDynamicGain1 = gain_1_parameters['Gamma Dynamic Gain']
	# GammaStaticGain1 = gain_1_parameters['Gamma Static Gain']
	# IaSpindleGain1 = gain_1_parameters['Ia Gain']
	# IISpindleGain1 = gain_1_parameters['II Gain']
	# IbGTOGain1 = gain_1_parameters['Ib Gain']

	# GammaDynamicGain2 = gain_2_parameters['Gamma Dynamic Gain']
	# GammaStaticGain2 = gain_2_parameters['Gamma Static Gain']
	# IaSpindleGain2 = gain_2_parameters['Ia Gain']
	# IISpindleGain2 = gain_2_parameters['II Gain']
	# IbGTOGain2 = gain_2_parameters['Ib Gain']

	# # Convert1 delays in ms to samples
	# IaAfferentDelayTimeStep1 = int(IaAfferentDelay1*SamplingRatio)   #Ia +1 II 30
	# IIAfferentDelayTimeStep1 = int(IIAfferentDelay1*SamplingRatio)
	# IbAfferentDelayTimeStep1 = int(IbAfferentDelay1*SamplingRatio)   #Ib 40
	# CorticalDelayTimeStep1 = int(CorticalDelay1*SamplingRatio)   #cortical 50
	# EfferentDelayTimeStep1 = int(EfferentDelay1*SamplingRatio)

	# IaAfferentDelayTimeStep2 = int(IaAfferentDelay2*SamplingRatio)   #Ia + II 30
	# IIAfferentDelayTimeStep2 = int(IIAfferentDelay2*SamplingRatio)
	# IbAfferentDelayTimeStep2 = int(IbAfferentDelay2*SamplingRatio)   #Ib 40
	# CorticalDelayTimeStep2 = int(CorticalDelay2*SamplingRatio)   #cortical 50
	# EfferentDelayTimeStep2 = int(EfferentDelay2*SamplingRatio)

	# FeedbackInput1, FeedbackInput2 = 0,0

	# Convert force trajectory to unit of newton
	# TargetForceTrajectory1 = TargetTrajectory1*CEE_1['MaximumForce']
	# TargetForceTrajectory2 = TargetTrajectory2*CEE_2['MaximumForce']
	# FeedforwardInput1 = TargetForceTrajectory1/CEE_1['MaximumForce']
	# FeedforwardInput2 = TargetForceTrajectory2/CEE_2['MaximumForce']
	def contractile_element_parameters(Length,Velocity,Acceleration,MaximumContractileElementLength,ForceLength,ForceVelocity,MaximumContractileElementForce):
		CE = initialize_dictionary(['Length','Velocity','Acceleration','MaximumLength','ForceLength','ForceVelocity','MaximumForce'],\
									[[Length],[Velocity],[Acceleration],MaximumContractileElementLength,ForceLength,ForceVelocity,MaximumContractileElementForce])
		return(CE)
	def bag_1_parameters(GammaDynamicGain,DynamicSpindleFrequency,Bag1Tension,Bag1TensionFirstDeriv):
		Bag1 = initialize_dictionary(['GammaDynamicGain','DynamicSpindleFrequency','Tension','TensionFirstDeriv'],\
										[GammaDynamicGain,DynamicSpindleFrequency,Bag1Tension,Bag1TensionFirstDeriv])
		return(Bag1)
	def bag_2_parameters(GammaStaticGain,StaticSpindleFrequency,Bag2Tension,Bag2TensionFirstDeriv):
			Bag2 = initialize_dictionary(['GammaStaticGain','StaticSpindleFrequency','Tension','TensionFirstDeriv'],\
											[GammaStaticGain,StaticSpindleFrequency,Bag2Tension,Bag2TensionFirstDeriv])
			return(Bag2)
	def chain_parameters(ChainTension,ChainTensionFirstDeriv):
		Chain = initialize_dictionary(['Tension','TensionFirstDeriv'],\
										[ChainTension, ChainTensionFirstDeriv])
		return(Chain)
	def slow_twitch_parameters(Y,fint,feff_dot,feff,ActivationFrequencySlow):
		SlowTwitch = initialize_dictionary(['Y','fint','feff_dot','feff','ActivationFrequency'],\
											[Y,fint,feff_dot,feff,ActivationFrequencySlow])
		return(SlowTwitch)
	def fast_twitch_parameters(Saf,fint,feff_dot,feff,ActivationFrequencyFast):
		FastTwitch = initialize_dictionary(['Saf','fint','feff_dot','feff','ActivationFrequency'],\
											[Saf,fint,feff_dot,feff,ActivationFrequencyFast])
		return(FastTwitch)
	def series_elastic_element_parameters(Length,Force):
		SEE = initialize_dictionary(['Length','Force', 'Temp'],\
									[[Length], [Force], []])
		return(SEE)
	def parallel_elastic_element_parameters(ForcePassive1,ForcePassive2):
		PEE = initialize_dictionary(['ForcePassive1','ForcePassive2'],\
									[ForcePassive1,ForcePassive2])
		return(PEE)
	# CE_1 = contractile_element_parameters(ContractileElementLength1,ContractileElementVelocity1,ContractileElementAcceleration1,\
	# 										MaximumContractileElementLength1,[],[],MaximumContractileElementForce1)
	# CE_2 = contractile_element_parameters(ContractileElementLength2,ContractileElementVelocity2,ContractileElementAcceleration2,\
	# 										MaximumContractileElementLength2,[],[],MaximumContractileElementForce2)
	# Bag1_1 = bag_1_parameters(GammaDynamicGain1,DynamicSpindleFrequency1,Bag1Tension1,Bag1TensionFirstDeriv1)
	# Bag1_2 = bag_1_parameters(GammaDynamicGain2,DynamicSpindleFrequency2,Bag1Tension2,Bag1TensionFirstDeriv2)
	# Bag2_1 = bag_2_parameters(GammaStaticGain1,StaticSpindleFrequency1,Bag2Tension1,Bag2TensionFirstDeriv1)
	# Bag2_2 = bag_2_parameters(GammaStaticGain2,StaticSpindleFrequency2,Bag2Tension2,Bag2TensionFirstDeriv2)
	# Chain_1 = chain_parameters(ChainTension1,ChainTensionFirstDeriv1)
	# Chain_2 = chain_parameters(ChainTension2,ChainTensionFirstDeriv2)
	# SlowTwitch_1 = slow_twitch_parameters(Y1,fint1,feff_dot1,feff1,ActivationFrequency1)
	# SlowTwitch_2 = slow_twitch_parameters(Y2,fint2,feff_dot2,feff2,ActivationFrequency2)
	# FastTwitch_1 = slow_twitch_parameters(Saf1,fint1,feff_dot1,feff1,ActivationFrequency1)
	# FastTwitch_2 = slow_twitch_parameters(Saf2,fint2,feff_dot2,feff2,ActivationFrequency2)
	# SEE_1 = series_elastic_element_parameters(SeriesElasticElementLength1,SeriesElasticElementForce1)
	# SEE_2 = series_elastic_element_parameters(SeriesElasticElementLength2,SeriesElasticElementForce2)
	# PEE_1 = parallel_elastic_element_parameters([],[])
	# PEE_2 = parallel_elastic_element_parameters([],[])
	# Initial Feedback Values
	FeedbackInput1, FeedbackInput2 = 0,0

	# Set Initial Values and Convert force trajectory to unit of newton
	Bag1_1,Bag2_1,Chain_1,SlowTwitch_1,FastTwitch_1,CE_1,SEE_1,PEE_1,Muscle_1,Input_1 \
		= return_initial_values(muscle_1_parameters,gain_1_parameters)
	TargetForceTrajectory1 = TargetTrajectory1*CE_1['MaximumForce']
	FeedforwardInput1 = TargetForceTrajectory1/CE_1['MaximumForce']
	Input_1['Feedforward'],Input_1['TargetTrajectory']=FeedforwardInput1,TargetForceTrajectory1
	Input_1['CorticalInput'],Input_1['Feedback'] = CorticalInput1,[FeedbackInput1]

	Bag1_2,Bag2_2,Chain_2,SlowTwitch_2,FastTwitch_2,CE_2,SEE_2,PEE_2,Muscle_2,Input_2 \
		= return_initial_values(muscle_2_parameters,gain_2_parameters)
	TargetForceTrajectory2 = TargetTrajectory2*CE_2['MaximumForce']
	FeedforwardInput2 = TargetForceTrajectory2/CE_2['MaximumForce']
	Input_2['Feedforward'],Input_2['TargetTrajectory']=FeedforwardInput2,TargetForceTrajectory2
	Input_2['CorticalInput'],Input_2['Feedback'] = CorticalInput2,[FeedbackInput2]
	
	def bag1_model(CE,Bag1):
		## Feedback system parameters
		Length, LengthFirstDeriv, LengthSecondDeriv = CE['Length'][-1], CE['Velocity'][-1], CE['Acceleration'][-1]
		GammaDynamicGain, DynamicSpindleFrequency = Bag1['GammaDynamicGain'], Bag1['DynamicSpindleFrequency']
		Bag1Tension, Bag1TensionFirstDeriv = Bag1['Tension'], Bag1['TensionFirstDeriv']

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

		Bag1TensionSecondDeriv = (K_SR/M)*(	(C*Bag1Beta*np.sign(LengthFirstDeriv-Bag1TensionFirstDeriv/K_SR) \
											* ((abs(LengthFirstDeriv-Bag1TensionFirstDeriv/K_SR))**a)  \
											* (Length-L0_SR-Bag1Tension/K_SR-R)) \
											+ K_PR*(Length-L0_SR-Bag1Tension/K_SR-L0_PR) \
											+ M*LengthSecondDeriv \
											+ Bag1Gamma \
											- Bag1Tension 	)
		Bag1TensionFirstDeriv = Bag1TensionSecondDeriv*(SamplingPeriod) + Bag1TensionFirstDeriv
		Bag1Tension = Bag1TensionFirstDeriv*(SamplingPeriod) + Bag1Tension

		Bag1['AfferentPotential'] = G*(Bag1Tension/K_SR-(LN_SR-L0_SR))
		Bag1['DynamicSpindleFrequency']=DynamicSpindleFrequency
		Bag1['Tension']=Bag1Tension
		Bag1['TensionFirstDeriv'] = Bag1TensionFirstDeriv
		return(Bag1)
	def bag2_model(CE,Bag2):
		## Feedback system parameters
		Length, LengthFirstDeriv, LengthSecondDeriv = CE['Length'][-1], CE['Velocity'][-1], CE['Acceleration'][-1]
		GammaStaticGain, StaticSpindleFrequency = Bag2['GammaStaticGain'], Bag2['StaticSpindleFrequency']
		Bag2Tension, Bag2TensionFirstDeriv = Bag2['Tension'], Bag2['TensionFirstDeriv']
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

		#nan_test(StaticSpindleFrequencyFirstDeriv,'StaticSpindleFrequencyFirstDeriv')
		#nan_test(StaticSpindleFrequency,'StaticSpindleFrequency')

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
		Bag2TensionFirstDeriv = Bag2TensionSecondDeriv*(SamplingPeriod) + Bag2TensionFirstDeriv
		Bag2Tension = Bag2TensionFirstDeriv*(SamplingPeriod) + Bag2Tension

		Bag2['PrimaryAfferentPotential'] = G*(Bag2Tension/K_SR-(LN_SR-L0_SR))
		Bag2['SecondaryAfferentPotential'] = G*(	X*(L_secondary/L0_SR)*(Bag2Tension/K_SR-(LN_SR-L0_SR)) \
							+(1-X)*(L_secondary/L0_PR)*(Length-Bag2Tension/K_SR-(L0_SR+LN_PR))	)
		Bag2['StaticSpindleFrequency'] = StaticSpindleFrequency
		Bag2['Tension']=Bag2Tension
		Bag2['TensionFirstDeriv']=Bag2TensionFirstDeriv
		return(Bag2)
	def chain_model(CE,Bag2,Chain):
		## Feedback system parameters
		Length, LengthFirstDeriv, LengthSecondDeriv = CE['Length'][-1], CE['Velocity'][-1], CE['Acceleration'][-1]
		GammaStaticGain, StaticSpindleFrequency = Bag2['GammaStaticGain'], Bag2['StaticSpindleFrequency']
		ChainTension, ChainTensionFirstDeriv = Chain['Tension'], Chain['TensionFirstDeriv']
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
		ChainTensionFirstDeriv = ChainTensionSecondDeriv*SamplingPeriod + ChainTensionFirstDeriv
		ChainTension = ChainTensionFirstDeriv*SamplingPeriod + ChainTension
		Chain['PrimaryAfferentPotential'] = G*(ChainTension/K_SR-(LN_SR-L0_SR))
		Chain['SecondaryAfferentPotential'] = G*(	X*(L_secondary/L0_SR)*(ChainTension/K_SR-(LN_SR-L0_SR)) \
							+ (1-X)*(L_secondary/L0_PR)*(Length-ChainTension/K_SR-(L0_SR+LN_PR))	)
		Chain['Tension'] = ChainTension
		Chain['TensionFirstDeriv'] = ChainTensionFirstDeriv
		return(Chain)
	def spindle_model(CE,Bag1,Bag2,Chain):
		S = 0.156

		Bag1 = bag1_model(CE,Bag1)
		Bag2 = bag2_model(CE,Bag2)
		Chain = chain_model(CE,Bag2,Chain)

		if Bag1['AfferentPotential'] < 0: Bag1['AfferentPotential'] = 0 
		if Bag2['PrimaryAfferentPotential'] < 0: Bag2['PrimaryAfferentPotential'] = 0
		if Chain['PrimaryAfferentPotential'] < 0: Chain['PrimaryAfferentPotential'] = 0
		if Bag2['SecondaryAfferentPotential'] < 0: Bag2['SecondaryAfferentPotential'] = 0
		if Chain['SecondaryAfferentPotential'] < 0: Chain['SecondaryAfferentPotential'] = 0	        

		if Bag1['AfferentPotential'] > (Bag2['PrimaryAfferentPotential']+Chain['PrimaryAfferentPotential']):
			Larger = Bag1['AfferentPotential']
			Smaller = Bag2['PrimaryAfferentPotential']+Chain['PrimaryAfferentPotential']
		else:
			Larger = Bag2['PrimaryAfferentPotential']+Chain['PrimaryAfferentPotential']
			Smaller = Bag1['AfferentPotential']

		PrimaryOutput = Larger + S * Smaller
		SecondaryOutput = Bag2['SecondaryAfferentPotential'] + Chain['SecondaryAfferentPotential']

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
		CE['ForceLength'].append(np.exp(-abs((Length**beta - 1)/omega)**rho))
	def concentric_force_velocity(CE):
		Length,Velocity = CE['Length'][-1],CE['Velocity'][-1]
		MaximumVelocity = -7.88
		cv0 = 5.88
		cv1 = 0
		ConcentricForceVelocity = (MaximumVelocity - Velocity)/(MaximumVelocity + (cv0 + cv1*Length)*Velocity)
		#nan_test(ConcentricForceVelocity,'ConcentricForceVelocity')
		return(ConcentricForceVelocity)
	def eccentric_force_velocity(CE):
		Length,Velocity = CE['Length'][-1],CE['Velocity'][-1]
		av0 = -4.7
		av1 = 8.41
		av2 = -5.34
		bv = 0.35
		EccentricForceVelocity = (bv - (av0 + av1*Length + av2*Length**2)*Velocity)/(bv+Velocity)
		#nan_test(EccentricForceVelocity,'EccentricForceVelocity')
		return(EccentricForceVelocity)
	def parallel_elastic_element_force_1(PEE,CE):
		Length,MaximumContractileElementLength = CE['Length'][-1],CE['MaximumLength']
		c1_pe1 = 23.0  #355, 67.1
		k1_pe1 = 0.046 #0.04, 0.056
		Lr1_pe1 = 1.17  #1.35, 1.41
		ParallelElasticElementForce1 = c1_pe1 * k1_pe1 * np.log(np.exp((Length/MaximumContractileElementLength - Lr1_pe1)/k1_pe1)+1)
		PEE['ForcePassive1'].append(ParallelElasticElementForce1)
	def parallel_elastic_element_force_2(PEE,CE):
		Length = CE['Length'][-1]
		c2_pe2 = -0.02 #0.01  -0.1
		k2_pe2 = -21
		Lr2_pe2 = 0.70 #0.79 0.59
		ParallelElasticElementForce2 = c2_pe2*np.exp((k2_pe2*(Length-Lr2_pe2))-1)
		ParallelElasticElementForce2=(ParallelElasticElementForce2<=0)*ParallelElasticElementForce2
		PEE['ForcePassive2'].append(ParallelElasticElementForce2)
	def normalized_series_elastic_element_force(SEE):
		LT = SEE['Length'][-1]
		cT_se = 27.8
		kT_se = 0.0047
		LrT_se = 0.964
		NormalizedSeriesElasticElementForce = cT_se * kT_se * np.log(np.exp((LT - LrT_se)/kT_se)+1)
		SEE['Force'].append(NormalizedSeriesElasticElementForce)
	def temp_ib_afferent_activity(SEE):
		SeriesElasticElementForce = SEE['Force'][-1]
		## GTO model
		GTOConstant1 = 60
		GTOConstant2 = 4
		Temp = GTOConstant1*np.log(SeriesElasticElementForce/GTOConstant2+1)
		return(Temp)
	def ib_input(i,SEE,Input):
		# put in num,den here!
		if i%SamplingRatio!=0:
			TempIb1 = Input['TempIb1'][-1]
			Input['TempIb1'].append(TempIb1)
			TempIb2 = Input['TempIb2'][-1]
		elif i==0 or i==SamplingRatio: 
			TempIb1 = temp_ib_afferent_activity(SEE)
			Input['TempIb1'].append(TempIb1)
			TempIb2 = TempIb1
		elif i==2*SamplingRatio: 
			TempIb1 = temp_ib_afferent_activity(SEE)
			Input['TempIb1'].append(TempIb1)
			TempIb2 = (Num[0]*TempIb1/Den[0])	
		elif i==3*SamplingRatio:
			TempIb1 = temp_ib_afferent_activity(SEE)
			Input['TempIb1'].append(TempIb1)
			TempIb2 = (Num[1]*Input['TempIb1'][-2] + Num[0]*Input['TempIb1'][-1] \
					- Den[1]*Input['TempIb2'][-1])/Den[0]
		else:
			TempIb1 = temp_ib_afferent_activity(SEE)
			Input['TempIb1'].append(TempIb1)
			TempIb2 = (Num[2]*Input['TempIb1'][-3] + Num[1]*Input['TempIb1'][-2] + Num[0]*Input['TempIb1'][-1] \
						- Den[2]*Input['TempIb2'][-2] - Den[1]*Input['TempIb2'][-1])/Den[0]
		Input['TempIb2'].append(TempIb2)
		Input['Ib'].append(TempIb2*(TempIb2>0))
	def total_input_with_feedback(i,Input,CE,SEE,delay_parameters,gain_parameters):
		IaSpindleGain = gain_parameters['Ia Gain']
		IISpindleGain = gain_parameters['II Gain']
		IbGTOGain = gain_parameters['Ib Gain']

		IaAfferentDelay = delay_parameters['Ia Delay']
		IIAfferentDelay = delay_parameters['II Delay']
		IbAfferentDelay = delay_parameters['Ib Delay']
		CorticalDelay = delay_parameters['Cortical Delay']
		# Convert delays in ms to samples
		IaAfferentDelayTimeStep = IaAfferentDelay*SamplingRatio   #Ia + II 30
		IIAfferentDelayTimeStep = IIAfferentDelay*SamplingRatio
		IbAfferentDelayTimeStep = IbAfferentDelay*SamplingRatio   #Ib 40
		CorticalDelayTimeStep = CorticalDelay*SamplingRatio   #cortical 50

		IaAfferentDelayTimeStep = int(IaAfferentDelayTimeStep)   #Ia + II 30
		IIAfferentDelayTimeStep = int(IIAfferentDelayTimeStep)
		IbAfferentDelayTimeStep = int(IbAfferentDelayTimeStep)   #Ib 40
		CorticalDelayTimeStep = int(CorticalDelayTimeStep)   #cortical 50

		if i%SamplingRatio==0:
			FeedforwardInput = Input['Feedforward'][i]
			if i >=IaAfferentDelayTimeStep-1: DelayedIaInput = Input['Ia'][i-IaAfferentDelayTimeStep+1]
			if i >=IbAfferentDelayTimeStep+1: DelayedIbInput = Input['Ib'][i-IbAfferentDelayTimeStep-1]
			if i >=IIAfferentDelayTimeStep+1: DelayedIIInput = Input['II'][i-IIAfferentDelayTimeStep-1]
			if i >=CorticalDelayTimeStep+1: DelayedTendonForce = SEE['Force'][i-CorticalDelayTimeStep-1]*CE['MaximumForce']

			if i in range(IaAfferentDelayTimeStep-1,IbAfferentDelayTimeStep+1):
				Input['Total'].append(DelayedIaInput/IaSpindleGain \
										+ FeedforwardInput) #input to the muscle
			elif i in range(IbAfferentDelayTimeStep+1, IIAfferentDelayTimeStep):
				Input['Total'].append(DelayedIaInput/IaSpindleGain \
										-DelayedIbInput/IbGTOGain \
										+FeedforwardInput)  #input to the muscle
			elif i in range(IIAfferentDelayTimeStep+1, CorticalDelayTimeStep+1):
				Input['Total'].append(DelayedIaInput/IaSpindleGain \
										+DelayedIIInput/IISpindleGain \
										-DelayedIbInput/IbGTOGain \
										+FeedforwardInput) #input to the muscle
			elif i >= CorticalDelayTimeStep+1:
				FeedbackInput = Input['Feedback'][-1]
				FeedbackInput = TransCorticalLoopConstant*(Input['TargetTrajectory'][i]-DelayedTendonForce)/CE['MaximumForce'] + FeedbackInput  # feedback input through cortical pathway
				Input['Total'].append(DelayedIaInput/IaSpindleGain \
										+DelayedIIInput/IISpindleGain \
										-DelayedIbInput/IbGTOGain \
										+FeedbackInput) #input to the muscle
			else:	
				Input['Total'].append(FeedforwardInput)
		else:
			Input['Total'].append(Input['Total'][-1])
	def total_input_without_feedback(i,Input,CE,delay_parameters,gain_parameters):
		IaSpindleGain = gain_parameters['Ia Gain']
		IISpindleGain = gain_parameters['II Gain']
		IbGTOGain = gain_parameters['Ib Gain']

		IaAfferentDelayTimeStep = IaAfferentDelay*SamplingRatio   #Ia + II 30
		IIAfferentDelayTimeStep = IIAfferentDelay*SamplingRatio
		IbAfferentDelayTimeStep = IbAfferentDelay*SamplingRatio

		IaAfferentDelayTimeStep = int(IaAfferentDelayTimeStep)   #Ia + II 30
		IIAfferentDelayTimeStep = int(IIAfferentDelayTimeStep)
		IbAfferentDelayTimeStep = int(IbAfferentDelayTimeStep)   #Ib 40

		if i%SamplingRatio==0:
			FeedforwardInput = Input['Feedforward'][i]
			if i >=IaAfferentDelayTimeStep-1: DelayedIaInput = Input['Ia'][i-IaAfferentDelayTimeStep+1]
			if i >=IbAfferentDelayTimeStep+1: DelayedIbInput = Input['Ib'][i-IbAfferentDelayTimeStep-1]
			if i >=IIAfferentDelayTimeStep+1: DelayedIIInput = Input['II'][i-IIAfferentDelayTimeStep-1]

			if i in range(IaAfferentDelayTimeStep-1,IbAfferentDelayTimeStep+1):
				Input['Total'].append(DelayedIaInput/IaSpindleGain \
										+ FeedforwardInput) #input to the muscle
			elif i in range(IbAfferentDelayTimeStep+1, IIAfferentDelayTimeStep):
				Input['Total'].append(DelayedIaInput/IaSpindleGain \
										-DelayedIbInput/IbGTOGain \
										+FeedforwardInput)  #input to the muscle
			elif i >= IIAfferentDelayTimeStep+1:
				Input['Total'].append(DelayedIaInput/IaSpindleGain \
										+DelayedIIInput/IISpindleGain \
										-DelayedIbInput/IbGTOGain \
										+FeedforwardInput) #input to the muscle
			else:
				Input['Total'].append(FeedforwardInput)
		else:
			Input['Total'].append(Input['Total'][-1])
		return(Input)
	def bound(value,min_value,max_value):
		if value > max_value:
			result = 1
		elif value < min_value:
			result = 0
		else:
			result = value
		return(result)
	def delay_total_input(i,Input,delay_parameters):
		EfferentDelayTimeStep = int(delay_parameters['Efferent Delay']*SamplingRatio)
		if i >= EfferentDelayTimeStep:
			LongInput = Input['Total'][i-EfferentDelayTimeStep]
		else:
			LongInput = 0
		return(LongInput)
	def activation_filter(Input):
		TU = 0.12*(Input['LongInput'][-1] < Input['EffectiveMuscleActivation'][-1]) + 0.03 # When true TU = 0.15, else TU = 0.03
		EffectiveMuscleActivationFirstDeriv = (Input['LongInput'][-1] - Input['EffectiveMuscleActivation'][-1])/TU
		Input['EffectiveMuscleActivation'].append(EffectiveMuscleActivationFirstDeriv*(SamplingPeriod) + Input['EffectiveMuscleActivation'][-1]) # effective neural drive
	# Input_1 = initialize_dictionary(['Feedforward','TargetTrajectory','CorticalInput','Feedback','EffectiveMuscleActivation',\
	# 								'Ia','II','Ib','TempIb1','TempIb2','Noise','FilteredNoise','Total','LongInput'],\
	# 								[FeedforwardInput1,TargetForceTrajectory1,CorticalInput1,[FeedbackInput1],[EffectiveMuscleActivation1],\
	# 								[],[],[],[],[],[],[],[],[]])
	# Input_2 = initialize_dictionary(['Feedforward','TargetTrajectory','CorticalInput','Feedback','EffectiveMuscleActivation',\
	# 								'Ia','II','Ib','TempIb1','TempIb2','Noise','FilteredNoise','Total','LongInput'],\
	# 								[FeedforwardInput2,TargetForceTrajectory2,CorticalInput2,[FeedbackInput2],[EffectiveMuscleActivation2],\
	# 								[],[],[],[],[],[],[],[],[]])
	# Muscle_1 = initialize_dictionary(['Length','Velocity','Acceleration','Muscle Force','Tendon Force','Activation Frequency'],\
	# 									[MuscleLength1,MuscleVelocity1,MuscleAcceleration1,[],[],[]])
	# Muscle_2 = initialize_dictionary(['Length','Velocity','Acceleration','Muscle Force','Tendon Force','Activation Frequency'],\
	# 									[MuscleLength2,MuscleVelocity2,MuscleAcceleration2,[],[],[]])
	StartTime = time.time()
	for i in range(len(Time)): #= 1:length(Time)
		if FeedbackOption == 'ff_only': # Feedforward input onl
			append_dictionary(Input_1,['Ia','II','Ib','Total'],[0,0,0,FeedforwardInput1[i]])
			append_dictionary(Input_2,['Ia','II','Ib','Total'],[0,0,0,FeedforwardInput2[i]])
		elif FeedbackOption == 'servo_control': # Servo control (feedforward + spindle and GTO)
			PrimaryOutput1,SecondaryOutput1 = spindle_model(CE_1,Bag1_1,Bag2_1,Chain_1)
			PrimaryOutput2,SecondaryOutput2 = spindle_model(CE_2,Bag1_2,Bag2_2,Chain_2)

			ib_input(i,SEE_1,Input_1)
			ib_input(i,SEE_2,Input_2)

			append_dictionary(Input_1,['Ia','II'],[PrimaryOutput1,SecondaryOutput1])
			append_dictionary(Input_2,['Ia','II'],[PrimaryOutput2,SecondaryOutput2])

			total_input_without_feedback(i,Input_1,CE_1,delay_1_parameters,gain_1_parameters)
			total_input_without_feedback(i,Input_2,CE_2,delay_2_parameters,gain_2_parameters)
		elif FeedbackOption == 'fb_control': # Feedback control (proprioceptive systems + supraspinal loop)
			PrimaryOutput1,SecondaryOutput1 = spindle_model(CE_1,Bag1_1,Bag2_1,Chain_1)
			PrimaryOutput2,SecondaryOutput2 = spindle_model(CE_2,Bag1_2,Bag2_2,Chain_2)
			
			ib_input(i,SEE_1,Input_1)
			ib_input(i,SEE_2,Input_2)
			
			append_dictionary(Input_1,['Ia','II'],[PrimaryOutput1,SecondaryOutput1])
			append_dictionary(Input_2,['Ia','II'],[PrimaryOutput2,SecondaryOutput2])

			total_input_with_feedback(i,Input_1,CE_1,SEE_1,delay_1_parameters,gain_1_parameters)
			total_input_with_feedback(i,Input_2,CE_2,SEE_2,delay_2_parameters,gain_2_parameters)
		# elif FeedbackOption == 'cortical_fb_only': # STILL NEED TO DO THIS ONE!
		# 	if i > CorticalDelayTimeStep-1:
		# 		FeedbackInput1 = TransCorticalLoopConstant*(TargetForceTrajectory1[i]-OutputForceTendon1[i-CorticalDelayTimeStep1-1])/MaximumContractileElementForce1 + FeedbackInput1 # feedback input through cortical pathway
		# 		FeedbackInput2 = TransCorticalLoopConstant*(TargetForceTrajectory2[i]-OutputForceTendon2[i-CorticalDelayTimeStep1-1])/MaximumContractileElementForce2 + FeedbackInput2 # feedback input through cortical pathway
		# 		Input_1.append(FeedbackInput1) #input to the muscle
		# 		Input_2.append(FeedbackInput2) #input to the muscle

		# 	else:
		# 		Input_1.append(FeedforwardInput1[i])
		# 		Input_2.append(FeedforwardInput2[i])

		# 	IaInput1.append(0)
		# 	IIInput1.append(0)
		# 	IbInput1.append(0)
		# 	IaInput2.append(0)
		# 	IIInput2.append(0)
		# 	IbInput2.append(0)

		# Feedforward only
		# Input[i] = input[i] #TargetForceTrajectory[i]/MaximumContractileElementForce
		## Noise + additional input

		Input_1['Total'][-1]=bound(Input_1['Total'][-1],0,1)
		Input_2['Total'][-1]=bound(Input_2['Total'][-1],0,1)

		if i > 4:
			Random1 = 2*(random.random()-0.5)*(np.sqrt(0.01*Input_1['Total'][i])*np.sqrt(3))
			Input_1['Noise'].append(Random1)
			Input_1['FilteredNoise'].append((np.dot(Input_1['Noise'][i-4:i+1],BButtersCoefficients[::-1]) \
													- np.dot(Input_1['FilteredNoise'][i-4:i],AButtersCoefficients[:0:-1]))\
													/AButtersCoefficients[0])
			Random2 = 2*(random.random()-0.5)*(np.sqrt(0.01*Input_2['Total'][i])*np.sqrt(3))
			Input_2['Noise'].append(Random2)
			Input_2['FilteredNoise'].append((np.dot(Input_2['Noise'][i-4:i+1],BButtersCoefficients[::-1]) \
													- np.dot(Input_2['FilteredNoise'][i-4:i],AButtersCoefficients[:0:-1]))\
													/AButtersCoefficients[0])

		else:
			append_dictionary(Input_1,['Noise','FilteredNoise'],[0,0])
			append_dictionary(Input_2,['Noise','FilteredNoise'],[0,0])

		Input_1['Total'][-1] = Input_1['Total'][-1] + Input_1['FilteredNoise'][-1] + Input_1['CorticalInput'][i]
		Input_2['Total'][-1] = Input_2['Total'][-1] + Input_2['FilteredNoise'][-1] + Input_2['CorticalInput'][i]

		Input_1['Total'][-1]=bound(Input_1['Total'][-1],0,1)
		Input_2['Total'][-1]=bound(Input_2['Total'][-1],0,1)

		# add delay along efferent pathway
		Input_1['LongInput'].append(delay_total_input(i,Input_1,delay_1_parameters))
		Input_2['LongInput'].append(delay_total_input(i,Input_2,delay_2_parameters))

		# add activation filter
		activation_filter(Input_1)
		activation_filter(Input_2)

		# SlowTwitch_1 = activation_frequency_slow(CE_1,SlowTwitch_1,EffectiveMuscleActivation1,SamplingFrequency) # not used
		# SlowTwitch_2 = activation_frequency_slow(CE_2,SlowTwitch_2,EffectiveMuscleActivation2,SamplingFrequency) # not used

		# force-velocity relationship
		CE_1['ForceVelocity'].append((CE_1['Velocity'][-1] <= 0)*concentric_force_velocity(CE_1) \
									+ (CE_1['Velocity'][-1] > 0)*eccentric_force_velocity(CE_1)) 
		CE_2['ForceVelocity'].append((CE_2['Velocity'][-1] <= 0)*concentric_force_velocity(CE_2) \
									+ (CE_2['Velocity'][-1] > 0)*eccentric_force_velocity(CE_2)) 
		# force-length relationship
		force_length(CE_1)
		ContractileElementForce1 = CE_1['ForceLength'][-1]*CE_1['ForceVelocity'][-1]
		force_length(CE_2)
		ContractileElementForce2 = CE_2['ForceLength'][-1]*CE_2['ForceVelocity'][-1]
		# viscous property
		ForceViscocity1 = MuscleViscosity * CE_1['Velocity'][-1]
		ForceViscocity2 = MuscleViscosity * CE_2['Velocity'][-1]
		# passive element 1
		parallel_elastic_element_force_1(PEE_1,CE_1)
		parallel_elastic_element_force_1(PEE_2,CE_2)
		# passive element 2
		parallel_elastic_element_force_2(PEE_1,CE_1)
		parallel_elastic_element_force_2(PEE_2,CE_2)


		# total force from contractile element
		ForceTotal1 = (Input_1['EffectiveMuscleActivation'][-1]*(ContractileElementForce1 + PEE_1['ForcePassive2'][-1]) + PEE_1['ForcePassive1'][-1] + ForceViscocity1)*CE_1['MaximumForce']
		ForceTotal1 = ForceTotal1*(ForceTotal1>=0.0)
		Muscle_1['Muscle Force'].append(ForceTotal1)
		ForceTotal2 = (Input_2['EffectiveMuscleActivation'][-1]*(ContractileElementForce2 + PEE_2['ForcePassive2'][-1]) + PEE_2['ForcePassive1'][-1] + ForceViscocity2)*CE_2['MaximumForce']
		ForceTotal2 = ForceTotal2*(ForceTotal2>=0.0)
		Muscle_2['Muscle Force'].append(ForceTotal2)
		#nan_test(ForceTotal,'ForceTotal')

		# force from series elastic element
		normalized_series_elastic_element_force(SEE_1)
		normalized_series_elastic_element_force(SEE_2)
		Muscle_1['Tendon Force'].append(SEE_1['Force'][-1]*CE_1['MaximumForce'])
		Muscle_2['Tendon Force'].append(SEE_2['Force'][-1]*CE_2['MaximumForce'])
		# Muscle_1['Activation Frequency'].append(SlowTwitch_1['ActivationFrequency'])
		# Muscle_2['Activation Frequency'].append(SlowTwitch_2['ActivationFrequency'])
		#nan_test(SeriesElasticElementForce,'SeriesElasticElementForce')

		# calculate muscle excursion acceleration based on the difference
		# between muscle force and tendon force
		if i < len(Time)-1:
			Muscle_1['Acceleration'].append((SEE_1['Force'][-1]*np.cos(PennationAngle1) - Muscle_1['Muscle Force'][-1]*(np.cos(PennationAngle1))**2)/(MuscleMass1) \
				+ (Muscle_1['Velocity'][-1])**2*np.tan(PennationAngle1)**2/(Muscle_1['Length'][-1]))
			Muscle_2['Acceleration'].append((SEE_2['Force'][-1]*np.cos(PennationAngle2) - Muscle_2['Muscle Force'][-1]*(np.cos(PennationAngle2))**2)/(MuscleMass2) \
				+ (Muscle_2['Velocity'][-1])**2*np.tan(PennationAngle2)**2/(Muscle_2['Length'][-1]))
			# integrate acceleration to get velocity
			Muscle_1['Velocity'].append((Muscle_1['Acceleration'][-1]+ \
				Muscle_1['Acceleration'][-2])/2*(SamplingPeriod)+Muscle_1['Velocity'][-1])
			Muscle_2['Velocity'].append((Muscle_2['Acceleration'][-1]+ \
				Muscle_2['Acceleration'][-2])/2*(SamplingPeriod)+Muscle_2['Velocity'][-1])
			# integrate velocity to get length
			Muscle_1['Length'].append((Muscle_1['Velocity'][-1]+ \
				Muscle_1['Velocity'][-2])/2*(SamplingPeriod)+Muscle_1['Length'][-1])
			Muscle_2['Length'].append((Muscle_2['Velocity'][-1]+ \
				Muscle_2['Velocity'][-2])/2*(SamplingPeriod)+Muscle_2['Length'][-1])

			# normalize each variable to optimal muscle length or tendon legnth
			CE_1['Acceleration'].append(Muscle_1['Acceleration'][-1]/(OptimalLength1/100))
			CE_1['Velocity'].append(Muscle_1['Velocity'][-1]/(OptimalLength1/100))
			CE_1['Length'].append(Muscle_1['Length'][-1]/(OptimalLength1/100))
			SEE_1['Length'].append((InitialMusculoTendonLength1 - CE_1['Length'][-1]*OptimalLength1*np.cos(PennationAngle1))/OptimalTendonLength1)
			CE_2['Acceleration'].append(Muscle_2['Acceleration'][-1]/(OptimalLength2/100))
			CE_2['Velocity'].append(Muscle_2['Velocity'][-1]/(OptimalLength2/100))
			CE_2['Length'].append(Muscle_2['Length'][-1]/(OptimalLength2/100))
			SEE_2['Length'].append((InitialMusculoTendonLength2 - CE_2['Length'][-1]*OptimalLength2*np.cos(PennationAngle2))/OptimalTendonLength2)

			# store data
			# OutputSeriesElasticElementLength1.append(SEE_1['Length'])
			# OutputContractileElementLength1.append(CE_1['Length'])
			# OutputContractileElementVelocity1.append(CE_1['Velocity'])
			# OutputContractileElementAcceleration1.append(CE_1['Acceleration'])
			# OutputActivationFrequency1.append(ActivationFrequency1)
			# OutputEffectiveMuscleActivation1.append(Input_1['EffectiveMuscleActivation'][-1])
			# OutputSeriesElasticElementLength2.append(SEE_2['Length'])
			# OutputContractileElementLength2.append(CE_2['Length'])
			# OutputContractileElementVelocity2.append(CE_2['Velocity'])
			# OutputContractileElementAcceleration2.append(CE_2['Acceleration'])
			# OutputActivationFrequency2.append(ActivationFrequency2)
			# OutputEffectiveMuscleActivation2.append(Input_2['EffectiveMuscleActivation'][-1])

		# OutputForceMuscle1.append(ForceTotal1)
		# OutputForceTendon1.append(SEE_1['Force'])
		# OutputForceLength1.append(ForceLength1)
		# OutputForceVelocity1.append(ForceVelocity1)
		# OutputForcePassive1_1.append(ForcePassive1_1)
		# OutputForcePassive2_1.append(ForcePassive2_1)
		# OutputForceMuscle2.append(ForceTotal2)
		# OutputForceTendon2.append(SEE_2['Force'])
		# OutputForceLength2.append(ForceLength2)
		# OutputForceVelocity2.append(ForceVelocity2)
		# OutputForcePassive1_2.append(ForcePassive1_2)
		# OutputForcePassive2_2.append(ForcePassive2_2)
		statusbar(i,len(Time),StartTime=StartTime,Title = 'afferented_muscle_model')

	# plt.figure()
	# plt.plot(Time,SEE_1['Force'])
	# plt.plot(Time,Input_1['TargetTrajectory'],'r')
	# plt.legend(['Output Force','Target Force'])
	# plt.xlabel('Time (sec)')
	# plt.ylabel('Force (N)')

	# plt.figure()
	# plt.plot(Time,SEE_2['Force'])
	# plt.plot(Time,Input_2['TargetTrajectory'],'r')
	# plt.legend(['Output Force','Target Force'])
	# plt.xlabel('Time (sec)')
	# plt.ylabel('Force (N)')

	# save data as output in structure format
	output1 = {	'Input' : Input_1, 'Muscle' : Muscle_1, 'CE' : CE_1, \
				'SEE' : SEE_1, 'PEE' : PEE_1}
	output2 = {	'Input' : Input_2, 'Muscle' : Muscle_2, 'CE' : CE_2, \
				'SEE' : SEE_2, 'PEE' : PEE_2}
	return(output1,output2)

muscle_1_parameters = {	"Pennation Angle":5*np.pi/180, 		"Muscle Mass":0.075,\
						"Optimal Length":10.1, 				"Tendon Length":23.5,\
					 	"Initial Muscle Length":10.1+0.41, 	"Initial Tendon Length":23.5+0.09}
muscle_2_parameters = {	"Pennation Angle":5*np.pi/180, 		"Muscle Mass":0.075,\
						"Optimal Length":10.1, 				"Tendon Length":23.5,\
					 	"Initial Muscle Length":10.1+0.41, 	"Initial Tendon Length":23.5+0.09}

# Define delays based on limb length and conduction velocity of each pathway
DistanceToTheSpinalCord = 0.8 # cm
EfferentConductionVelocity = 48.5 #m/s (ref. Elias et al. 2014) S-type = 44-51, FR = 51-52, FF 52-53
IaConductionVelocity = 64.5 #m/s (ref. Elias et al. 2014)
IIConductionVelocity = 32.5 #m/s (ref. Elias et al. 2014)
IbConductionVelocity = 59 #m/s (ref. Elias et al. 2014)
SynapticDelay = 2 #ms (ref. Kandel)

delay_1_parameters = {"Efferent Delay": round(DistanceToTheSpinalCord/EfferentConductionVelocity*1000), \
						"Ia Delay" : round(DistanceToTheSpinalCord/IaConductionVelocity*1000) + SynapticDelay, \
						"II Delay" : round(DistanceToTheSpinalCord/IIConductionVelocity*1000) + 2*SynapticDelay, \
						"Ib Delay" : round(DistanceToTheSpinalCord/IbConductionVelocity*1000) + 2*SynapticDelay, \
						"Cortical Delay" : 50 }
delay_2_parameters = {"Efferent Delay": round(DistanceToTheSpinalCord/EfferentConductionVelocity*1000), \
						"Ia Delay" : round(DistanceToTheSpinalCord/IaConductionVelocity*1000) + SynapticDelay, \
						"II Delay" : round(DistanceToTheSpinalCord/IIConductionVelocity*1000) + 2*SynapticDelay, \
						"Ib Delay" : round(DistanceToTheSpinalCord/IbConductionVelocity*1000) + 2*SynapticDelay, \
						"Cortical Delay" : 50 }

# Define gain parameters for each neural pathway
gain_1_parameters = {	"Gamma Dynamic Gain" : 70, \
						"Gamma Static Gain" : 70, \
						"Ia Gain" : 800, \
						"II Gain" : 3000, \
						"Ib Gain" : 1000 }
gain_2_parameters = {	"Gamma Dynamic Gain" : 70, \
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
TrajectoryType1 = 'constant' #'triangle','sinwave','trapezoid' (Still need to do Trap)
TrajectoryType2 = 'constant' #'triangle','sinwave','trapezoid' (Still need to do Trap)
TargetTrajectory1 = generate_target_force_trajectory(TrajectoryType1,Time,0.3,0.0,0.0) # in unit of #MVC
TargetTrajectory2 = generate_target_force_trajectory(TrajectoryType2,Time,0.3,0.0,0.0) # in unit of #MVC

# Define additional input to muscle (e.g., 20 Hz oscillatory input)
# CorticalInputAmplitude = (1/10)*TargetTrajectory # in unit of #MVC
# CorticalInputFrequency = 20 # in unit of Hz
# CorticalInput = CorticalInputAmplitude*np.sin(2*np.pi*CorticalInputFrequency*Time)
CorticalInput1 = np.zeros(len(TargetTrajectory1))
CorticalInput2 = np.zeros(len(TargetTrajectory2))

FeedbackType = 'fb_control'
if FeedbackType == 'ff_only':
    FeedbackTypeString = 'FF'
    TargetTrajectory1 = 1.028*TargetTrajectory1
    TargetTrajectory2 = 1.028*TargetTrajectory2
elif FeedbackType == 'servo_control':
    FeedbackTypeString = 'SC'
    TargetTrajectory1 = 0.9285*TargetTrajectory1
    TargetTrajectory2 = 0.9285*TargetTrajectory2
elif FeedbackType == 'fb_control':
    FeedbackTypeString = 'FB'
elif FeedbackType == 'cortical_fb_only':
    FeedbackTypeString = 'CC'

# Output from the muscle model
NumberOfTrials = 1
Output1 = []
Output2 = []
PXX1 = []
PXX2 = []
Range1 = []
Range2 = []
for i in range(NumberOfTrials): #trialN = 1:10
	#rng('shufle')
	output1,output2 = afferented_2_muscles_model(muscle_1_parameters,muscle_2_parameters,\
												delay_1_parameters,delay_2_parameters,\
												gain_1_parameters,gain_2_parameters,
												TargetTrajectory1,TargetTrajectory2,\
												CorticalInput1,CorticalInput2,\
												FeedbackOption = FeedbackType)    
	Output1.append(output1)
	Output2.append(output2)
	f1,pxx1 = welch(output1['SEE']['Force'][-int(10*SamplingFrequency)-1:]-np.average(output1['SEE']['Force'][int(-10*SamplingFrequency-1):]),window = gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5)),\
		noverlap = SamplingFrequency,nperseg = len(gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5))), fs = SamplingFrequency) 
	f2,pxx2 = welch(output2['SEE']['Force'][-int(10*SamplingFrequency)-1:]-np.average(output2['SEE']['Force'][int(-10*SamplingFrequency-1):]),window = gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5)),\
		noverlap = SamplingFrequency,nperseg = len(gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5))), fs = SamplingFrequency) 
	pxx1 = smooth(pxx1,5)
	pxx2 = smooth(pxx2,5)
	PXX1.append(np.array(pxx1,ndmin=2))
	PXX2.append(np.array(pxx2,ndmin=2))
	# Range1.append(max(output1['Ia Input'][-10*SamplingFrequency:-1])-min(output1['Ia Input'][-10*SamplingFrequency-1:]))
	# Range2.append(max(output2['Ia Input'][-10*SamplingFrequency:-1])-min(output2['Ia Input'][-10*SamplingFrequency-1:]))
	print("\n")

PXX1 = np.concatenate(PXX1,axis=0)
plt.figure()
plt.plot(f1[:131],np.average(PXX1[:,:131],axis = 0))

PXX2 = np.concatenate(PXX2,axis=0)
plt.figure()
plt.plot(f2[:131],np.average(PXX2[:,:131],axis = 0))

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