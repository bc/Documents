import numpy as np
from scipy.signal import sawtooth,square,gaussian,welch
import matplotlib.pyplot as plt
import time
import ipdb
import danpy.stuff as dp
from amm_functions import smooth, generate_target_force_trajectory

def afferented_muscle_model(muscle_1_parameters,muscle_2_parameters,\
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
	assert FeedbackOption in ['ff_only','servo_control','fb_control','cortical_fb_only'],\
	"FeedbackOption must be either 'ff_only', 'servo_control', 'fb_control', or 'cortical_fb_only'"
	assert type(muscle_1_parameters)==dict, "muscle_1_parameters must be a dictionary"
	assert len(muscle_1_parameters)==6, "dict muscle_1_parameters can only have 6 entries"
	assert 'Pennation Angle' in muscle_1_parameters, "'Pennation Angle' missing in muscle_1_parameters"
	assert 'Muscle Mass' in muscle_1_parameters, "'Muscle Mass' missing in muscle_1_parameters"
	assert 'Optimal Length' in muscle_1_parameters, "'Optimal Length' missing in muscle_1_parameters"
	assert 'Tendon Length' in muscle_1_parameters, "'Tendon Length' missing in muscle_1_parameters"
	assert 'Initial Muscle Length' in muscle_1_parameters, "'Initial Muscle Length' missing in muscle_1_parameters"
	assert 'Initial Tendon Length' in muscle_1_parameters, "'Initial Tendon Length' missing in muscle_1_parameters"
	assert type(delay_parameters)==dict, "delay_parameters must be a dictionary"
	assert len(delay_parameters)==5, "dict delay_parameters can only have 5 entries"
	assert 'Efferent Delay' in delay_parameters, "'Efferent Delay' missing in delay_parameters"
	assert 'Ia Delay' in delay_parameters, "'Ia Delay' missing in delay_parameters"
	assert 'II Delay' in delay_parameters, "'II Delay' missing in delay_parameters"
	assert 'Ib Delay' in delay_parameters, "'Ib Delay' missing in delay_parameters"
	assert 'Cortical Delay' in delay_parameters, "'Cortical Delay' missing in delay_parameters"
	assert type(gain_parameters)==dict, "gain_parameters must be a dictionary"
	assert len(gain_parameters)==5, "dict gain_parameters can only have 5 entries"
	assert 'Gamma Dynamic Gain' in gain_parameters, "'Gamma Dynamic Gain' missing in gain_parameters"
	assert 'Gamma Static Gain' in gain_parameters, "'Gamma Static Gain' missing in gain_parameters"
	assert 'Ia Gain' in gain_parameters, "'Ia Gain' missing in gain_parameters"
	assert 'II Gain' in gain_parameters, "'II Gain' missing in gain_parameters"
	assert 'Ib Gain' in gain_parameters, "'Ib Gain' missing in gain_parameters"

	assert type(muscle_2_parameters)==dict, "muscle_2_parameters must be a dictionary"
	assert len(muscle_2_parameters)==6, "dict muscle_2_parameters can only have 6 entries"
	assert 'Pennation Angle' in muscle_2_parameters, "'Pennation Angle' missing in muscle_2_parameters"
	assert 'Muscle Mass' in muscle_2_parameters, "'Muscle Mass' missing in muscle_2_parameters"
	assert 'Optimal Length' in muscle_2_parameters, "'Optimal Length' missing in muscle_2_parameters"
	assert 'Tendon Length' in muscle_2_parameters, "'Tendon Length' missing in muscle_2_parameters"
	assert 'Initial Muscle Length' in muscle_2_parameters, "'Initial Muscle Length' missing in muscle_2_parameters"
	assert 'Initial Tendon Length' in muscle_2_parameters, "'Initial Tendon Length' missing in muscle_2_parameters"
	assert type(delay_parameters)==dict, "delay_parameters must be a dictionary"
	assert len(delay_parameters)==5, "dict delay_parameters can only have 5 entries"
	assert 'Efferent Delay' in delay_parameters, "'Efferent Delay' missing in delay_parameters"
	assert 'Ia Delay' in delay_parameters, "'Ia Delay' missing in delay_parameters"
	assert 'II Delay' in delay_parameters, "'II Delay' missing in delay_parameters"
	assert 'Ib Delay' in delay_parameters, "'Ib Delay' missing in delay_parameters"
	assert 'Cortical Delay' in delay_parameters, "'Cortical Delay' missing in delay_parameters"
	assert type(gain_parameters)==dict, "gain_parameters must be a dictionary"
	assert len(gain_parameters)==5, "dict gain_parameters can only have 5 entries"
	assert 'Gamma Dynamic Gain' in gain_parameters, "'Gamma Dynamic Gain' missing in gain_parameters"
	assert 'Gamma Static Gain' in gain_parameters, "'Gamma Static Gain' missing in gain_parameters"
	assert 'Ia Gain' in gain_parameters, "'Ia Gain' missing in gain_parameters"
	assert 'II Gain' in gain_parameters, "'II Gain' missing in gain_parameters"
	assert 'Ib Gain' in gain_parameters, "'Ib Gain' missing in gain_parameters"

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
	OptimalTendonLength1 = TendonLength*1.05
	InitialMuscleLength1 = muscle_1_parameters['Initial Muscle Length']
	InitialTendonLength1 = muscle_1_parameters['Initial Tendon Length']
	InitialMusculoTendonLength1 = InitialMuscleLength1*np.cos(PennationAngle1)+InitialTendonLength1

	PennationAngle2 = muscle_2_parameters['Pennation Angle']
	MuscleMass2 = muscle_2_parameters['Muscle Mass']
	OptimalLength2 = muscle_2_parameters['Optimal Length']
	TendonLength2 = muscle_2_parameters['Tendon Length']
	OptimalTendonLength2 = TendonLength*1.05
	InitialMuscleLength2 = muscle_2_parameters['Initial Muscle Length']
	InitialTendonLength2 = muscle_2_parameters['Initial Tendon Length']
	InitialMusculoTendonLength2 = InitialMuscleLength2*np.cos(PennationAngle2)+InitialTendonLength2

	def return_initial_values(PennationAngle,MuscleMass,OptimalLength,OptimalTendonLength,\
		InitialMuscleLength,InitialTendonLength,InitialMusculoTendonLength):

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

		# dp.compare(PassiveMuscleForce,0.0260)
		# dp.compare(NormalizedSeriesElasticLength,0.9677)

		MaximumMusculoTendonLength = float(OptimalLength * np.cos(PennationAngle) + TendonLength + 0.5)
		NormalizedMaximumFascicleLength = float((MaximumMusculoTendonLength - SeriesElasticLength)/OptimalLength)
		MaximumContractileElementLength = float(NormalizedMaximumFascicleLength/np.cos(PennationAngle)) # Contractile Element = Muscle

		# dp.compare(MaximumMusculoTendonLength,34.0616)
		# dp.compare(NormalizedSeriesElasticLength,0.9677)

		# dp.nan_test(PassiveMuscleForce,'PassiveMuscleForce')
		# dp.nan_test(NormalizedSeriesElasticLength,'NormalizedSeriesElasticLength')
		# dp.nan_test(MaximumMusculoTendonLength,'MaximumMusculoTendonLength')
		# dp.nan_test(SeriesElasticLength,'SeriesElasticLength')
		# dp.nan_test(NormalizedMaximumFascicleLength,'NormalizedMaximumFascicleLength')
		# dp.nan_test(MaximumContractileElementLength,'MaximumContractileElementLength')

		InitialLength = float((InitialMusculoTendonLength+OptimalTendonLength\
								*((TendonRateConstant/MuscleRateConstant)*RestingMuscleLength-NormalizedRestingTendonLength\
								- TendonRateConstant*np.log((MuscleConstant/TendonConstant)*(MuscleRateConstant/TendonRateConstant))))\
								/(100*(1+TendonRateConstant/MuscleRateConstant*OptimalTendonLength/MaximumContractileElementLength*(1/OptimalLength))*np.cos(PennationAngle)))
		ContractileElementLength = float(InitialLength/(OptimalLength/100)) # initializing CE Length
		SeriesElasticElementLength = float((InitialMusculoTendonLength - InitialLength*100)/OptimalTendonLength) # MTU - initial muscle = initial SEE length

		# dp.compare(InitialLength,0.1066)
		# dp.nan_test(InitialLength,'InitialLength')
		# dp.nan_test(ContractileElementLength,'ContractileElementLength')
		# dp.nan_test(SeriesElasticElementLength,'SeriesElasticElementLength')

		# calculate maximal force output of the muscle based on PCSA
		MuscleDensity = 1.06 #g/cm**3
		PCSA = (MuscleMass*1000)/MuscleDensity/OptimalLength #physionp.logical cross-sectional area
		SpecificTension = 31.4 # WHERE DOES THIS COME FROM AND WHAT DOES IT MEAN? N/
		MaximumContractileElementForce = PCSA * SpecificTension

		# parameter initialization
		ContractileElementVelocity = 0
		ContractileElementAcceleration = 0
		MuscleAcceleration = [0]
		MuscleVelocity = [0]
		MuscleLength = [ContractileElementLength*OptimalLength/100]
		SeriesElasticElementForce = 0.105

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
		Bag1Tension = 0
		Bag2TensionSecondDeriv = 0
		Bag2TensionFirstDeriv = 0
		Bag2Tension = 0
		ChainTensionSecondDeriv = 0
		ChainTensionFirstDeriv = 0
		ChainTension = 0

		Input, IaInput, IIInput = [],[],[]
		IbInput, x, TemporaryIbInput = [],[],[]
		Noise, FilteredNoise, LongInput = [],[],[]
		OutputForceMuscle,OutputForceTendon,OutputForceLength = [], [], []
		OutputForceVelocity,OutputForcePassive1,OutputForcePassive2 = [], [], []
		OutputSeriesElasticElementLength,OutputContractileElementVelocity = [(InitialMusculoTendonLength - InitialLength*100)/OptimalTendonLength], []
		OutputContractileElementLength = [InitialLength/(OptimalLength/100)]
		OutputContractileElementAcceleration,OutputActivationFrequency,OutputEffectiveMuscleActivation = [], [], []

		return(PassiveMuscleForce,NormalizedSeriesElasticLength,SeriesElasticLength,\
				MaximumMusculoTendonLength,NormalizedMaximumFascicleLength,MaximumContractileElementLength, \
				InitialLength,ContractileElementLength,SeriesElasticElementLength, \
				MaximumContractileElementForce,ContractileElementVelocity,ContractileElementAcceleration, \
				MuscleAcceleration,MuscleVelocity,MuscleLength, \
				SeriesElasticElementForce,Y_dot,Y, \
				Saf_dot,Saf,fint_dot, \
				fint,feff_dot,feff, \
				ActivationFrequency,EffectiveMuscleActivation,DynamicSpindleFrequency, \
				StaticSpindleFrequency,Bag1TensionSecondDeriv,Bag1TensionFirstDeriv, \
				Bag1Tension,Bag2TensionSecondDeriv,Bag2TensionFirstDeriv, \
				Bag2Tension,ChainTensionSecondDeriv,ChainTensionFirstDeriv, \
				ChainTension,Input, IaInput, \
				IIInput, IbInput, x, \
				TemporaryIbInput, Noise, FilteredNoise, \
				LongInput,OutputForceMuscle,OutputForceTendon,\
				OutputForceLength,OutputForceVelocity,OutputForcePassive1,\
				OutputForcePassive2,OutputSeriesElasticElementLength,OutputContractileElementVelocity,\
				OutputContractileElementLength,OutputContractileElementAcceleration,OutputActivationFrequency,\
				OutputEffectiveMuscleActivation)

	PassiveMuscleForce1,NormalizedSeriesElasticLength1,SeriesElasticLength1,\
			MaximumMusculoTendonLength1,NormalizedMaximumFascicleLength1,MaximumContractileElementLength1, \
			InitialLength1,ContractileElementLength1,SeriesElasticElementLength1, \
			MaximumContractileElementForce1,ContractileElementVelocity1,ContractileElementAcceleration1, \
			MuscleAcceleration1,MuscleVelocity1,MuscleLength1, \
			SeriesElasticElementForce1,Y_dot1,Y1, \
			Saf_dot1,Saf1,fint_dot1, \
			fint1,feff_dot1,feff1, \
			ActivationFrequency1,EffectiveMuscleActivation1,DynamicSpindleFrequency1, \
			StaticSpindleFrequency1,Bag1TensionSecondDeriv1,Bag1TensionFirstDeriv1, \
			Bag1Tension1,Bag2TensionSecondDeriv1,Bag2TensionFirstDeriv1, \
			Bag2Tension1,ChainTensionSecondDeriv1,ChainTensionFirstDeriv1, \
			ChainTension1, Input1, IaInput1, \
			IIInput1, IbInput1, x1, \
			TemporaryIbInput1, Noise1, FilteredNoise1, \
			LongInput1,OutputForceMuscle1,OutputForceTendon1,\
			OutputForceLength1,OutputForceVelocity1,OutputForcePassive1_1,\
			OutputForcePassive2_1,OutputSeriesElasticElementLength1,OutputContractileElementVelocity1,\
			OutputContractileElementLength1,OutputContractileElementAcceleration1,OutputActivationFrequency1,\
			OutputEffectiveMuscleActivation1 \
				= return_initial_values(PennationAngle1,MuscleMass1,OptimalLength1,OptimalTendonLength1,\
										InitialMuscleLength1,InitialTendonLength1,InitialMusculoTendonLength1)
	PassiveMuscleForce2,NormalizedSeriesElasticLength2,SeriesElasticLength2,\
			MaximumMusculoTendonLength2,NormalizedMaximumFascicleLength2,MaximumContractileElementLength2, \
			InitialLength2,ContractileElementLength2,SeriesElasticElementLength2, \
			MaximumContractileElementForce2,ContractileElementVelocity2,ContractileElementAcceleration2, \
			MuscleAcceleration2,MuscleVelocity2,MuscleLength2, \
			SeriesElasticElementForce2,Y_dot2,Y2, \
			Saf_dot2,Saf2,fint_dot2, \
			fint2,feff_dot2,feff2, \
			ActivationFrequency2,EffectiveMuscleActivation2,DynamicSpindleFrequency2, \
			StaticSpindleFrequency2,Bag2TensionSecondDeriv2,Bag2TensionFirstDeriv2, \
			Bag2Tension2,Bag2TensionSecondDeriv2,Bag2TensionFirstDeriv2, \
			Bag2Tension2,ChainTensionSecondDeriv2,ChainTensionFirstDeriv2, \
			ChainTension2, Input2, IaInput2, \
			IIInput2, IbInput2, x2, \
			TemporaryIbInput2, Noise2, FilteredNoise2, \
			LongInput2,OutputForceMuscle2,OutputForceTendon2,\
			OutputForceLength2,OutputForceVelocity2,OutputForcePassive1_2,\
			OutputForcePassive2_2,OutputSeriesElasticElementLength2,OutputContractileElementVelocity2,\
			OutputContractileElementLength2,OutputContractileElementAcceleration2,OutputActivationFrequency2,\
			OutputEffectiveMuscleActivation2 \
				= return_initial_values(PennationAngle2,MuscleMass2,OptimalLength2,OptimalTendonLength2,\
										InitialMuscleLength2,InitialTendonLength2,InitialMusculoTendonLength2)
	#dp.nan_test(PCSA,'PCSA')

	# visnp.cosity of muscle (ref: Elias et al. 2014)
	MuscleViscosity = 0.005 #0.001

	# assign delays
	EfferentDelay1 = delay_1_parameters['Efferent Delay']
	IaAfferentDelay1 = delay_1_parameters['Ia Delay']
	IIAfferentDelay1 = delay_1_parameters['II Delay']
	IbAfferentDelay1 = delay_1_parameters['Ib Delay']
	CorticalDelay1 = delay_1_parameters['Cortical Delay']

	EfferentDelay2 = delay_2_parameters['Efferent Delay']
	IaAfferentDelay2 = delay_2_parameters['Ia Delay']
	IIAfferentDelay2 = delay_2_parameters['II Delay']
	IbAfferentDelay2 = delay_2_parameters['Ib Delay']
	CorticalDelay2 = delay_2_parameters['Cortical Delay']

	SamplingFrequency = 10000
	FeedbackSamplingFrequency = 1000
	SamplingRatio = SamplingFrequency/FeedbackSamplingFrequency
	SamplingPeriod = 1/SamplingFrequency
	Time = np.arange(0,len(TargetTrajectory)*SamplingPeriod,SamplingPeriod) # 0:SamplingPeriod:(length(TargetTrajectory)-1)/SamplingFrequency

	# filter parameters for Noise
	BButtersCoefficients,AButtersCoefficients = signal.butter(4,100/(SamplingFrequency/2),'low')

	# dp.compare(BButtersCoefficients,np.array([0.0898e-05, 0.3594e-05, 0.5391e-05, 0.3594e-05, 0.0898e-05]))
	# dp.compare(AButtersCoefficients,np.array([1.0000, -3.8358,  5.5208, -3.5335,  0.8486]))

	## GTO model
	GTOConstant1 = 60
	GTOConstant2 = 4

	Num,Den = [1.7,2.58,0.4],[1,2.2,0.4]
	ContinuousTransferFunction = control.tf(Num,Den)
	DiscreteTransferFunction = control.matlab.c2d(ContinuousTransferFunction,1/FeedbackSamplingFrequency)
	Num,Den = control.matlab.tfdata(DiscreteTransferFunction)
	Num,Den = Num[0][0],Den[0][0]

	# dp.compare(Num,np.array([1.7000,-3.3974, 1.6974]))
	# dp.compare(Den,np.array([1.0000,-1.9978,0.9978]))

	## Transcortical loop
	TransCorticalLoopConstant = 0.01

	## Gains
	GammaDynamicGain1 = gain_1_parameters['Gamma Dynamic Gain']
	GammaStaticGain1 = gain_1_parameters['Gamma Static Gain']
	IaSpindleGain1 = gain_1_parameters['1Ia Gain']
	IISpindleGain1 = gain_1_parameters['II Gain']
1	IbGTOGain1 =1 gain_1_parameters['Ib Gain']1

1	GammaDynamicGain2 = gain_2_parameters1['Gamma Dynamic Gain1']1
	GammaStaticGain2 =1 gain_2_parameters['Gamma Static1 Gain']
	IaSpindleGain2 = gain_2_parameters['Ia Gain']
	IISpindleGain21 = gain_2_parameters1['II Gain']
	IbGTOGain21 = gain_2_parameters1['Ib Gain'1]

1	# Convert1 delays in ms1 to samples
1	IaAfferentDelayTimeStep1 =1 int(IaAfferentDelay1*SamplingRatio1)   #Ia +1 II 30
	IIAfferentDelayTimeStep1 = int(IIAfferentDelay1*SamplingRatio)
	IbAfferentDelayTimeStep1 = int(IbAfferentDelay1*SamplingRatio)   #Ib 40
	CorticalDelayTimeStep1 = int(CorticalDelay1*SamplingRatio)   #cortical 50
	EfferentDelayTimeStep1 = int(EfferentDelay1*SamplingRatio)

	IaAfferentDelayTimeStep2 = int(IaAfferentDelay2*SamplingRatio)   #Ia + II 30
	IIAfferentDelayTimeStep2 = int(IIAfferentDelay2*SamplingRatio)
	IbAfferentDelayTimeStep2 = int(IbAfferentDelay2*SamplingRatio)   #Ib 40
	CorticalDelayTimeStep2 = int(CorticalDelay2*SamplingRatio)   #cortical 50
	EfferentDelayTimeStep2 = int(EfferentDelay2*SamplingRatio)

	# assert IaAfferentDelayTimeStep1 == IaAfferentDelay1*SamplingRatio, "Delay Time Step Mismatch!"
	# assert IIAfferentDelayTimeStep1 == IIAfferentDelay1*SamplingRatio, "Delay Time Step Mismatch!"
	# assert IbAfferentDelayTimeStep1 == IbAfferentDelay1*SamplingRatio, "Delay Time Step Mismatch!"
	# assert CorticalDelayTimeStep1 == CorticalDelay1*SamplingRatio, "Delay Time Step Mismatch!"
	# assert EfferentDelayTimeStep1 == EfferentDelay1*SamplingRatio, "Delay Time Step Mismatch!"

	# assert IaAfferentDelayTimeStep2 == IaAfferentDelay2*SamplingRatio, "Delay Time Step Mismatch!"
	# assert IIAfferentDelayTimeStep2 == IIAfferentDelay2*SamplingRatio, "Delay Time Step Mismatch!"
	# assert IbAfferentDelayTimeStep2 == IbAfferentDelay2*SamplingRatio, "Delay Time Step Mismatch!"
	# assert CorticalDelayTimeStep2 == CorticalDelay2*SamplingRatio, "Delay Time Step Mismatch!"
	# assert EfferentDelayTimeStep2 == EfferentDelay2*SamplingRatio, "Delay Time Step Mismatch!"

	FeedbackInput1, FeedbackInput2 = 0,0
	Count = 0

	# Convert force trajectory to unit of newton
	TargetForceTrajectory1 = TargetTrajectory1*MaximumContractileElementForce1
	TargetForceTrajectory2 = TargetTrajectory2*MaximumContractileElementForce2
	FeedforwardInput1 = TargetForceTrajectory1/MaximumContractileElementForce1
	FeedforwardInput2 = TargetForceTrajectory2/MaximumContractileElementForce2

	def bag1_model(Length,LengthFirstDeriv,LengthSecondDeriv,DynamicSpindleFrequency,Bag1TensionFirstDeriv,Bag1Tension):
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

		#dp.nan_test(DynamicSpindleFrequencyFirstDeriv,'DynamicSpindleFrequencyFirstDeriv')
		#dp.nan_test(DynamicSpindleFrequency,'DynamicSpindleFrequency')

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

		#dp.nan_test(Bag1TensionSecondDeriv,'Bag1TensionSecondDeriv')
		#dp.nan_test(Bag1TensionFirstDeriv,'Bag1TensionFirstDeriv')
		#dp.nan_test(Bag1Tension,'Bag1Tension')

		Bag1AfferentPotential = G*(Bag1Tension/K_SR-(LN_SR-L0_SR))

		#dp.nan_test(Bag1AfferentPotential,'Bag1AfferentPotential')

		return(Bag1AfferentPotential,DynamicSpindleFrequency,Bag1TensionFirstDeriv,Bag1Tension)
	def bag2_model(Length,LengthFirstDeriv,LengthSecondDeriv,StaticSpindleFrequency,Bag2TensionFirstDeriv,Bag2Tension):
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

		#dp.nan_test(StaticSpindleFrequencyFirstDeriv,'StaticSpindleFrequencyFirstDeriv')
		#dp.nan_test(StaticSpindleFrequency,'StaticSpindleFrequency')

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

		#dp.nan_test(Bag2TensionSecondDeriv,'Bag2TensionSecondDeriv')
		#dp.nan_test(Bag2TensionFirstDeriv,'Bag2TensionFirstDeriv')
		#dp.nan_test(Bag2Tension,'Bag2Tension')

		Bag2PrimaryAfferentPotential = G*(Bag2Tension/K_SR-(LN_SR-L0_SR))
		Bag2SecondaryAfferentPotential = G*(	X*(L_secondary/L0_SR)*(Bag2Tension/K_SR-(LN_SR-L0_SR)) \
							+(1-X)*(L_secondary/L0_PR)*(Length-Bag2Tension/K_SR-(L0_SR+LN_PR))	)

		#dp.nan_test(Bag2PrimaryAfferentPotential,'Bag2PrimaryAfferentPotential')
		#dp.nan_test(Bag2SecondaryAfferentPotential,'Bag2SecondaryAfferentPotential')

		return(Bag2PrimaryAfferentPotential,Bag2SecondaryAfferentPotential,StaticSpindleFrequency,Bag2TensionFirstDeriv,Bag2Tension)
	def chain_model(Length,LengthFirstDeriv,LengthSecondDeriv,ChainTensionFirstDeriv,ChainTension):
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

		#dp.nan_test(ChainStaticSpindleFrequency,'ChainStaticSpindleFrequency')

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

		#dp.nan_test(ChainTensionSecondDeriv,'ChainTensionSecondDeriv')
		#dp.nan_test(ChainTensionFirstDeriv,'ChainTensionFirstDeriv')
		#dp.nan_test(ChainTension,'ChainTension')

		ChainPrimaryAfferentPotential = G*(ChainTension/K_SR-(LN_SR-L0_SR))
		ChainSecondaryAfferentPotential = G*(	X*(L_secondary/L0_SR)*(ChainTension/K_SR-(LN_SR-L0_SR)) \
							+ (1-X)*(L_secondary/L0_PR)*(Length-ChainTension/K_SR-(L0_SR+LN_PR))	)

		#dp.nan_test(ChainPrimaryAfferentPotential,'ChainPrimaryAfferentPotential')
		#dp.nan_test(ChainSecondaryAfferentPotential,'ChainSecondaryAfferentPotential')

		return(ChainPrimaryAfferentPotential,ChainSecondaryAfferentPotential,ChainTensionFirstDeriv,ChainTension)
	def spindle_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,\
		DynamicSpindleFrequency,Bag1TensionFirstDeriv,Bag1Tension,\
		StaticSpindleFrequency,Bag2TensionFirstDeriv,Bag2Tension,\
		ChainTensionFirstDeriv,ChainTension):

		S = 0.156

		Bag1AfferentPotential,DynamicSpindleFrequency,Bag1TensionFirstDeriv,Bag1Tension = bag1_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,DynamicSpindleFrequency,Bag1TensionFirstDeriv,Bag1Tension)
		Bag2PrimaryAfferentPotential,Bag2SecondaryAfferentPotential,StaticSpindleFrequency,Bag2TensionFirstDeriv,Bag2Tension = bag2_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,StaticSpindleFrequency,Bag2TensionFirstDeriv,Bag2Tension)
		ChainPrimaryAfferentPotential,ChainSecondaryAfferentPotential,ChainTensionFirstDeriv,ChainTension = chain_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,ChainTensionFirstDeriv,ChainTension)

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

		#dp.nan_test(PrimaryOutput,'PrimaryOutput')
		#dp.nan_test(SecondaryOutput,'SecondaryOutput')

		return(PrimaryOutput,SecondaryOutput)
	def activation_frequency_slow(Activation,Length,LengthFirstDeriv,Y,fint,feff_dot,feff,SamplingFrequency,ActivationFrequencySlow):
	
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

		#dp.nan_test(Y_dot,'Y_dot')
		#dp.nan_test(Y,'Y')

		fenv = (fmax-fmin)/(1-Uth)*Activation+fmin-(fmax-fmin)*Uth
		fenv = fenv/f_half

		#dp.nan_test(fenv,'fenv')

		if feff_dot >= 0:
		    Tf = Tf1 * Length**2 + Tf2 * fenv
		else:
		    Tf = (Tf3 + Tf4*ActivationFrequencySlow)/Length

		#dp.nan_test(Tf,'Tf')

		fint_dot = (fenv - fint)/Tf
		fint = fint_dot*SamplingPeriod + fint
		feff_dot = (fint - feff)/Tf
		feff = feff_dot*SamplingPeriod + feff

		#dp.nan_test(fint_dot,'fint_dot')
		#dp.nan_test(fint,'fint')
		#dp.nan_test(feff_dot,'feff_dot')
		#dp.nan_test(feff,'feff')

		nf = nf0 + nf1*(1/Length-1)
		ActivationFrequencySlow = 1 - np.exp(-(Y*feff/(af*nf))**nf)

		#dp.nan_test(nf,'nf')
		#dp.nan_test(ActivationFrequencySlow,'ActivationFrequencySlow')

		return(ActivationFrequencySlow,Y,fint,feff_dot,feff)
	def activation_frequency_fast(Activation,Length,Saf,fint,feff_dot,feff,SamplingFrequency,ActivationFrequencyFast):
		
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

		#dp.nan_test(fenv,'fenv')

		if feff_dot >= 0:
			Tf = Tf1 * Length**2 + Tf2 * fenv
		elif feff_dot < 0:
			Tf = (Tf3 + Tf4*ActivationFrequencyFast)/Length

		#dp.nan_test(Tf,'Tf')

		if feff < 0.1:
			AS = as1
		elif feff >= 0.1:
			AS = as2

		#dp.nan_test(AS,'AS')

		Saf_dot = (AS - Saf)/Ts
		Saf = Saf_dot*SamplingPeriod + Saf
		fint_dot = (fenv - fint)/Tf
		fint = fint_dot*SamplingPeriod + fint
		feff_dot = (fint - feff)/Tf
		feff = feff_dot*SamplingPeriod + feff
		nf = nf0 + nf1*(1/Length-1)
		ActivationFrequencyFast = 1 - np.exp(-(S_af*feff/(af*nf))**nf)

		#dp.nan_test(Saf_dot,'Saf_dot')
		#dp.nan_test(Saf,'Saf')
		#dp.nan_test(fint_dot,'fint_dot')
		#dp.nan_test(fint,'fint')
		#dp.nan_test(feff_dot,'feff_dot')
		#dp.nan_test(feff,'feff')
		#dp.nan_test(nf,'nf')
		#dp.nan_test(ActivationFrequencyFast,'ActivationFrequencyFast')

		return(ActivationFrequencyFast,fint,feff_dot,feff,Saf)
	def force_length(Length):
		beta = 2.3
		omega = 1.12
		rho = 1.62
		ForceLength = np.exp(-abs((Length**beta - 1)/omega)**rho)
		#dp.nan_test(ForceLength,'ForceLength')
		return(ForceLength)
	def concentric_force_velocity(Length,V):
		Vmax = -7.88
		cv0 = 5.88
		cv1 = 0
		ConcentricForceVelocity = (Vmax - V)/(Vmax + (cv0 + cv1*Length)*V)
		#dp.nan_test(ConcentricForceVelocity,'ConcentricForceVelocity')
		return(ConcentricForceVelocity)
	def eccentric_force_velocity(Length,V):
		av0 = -4.7
		av1 = 8.41
		av2 = -5.34
		bv = 0.35
		EccentricForceVelocity = (bv - (av0 + av1*Length + av2*Length**2)*V)/(bv+V)
		#dp.nan_test(EccentricForceVelocity,'EccentricForceVelocity')
		return(EccentricForceVelocity)
	def parallel_elastic_element_force_1(Length,MaximumContractileElementLength):
		c1_pe1 = 23.0  #355, 67.1
		k1_pe1 = 0.046 #0.04, 0.056
		Lr1_pe1 = 1.17  #1.35, 1.41
		ParallelElasticElementForce1 = c1_pe1 * k1_pe1 * np.log(np.exp((Length/MaximumContractileElementLength - Lr1_pe1)/k1_pe1)+1)
		#dp.nan_test(ParallelElasticElementForce1,'ParallelElasticElementForce1')
		return(ParallelElasticElementForce1)
	def parallel_elastic_element_force_2(Length):
		c2_pe2 = -0.02 #0.01  -0.1
		k2_pe2 = -21
		Lr2_pe2 = 0.70 #0.79 0.59
		ParallelElasticElementForce2 = c2_pe2*np.exp((k2_pe2*(Length-Lr2_pe2))-1)
		#dp.nan_test(ParallelElasticElementForce2,'ParallelElasticElementForce2')
		return(ParallelElasticElementForce2)
	def normalized_series_elastic_element_force(LT):
		cT_se = 27.8
		kT_se = 0.0047
		LrT_se = 0.964
		NormalizedSeriesElasticElementForce = cT_se * kT_se * np.log(np.exp((LT - LrT_se)/kT_se)+1)
		#dp.nan_test(NormalizedSeriesElasticElementForce,'NormalizedSeriesElasticElementForce')
		return(NormalizedSeriesElasticElementForce)
	
	StartTime = time.time()
	for i in range(len(Time)): #= 1:length(Time)
		if FeedbackOption == 'ff_only': # Feedforward input onl
			Input1.append(FeedforwardInput1[i])
			IaInput1.append(0)
			IIInput1.append(0)
			IbInput1.append(0)

			Input2.append(FeedforwardInput2[i])
			IaInput2.append(0)
			IIInput2.append(0)
			IbInput2.append(0)
		elif FeedbackOption == 'servo_control': # Servo control (feedforward + spindle and GTO)
			PrimaryOutput1,SecondaryOutput1 = spindle_model(ContractileElementLength1,ContractileElementVelocity1,ContractileElementAcceleration1,\
				DynamicSpindleFrequency1,Bag1TensionFirstDeriv1,Bag1Tension1,\
				StaticSpindleFrequency1,Bag2TensionFirstDeriv1,Bag2Tension1,\
				ChainTensionFirstDeriv1,ChainTension1)
			PrimaryOutput2,SecondaryOutput2 = spindle_model(ContractileElementLength2,ContractileElementVelocity2,ContractileElementAcceleration2,\
				DynamicSpindleFrequency2,Bag1TensionFirstDeriv2,Bag1Tension2,\
				StaticSpindleFrequency2,Bag2TensionFirstDeriv2,Bag2Tension2,\
				ChainTensionFirstDeriv2,ChainTension2)
			IaInput1.append(PrimaryOutput1)
			IIInput1.append(SecondaryOutput1)
			IaInput2.append(PrimaryOutput2)
			IIInput2.append(SecondaryOutput2)
			if i == Count:
				# Get spindle primary and secondary afferent activity

				# Get Ib activity
				x1.append(GTOConstant1*np.log(SeriesElasticElementForce1/GTOConstant2+1))
				x2.append(GTOConstant1*np.log(SeriesElasticElementForce2/GTOConstant2+1))
				if i == 0:
					TemporaryIbInput1.append(x1[-1])
					TemporaryIbInput2.append(x2[-1])
				elif i == 1*SamplingRatio:
					TemporaryIbInput1.append(x1[-1])
					TemporaryIbInput2.append(x2[-1])
				elif i == 2*SamplingRatio:
					TemporaryIbInput1.append((Num[0]*x1[-1])/Den[0])
					TemporaryIbInput2.append((Num[0]*x2[-1])/Den[0])
				elif i == 3*SamplingRatio:
					TemporaryIbInput1.append((Num[1]*x1[-2] + Num[0]*x1[-1] \
						- Den[1]*TemporaryIbInput1[-2])/Den[0])
					TemporaryIbInput2.append((Num[1]*x2[-2] + Num[0]*x2[-1] \
						- Den[1]*TemporaryIbInput2[-2])/Den[0])
				elif i >= 3*SamplingRatio:
					TemporaryIbInput1.append((Num[2]*x1[-3] + Num[1]*x1[-2] + Num[0]*x1[-1] \
						- Den[2]*TemporaryIbInput1[-3] - Den[1]*TemporaryIbInput1[-2])/Den[0])
					TemporaryIbInput2.append((Num[2]*x2[-3] + Num[1]*x2[-2] + Num[0]*x2[-1] \
						- Den[2]*TemporaryIbInput2[-3] - Den[1]*TemporaryIbInput2[-2])/Den[0])

				IbInput1.append(TemporaryIbInput1[-1]*(TemporaryIbInput1[-1]>0))
				IbInput2.append(TemporaryIbInput2[-1]*(TemporaryIbInput2[-1]>0))
			
				if i in range(IaAfferentDelayTimeStep1-1,IbAfferentDelayTimeStep1):
					Input1.append(IaInput1[i-IaAfferentDelayTimeStep1-1]/IaSpindleGain1 \
								+FeedforwardInput1[i]) #input to the muscle
				elif i in range(IbAfferentDelayTimeStep1, IIAfferentDelayTimeStep1):
					Input1.append(IaInput1[i-IaAfferentDelayTimeStep1-1]/IaSpindleGain1 \
								-IbInput1[i-IbAfferentDelayTimeStep1-1]/IbGTOGain1 \
								+FeedforwardInput1[i]) #input to the muscle
				elif i in range(IIAfferentDelayTimeStep1, CorticalDelayTimeStep1):
					Input1.append(IaInput1[i-IaAfferentDelayTimeStep1-1]/IaSpindleGain1 \
								+IIInput1[i-IIAfferentDelayTimeStep1-1]/IISpindleGain1 \
								-IbInput1[i-IbAfferentDelayTimeStep1-1]/IbGTOGain1 \
								+FeedforwardInput1[i]) #input to the muscle
				else:	
					Input1.append(FeedforwardInput1[i])

				if i in range(IaAfferentDelayTimeStep2-1,IbAfferentDelayTimeStep2):
					Input2.append(IaInput2[i-IaAfferentDelayTimeStep2-1]/IaSpindleGain2 \
								+FeedforwardInput2[i]) #input to the muscle
				elif i in range(IbAfferentDelayTimeStep2, IIAfferentDelayTimeStep2):
					Input2.append(IaInput2[i-IaAfferentDelayTimeStep2-1]/IaSpindleGain2 \
								-IbInput2[i-IbAfferentDelayTimeStep2-1]/IbGTOGain2 \
								+FeedforwardInput2[i]) #input to the muscle
				elif i in range(IIAfferentDelayTimeStep2, CorticalDelayTimeStep2):
					Input2.append(IaInput2[i-IaAfferentDelayTimeStep2-1]/IaSpindleGain2 \
								+IIInput2[i-IIAfferentDelayTimeStep2-1]/IISpindleGain2 \
								-IbInput2[i-IbAfferentDelayTimeStep2-1]/IbGTOGain2 \
								+FeedforwardInput2[i]) #input to the muscle
				else:	
					Input2.append(FeedforwardInput2[i])
				
				Count = SamplingRatio+Count
			else:
				IbInput1.append(IbInput1[-1])
				Input1.append(Input1[-1])
				IbInput2.append(IbInput2[-1])
				Input2.append(Input2[-1])
		
		elif FeedbackOption == 'fb_control': # Feedback control (proprioceptive systems + supraspinal loop)
			PrimaryOutput1,SecondaryOutput1 = spindle_model(ContractileElementLength1,ContractileElementVelocity1,ContractileElementAcceleration1,\
				DynamicSpindleFrequency1,Bag1TensionFirstDeriv1,Bag1Tension1,\
				StaticSpindleFrequency1,Bag2TensionFirstDeriv1,Bag2Tension1,\
				ChainTensionFirstDeriv1,ChainTension1)
			PrimaryOutput2,SecondaryOutput2 = spindle_model(ContractileElementLength2,ContractileElementVelocity2,ContractileElementAcceleration2,\
				DynamicSpindleFrequency2,Bag1TensionFirstDeriv2,Bag1Tension2,\
				StaticSpindleFrequency2,Bag2TensionFirstDeriv2,Bag2Tension2,\
				ChainTensionFirstDeriv2,ChainTension2)
			IaInput1.append(PrimaryOutput1)
			IIInput1.append(SecondaryOutput1)
			IaInput2.append(PrimaryOutput2)
			IIInput2.append(SecondaryOutput2)
			# if i == 0: dp.compare([PrimaryOutput,SecondaryOutput],[0,7.5646])

			if i == Count:
				# Get spindle primary and secondary afferent activity
				
				# Get Ib activity
				x1.append(GTOConstant1*np.log((SeriesElasticElementForce1/GTOConstant2+1)))
				x2.append(GTOConstant1*np.log((SeriesElasticElementForce2/GTOConstant2+1)))
				if i == 0:
					TemporaryIbInput1.append(x1[-1])
					TemporaryIbInput2.append(x2[-1])
				elif i == 1*SamplingRatio:
					TemporaryIbInput1.append(x1[-1])
					TemporaryIbInput2.append(x2[-1])
				elif i == 2*SamplingRatio:
					TemporaryIbInput1.append((Num[0]*x1[-1])/Den[0])
					TemporaryIbInput2.append((Num[0]*x2[-1])/Den[0])
				elif i == 3*SamplingRatio:
					TemporaryIbInput1.append((Num[1]*x1[-2] + Num[0]*x1[-1] \
						- Den[1]*TemporaryIbInput1[-2])/Den[0])
					TemporaryIbInput2.append((Num[1]*x2[-2] + Num[0]*x2[-1] \
						- Den[1]*TemporaryIbInput2[-2])/Den[0])
				elif i >= 3*SamplingRatio:
					TemporaryIbInput1.append((Num[2]*x1[-3] + Num[1]*x1[-2] + Num[0]*x1[-1] \
						- Den[2]*TemporaryIbInput1[-2] - Den[1]*TemporaryIbInput1[-1])/Den[0])
					TemporaryIbInput2.append((Num[2]*x2[-3] + Num[1]*x2[-2] + Num[0]*x2[-1] \
						- Den[2]*TemporaryIbInput2[-2] - Den[1]*TemporaryIbInput2[-1])/Den[0])

				IbInput1.append(TemporaryIbInput1[-1]*(TemporaryIbInput1[-1]>0))
				IbInput2.append(TemporaryIbInput2[-1]*(TemporaryIbInput2[-1]>0))

				# if i == 0: dp.compare(float(IbInput[-1]),1.554682907070062)
				# if i == 2*SamplingRatio: dp.compare(float(IbInput[-1]),1.554682907070062)
				# if i == 3*SamplingRatio: dp.compare(float(IbInput[-1]),1.554682907070062)
				# if i == 4*SamplingRatio: dp.compare(float(IbInput[-1]),1.554682907070062)
				if i in range(IaAfferentDelayTimeStep1-1,IbAfferentDelayTimeStep1):
					Input1.append(IaInput1[i-IaAfferentDelayTimeStep1-1]/IaSpindleGain1 \
								+FeedforwardInput1[i]) #input to the muscle
				elif i in range(IbAfferentDelayTimeStep1, IIAfferentDelayTimeStep1):
					Input1.append(IaInput1[i-IaAfferentDelayTimeStep1-1]/IaSpindleGain1 \
								-IbInput1[i-IbAfferentDelayTimeStep1-1]/IbGTOGain1 \
								+FeedforwardInput1[i]) #input to the muscle
				elif i in range(IIAfferentDelayTimeStep1, CorticalDelayTimeStep1):
					Input1.append(IaInput1[i-IaAfferentDelayTimeStep1-1]/IaSpindleGain1 \
								+IIInput1[i-IIAfferentDelayTimeStep1-1]/IISpindleGain1 \
								-IbInput1[i-IbAfferentDelayTimeStep1-1]/IbGTOGain1 \
								+FeedforwardInput1[i]) #input to the muscle
				elif i > CorticalDelayTimeStep1-1:
					FeedbackInput1 = TransCorticalLoopConstant1*(TargetForceTrajectory1[i]-OutputForceTendon1[i-CorticalDelayTimeStep1-1])/MaximumContractileElementForce1 + FeedbackInput1  # feedback input through cortical pathway
					Input1.append(IaInput1[i-IaAfferentDelayTimeStep1]/IaSpindleGain1 \
								+IIInput1[i-IIAfferentDelayTimeStep1]/IISpindleGain1 \
								-IbInput1[i-IbAfferentDelayTimeStep1]/IbGTOGain1 \
								+FeedbackInput1	) #input to the muscle
				else:	
					Input1.append(FeedforwardInput1[i])

				if i in range(IaAfferentDelayTimeStep2-1,IbAfferentDelayTimeStep2):
					Input2.append(IaInput2[i-IaAfferentDelayTimeStep2-1]/IaSpindleGain2 \
								+FeedforwardInput2[i]) #input to the muscle
				elif i in range(IbAfferentDelayTimeStep2, IIAfferentDelayTimeStep2):
					Input2.append(IaInput2[i-IaAfferentDelayTimeStep2-1]/IaSpindleGain2 \
								-IbInput2[i-IbAfferentDelayTimeStep2-1]/IbGTOGain2 \
								+FeedforwardInput2[i]) #input to the muscle
				elif i in range(IIAfferentDelayTimeStep2, CorticalDelayTimeStep2):
					Input2.append(IaInput2[i-IaAfferentDelayTimeStep2-1]/IaSpindleGain2 \
								+IIInput2[i-IIAfferentDelayTimeStep2-1]/IISpindleGain2 \
								-IbInput2[i-IbAfferentDelayTimeStep2-1]/IbGTOGain2 \
								+FeedforwardInput2[i]) #input to the muscle
				elif i > CorticalDelayTimeStep1-1:
					FeedbackInput2 = TransCorticalLoopConstant2*(TargetForceTrajectory2[i]-OutputForceTendon2[i-CorticalDelayTimeStep2-1])/MaximumContractileElementForce2 + FeedbackInput2  # feedback input through cortical pathway
					Input2.append(IaInput2[i-IaAfferentDelayTimeStep2]/IaSpindleGain2 \
								+IIInput2[i-IIAfferentDelayTimeStep2]/IISpindleGain2 \
								-IbInput2[i-IbAfferentDelayTimeStep2]/IbGTOGain2 \
								+FeedbackInput2	) #input to the muscle
				else:	
					Input2.append(FeedforwardInput2[i])
				
				Count = SamplingRatio+Count
			else:
				IbInput1.append(IbInput1[-1])
				Input1.append(Input1[-1])
				IbInput2.append(IbInput2[-1])
				Input2.append(Input2[-1])

		elif FeedbackOption == 'cortical_fb_only': # STILL NEED TO DO THIS ONE!
			if i > CorticalDelayTimeStep-1:
				FeedbackInput = TransCorticalLoopConstant*(TargetForceTrajectory[i]-OutputForceTendon[i-CorticalDelayTimeStep-1])/MaximumContractileElementForce + FeedbackInput  # feedback input through cortical pathway
				Input.append(FeedbackInput) #input to the muscle
			else:
				Input.append(FeedforwardInput[i])

			IaInput.append(0)
			IIInput.append(0)
			IbInput.append(0)

		\ended here on 10/14/16

		# Feedforward only
		# Input[i] = input[i] #TargetForceTrajectory[i]/MaximumContractileElementForce
		## Noise + additional input
		if Input[-1] < 0:
			Input[-1] = 0
		elif Input[-1] > 1:
			Input[-1] = 1

		#dp.nan_test(Input[-1],'Input')

		#dp.nan_test(CorticalInput[-1],'CorticalInput')

		if i > 4:
			Noise.append(2*(random.random()-0.5)*(np.sqrt(0.01*Input[i])*np.sqrt(3)))
			FilteredNoise.append((BButtersCoefficients[4]*Noise[i-4] + BButtersCoefficients[3]*Noise[i-3] + BButtersCoefficients[2]*Noise[i-2] + BButtersCoefficients[1]*Noise[i-1] + BButtersCoefficients[0]*Noise[i] \
									- AButtersCoefficients[4]*FilteredNoise[i-4] - AButtersCoefficients[3]*FilteredNoise[i-3] - AButtersCoefficients[2]*FilteredNoise[i-2] - AButtersCoefficients[1]*FilteredNoise[i-1])/AButtersCoefficients[0])
			#dp.nan_test(Noise[i],'Noise')
			#dp.nan_test(FilteredNoise[-1],'FilteredNoise')
		else:
			Noise.append(0)
			FilteredNoise.append(0)
			#dp.nan_test(Noise[i],'Noise')
			#dp.nan_test(FilteredNoise[-1],'FilteredNoise')

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

		#dp.nan_test(LongInput[i],'LongInput')

		# add activation filter
		TU = 0.12*(LongInput[-1] < EffectiveMuscleActivation) + 0.03 # When true TU = 0.15, else TU = 0.03

		EffectiveMuscleActivationFirstDeriv = (LongInput[-1] - EffectiveMuscleActivation)/TU
		EffectiveMuscleActivation = EffectiveMuscleActivationFirstDeriv*(SamplingPeriod) + EffectiveMuscleActivation # effective neural drive

		#dp.nan_test(EffectiveMuscleActivationFirstDeriv,'EffectiveMuscleActivationFirstDeriv')
		#dp.nan_test(EffectiveMuscleActivation,'EffectiveMuscleActivation')

		# ActivationFrequency,Y,fint,feff_dot,feff = activation_frequency_slow(EffectiveMuscleActivation,ContractileElementLength,ContractileElementVelocity,Y,fint,feff_dot,feff,SamplingFrequency,ActivationFrequency) # not used

		# force-velocity relationship
		ForceVelocity = (ContractileElementVelocity <= 0)*concentric_force_velocity(ContractileElementLength,ContractileElementVelocity) \
							+ (ContractileElementVelocity > 0)*eccentric_force_velocity(ContractileElementLength,ContractileElementVelocity) 
		# force-length relationship
		ForceLength = force_length(ContractileElementLength)
		ContractileElementForce = ForceLength*ForceVelocity
		# viscous property
		ForceViscocity = MuscleViscosity * ContractileElementVelocity
		# passive element 1
		ForcePassive1 = parallel_elastic_element_force_1(ContractileElementLength/MaximumContractileElementLength,MaximumContractileElementLength)
		# passive element 2
		ForcePassive2 = parallel_elastic_element_force_2(ContractileElementLength)
		ForcePassive2 = (ForcePassive2 <= 0)*ForcePassive2

		# total force from contractile element
		ForceTotal = (EffectiveMuscleActivation*(ContractileElementForce + ForcePassive2) + ForcePassive1 + ForceViscocity)*MaximumContractileElementForce
		ForceTotal = ForceTotal*(ForceTotal>=0.0)

		#dp.nan_test(ForceTotal,'ForceTotal')

		# force from series elastic element
		SeriesElasticElementForce = normalized_series_elastic_element_force(SeriesElasticElementLength) * MaximumContractileElementForce

		#dp.nan_test(SeriesElasticElementForce,'SeriesElasticElementForce')

		# calculate muscle excursion acceleration based on the difference
		# between muscle force and tendon force
		if i < len(Time)-1:
			MuscleAcceleration.append((SeriesElasticElementForce*np.cos(PennationAngle) - ForceTotal*(np.cos(PennationAngle))**2)/(MuscleMass) \
				+ (MuscleVelocity[-1])**2*np.tan(PennationAngle)**2/(MuscleLength[-1]))
			# integrate acceleration to get velocity
			MuscleVelocity.append((MuscleAcceleration[-1]+ \
				MuscleAcceleration[-2])/2*(SamplingPeriod)+MuscleVelocity[-1])
			# integrate velocity to get length
			MuscleLength.append((MuscleVelocity[-1]+ \
				MuscleVelocity[-2])/2*(SamplingPeriod)+MuscleLength[-1])

			#dp.nan_test(MuscleAcceleration[i+1],'MuscleAcceleration')
			#dp.nan_test(MuscleVelocity[i+1],'MuscleVelocity')
			#dp.nan_test(MuscleLength[i+1],'MuscleLength')

			# normalize each variable to optimal muscle length or tendon legnth
			ContractileElementAcceleration = MuscleAcceleration[-1]/(OptimalLength/100)
			ContractileElementVelocity = MuscleVelocity[-1]/(OptimalLength/100)
			ContractileElementLength = MuscleLength[-1]/(OptimalLength/100)
			SeriesElasticElementLength = (InitialMusculoTendonLength - ContractileElementLength*OptimalLength*np.cos(PennationAngle))/OptimalTendonLength

			#dp.nan_test(ContractileElementAcceleration,'ContractileElementAcceleration')
			#dp.nan_test(ContractileElementVelocity,'ContractileElementVelocity')
			#dp.nan_test(ContractileElementLength,'ContractileElementLength')
			#dp.nan_test(SeriesElasticElementLength,'SeriesElasticElementLength')

			# store data
			OutputSeriesElasticElementLength.append(SeriesElasticElementLength)
			OutputContractileElementLength.append(ContractileElementLength)
			OutputContractileElementVelocity.append(ContractileElementVelocity)
			OutputContractileElementAcceleration.append(ContractileElementAcceleration)
			OutputActivationFrequency.append(ActivationFrequency)
			OutputEffectiveMuscleActivation.append(EffectiveMuscleActivation)

		OutputForceMuscle.append(ForceTotal)
		OutputForceTendon.append(SeriesElasticElementForce)
		OutputForceLength.append(ForceLength)
		OutputForceVelocity.append(ForceVelocity)
		OutputForcePassive1.append(ForcePassive1)
		OutputForcePassive2.append(ForcePassive2)
		dp.statusbar(i,len(Time),StartTime=StartTime,Title = 'afferented_muscle_model')

	plt.figure()
	plt.plot(Time,OutputForceTendon)
	plt.plot(Time,TargetForceTrajectory,'r')
	plt.legend(['Output Force','Target Force'])
	plt.xlabel('Time (sec)')
	plt.ylabel('Force (N)')

	# save data as output in structure format
	output = {	'Force' : OutputForceMuscle, 										'ForceTendon' : OutputForceTendon, \
				'FL' : OutputForceLength, 											'FV' : OutputForceVelocity, \
				'PE Force 1' : OutputForcePassive1, 								'PE Force 2' : OutputForcePassive2, \
				'Target' : TargetForceTrajectory, 									'ContractileElementLength' : OutputContractileElementLength, \
				'ContractileElementVelocity' : OutputContractileElementVelocity, 	'ContractileElementAcceleration' : OutputContractileElementAcceleration, \
				'SeriesElasticElementLength' : OutputSeriesElasticElementLength, 	'Activation Frequency' : OutputActivationFrequency, \
				'Input' : Input, 													'Noise' : Noise, \
				'FilteredNoise' : FilteredNoise, 									'U' : OutputEffectiveMuscleActivation, \
				'Ia Input' : IaInput,												'II Input' : IIInput, \
				'Ib Input' : IbInput 	} #, 'CorticalInput' : CorticalInput }

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
NumberOfTrials = 10
Output = []
PXX = []
Range = []
for i in range(NumberOfTrials): #trialN = 1:10
	#rng('shufle')
	output = afferented_muscle_model(muscle_parameters,delay_parameters,gain_parameters,TargetTrajectory,CorticalInput,FeedbackOption = FeedbackOption)    
	Output.append(output)
	f,pxx = welch(output['ForceTendon'][-int(10*SamplingFrequency)-1:]-np.average(output['ForceTendon'][int(-10*SamplingFrequency-1):]),window = gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5)),\
		noverlap = SamplingFrequency,nperseg = len(gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5))), fs = SamplingFrequency) 
	pxx = smooth(pxx,5)
	PXX.append(np.array(pxx,ndmin=2))
	Range.append(max(output['Ia Input'][-10*SamplingFrequency:-1])-min(output['Ia Input'][-10*SamplingFrequency-1:]))
	print("\n")

PXX = np.concatenate(PXX,axis=0)
plt.figure()
plt.plot(f[:131],np.average(PXX[:,:131],axis = 0))

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