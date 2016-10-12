import numpy as np
from scipy.signal import sawtooth,square,gaussian,welch
import matplotlib.pyplot as plt
import warnings
from math import isnan
import time

def statusbar(i,N,**kwargs):
    """
    i is the current iteration (must be an int) and N is the length of 
    the range (must be an int). i must also be in [0,N). 
    
    ~~~~~~~~~~~~~~
    **kwargs
    ~~~~~~~~~~~~~~
    
    StartTime should equal time.time() and should be defined before your
    loop to ensure that you get an accurate representation of elapsed time.

    Title should be a str that will be displayed before the statusbar. Title
    should be no longer than 25 characters.

    ~~~~~~~~~~~~~~

    NOTE: you should place a print('\n') after the loop to ensure you
    begin printing on the next line.

    """
    import time
    StartTime = kwargs.get("StartTime",False)
    Title = kwargs.get("Title",'')

    assert type(i)==int, "i must be an int"
    assert type(N)==int, "N must be an int"
    assert N>i, "N must be greater than i"
    assert N>0, "N must be a positive integer"
    assert i>=0, "i must not be negative (can be zero)"
    assert type(Title) == str, "Title should be a string"
    assert len(Title) <= 25, "Title should be less than 25 characters"
    if Title != '': Title = ' '*(25-len(Title)) + Title + ': '
    statusbar = Title +'[' + '\u25a0'*int((i+1)/(N/50)) + '\u25a1'*(50-int((i+1)/(N/50))) + '] '
    if StartTime != False:
        print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete, ' + '{0:1.1f}'.format(time.time() - StartTime) + 'sec        \r', end='')
    else:
        print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete           \r',end = '')
def nan_test(x,xname):
	assert type(xname)==str,'xname must be a string'
	if type(x)==np.ndarray or type(x)==list: 
		assert np.shape(x)[0]==1, 'x must be a 1 x N array/list'
		assert any([isnan(elem) for elem in x])==False, xname + " contains NaN Values!"
	else:
		assert isnan(x)==False, xname + " is NaN!"
def smooth(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return(np.convolve(interval, window, 'same'))
def generate_target_force_trajectory(TrajectoryType,Time,Amplitude,AmplitudeModulation,FrequencyModulation):
	if TrajectoryType == 'triangle':
	    delay = (1/FrequencyModulation)/4
	    triangle = 2*sawtooth(2*np.pi*(Time-Delay),0.5)-1
	    TargetForceTrajectory = AmplitudeModulation*triangle+Amplitude
	elif TrajectoryType == 'sinewave':
	    TargetForceTrajectory = AmplitudeModulation*np.sin(2*np.pi*FrequencyModulation*Time)+Amplitude
	elif TrajectoryType == 'constant':
		SamplingFrequency = 10000
		TargetForceTrajectory = np.concatenate((np.zeros(SamplingFrequency), (Amplitude)*1/2*Time[:2*SamplingFrequency], Amplitude*np.ones(len(Time)-3*SamplingFrequency), Amplitude*np.ones(5000)))
		TargetForceTrajectory = smooth(TargetForceTrajectory,500)
		TargetForceTrajectory = TargetForceTrajectory[:-5000]
	    #TargetForceTrajectory = [Amplitude*np.ones(1,len(Time))]
	"""
	elif strcmp(type,'trapezoid')
		SamplingFrequency = 10000
		offset = 0.01 % in %MVC
		slope = Amplitude/5 % %MVC/s
		%Amplitude = 30 % in %MVC
		holdtime = 0
		interval = 5 % in seconds
		ramp_up = offset/interval*(Time(1:interval*SamplingFrequency))
		delay_rise_start = SamplingFrequency*interval
		delay_rise_end = SamplingFrequency*(Amplitude/slope)+delay_rise_start
		delay_decay_start = delay_rise_end + SamplingFrequency*holdtime
		delay_decay_end = delay_rise_end + SamplingFrequency*holdtime + SamplingFrequency*(Amplitude/slope)
		y1 = slope*(Time(delay_rise_start:delay_rise_end)-interval)
		y2 = -slope*(Time(delay_decay_start:delay_decay_end)-delay_decay_end/SamplingFrequency)
		%zeros(1,interval*SamplingFrequency)
		Trapezoidal_Pulse = [zeros(1,interval*SamplingFrequency) y1 ...
			Amplitude*ones(1,holdtime*SamplingFrequency) y2 zeros(1,interval*SamplingFrequency)]+offset
		N = floor(Time(end)/(len(Trapezoidal_Pulse)/SamplingFrequency))
		TargetForceTrajectory = repmat(Trapezoidal_Pulse,1,N)
		TargetForceTrajectory = [TargetForceTrajectory zeros(1,len(Time)-len(TargetForceTrajectory))+offset]
		TargetForceTrajectory(1:interval*SamplingFrequency) = ramp_up

		# 6%MVC => 0.0667
		# 3%MVC => 0.05
		# 1.5%MVC => 0.033
		# 3%MVC + 5s hold => 0.04
		# 3%MVC + 10s hold => 
	"""
	return(TargetForceTrajectory)
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
	control - proprioceptive system + supraspinal loop) and 'option_4'

	"""

	FeedbackOption = kwargs.get("FeedbackOption",'ff_only')
	assert FeedbackOption in ['ff_only','servo_control','fb_control','option_4'],\
		"FeedbackOption must be either 'ff_only', 'servo_control', 'fb_control', or 'option_4'"
	assert type(muscle_parameters)==dict, "muscle_parameters must be a dictionary"
	assert len(muscle_parameters)==6, "dict muscle_parameters can only have 6 entries"
	assert 'Pennation Angle' in muscle_parameters, "'Pennation Angle' missing in muscle_parameters"
	assert 'Muscle Mass' in muscle_parameters, "'Muscle Mass' missing in muscle_parameters"
	assert 'Optimal Length' in muscle_parameters, "'Optimal Length' missing in muscle_parameters"
	assert 'Tendon Length' in muscle_parameters, "'Tendon Length' missing in muscle_parameters"
	assert 'Initial Muscle Length' in muscle_parameters, "'Initial Muscle Length' missing in muscle_parameters"
	assert 'Initial Tendon Length' in muscle_parameters, "'Initial Tendon Length' missing in muscle_parameters"
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
	PassiveMuscleForce = MuscleConstant * MuscleRateConstant * np.log(np.exp((1 - RestingMuscleLength)/MuscleRateConstant)+1)
	NormalizedSeriesElasticLength = TendonRateConstant*np.log(np.exp(PassiveMuscleForce/TendonConstant/TendonRateConstant)+1)+NormalizedRestingTendonLength        
	SeriesElasticLength = TendonLength * NormalizedSeriesElasticLength

	MaximumMusculoTendonLength = OptimalLength * np.cos(PennationAngle) + TendonLength + 0.5
	NormalizedMaximumFascicleLength = (MaximumMusculoTendonLength - SeriesElasticLength)/OptimalLength
	MaximumContractileElementLength = NormalizedMaximumFascicleLength/np.cos(PennationAngle) # Contractile Element = Muscle

	#nan_test(PassiveMuscleForce,'PassiveMuscleForce')
	#nan_test(NormalizedSeriesElasticLength,'NormalizedSeriesElasticLength')
	#nan_test(MaximumMusculoTendonLength,'MaximumMusculoTendonLength')
	#nan_test(SeriesElasticLength,'SeriesElasticLength')
	#nan_test(NormalizedMaximumFascicleLength,'NormalizedMaximumFascicleLength')
	#nan_test(MaximumContractileElementLength,'MaximumContractileElementLength')

	InitialLength = (InitialMusculoTendonLength-(-OptimalTendonLength*(TendonRateConstant/MuscleRateConstant*RestingMuscleLength-NormalizedRestingTendonLength-TendonRateConstant*np.log(MuscleConstant/TendonConstant*MuscleRateConstant/TendonRateConstant))))/(100*(1+TendonRateConstant/MuscleRateConstant*OptimalTendonLength/MaximumContractileElementLength*1/OptimalLength)*np.cos(PennationAngle))
	ContractileElementLength = InitialLength/(OptimalLength/100) # initializing CE Length
	SeriesElasticElementLength = (InitialMusculoTendonLength - InitialLength*100)/OptimalTendonLength # MTU - initial muscle = initial SEE length

	#nan_test(InitialLength,'InitialLength')
	#nan_test(ContractileElementLength,'ContractileElementLength')
	#nan_test(SeriesElasticElementLength,'SeriesElasticElementLength')

	# calculate maximal force output of the muscle based on PCSA
	MuscleDensity = 1.06 #g/cm**3
	PCSA = (MuscleMass*1000)/MuscleDensity/OptimalLength #physionp.logical cross-sectional area
	SpecificTension = 31.4 # WHERE DOES THIS COME FROM AND WHAT DOES IT MEAN? N/
	MaximumContractileElementForce = PCSA * SpecificTension

	#nan_test(PCSA,'PCSA')

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
	Time = np.arange(0,len(TargetTrajectory)/SamplingFrequency,1/SamplingFrequency) # 0:1/SamplingFrequency:(length(TargetTrajectory)-1)/SamplingFrequency

	# parameter initialization
	empty = np.zeros(len(Time))
	ContractileElementVelocity = 0
	ContractileElementAcceleration = 0
	MuscleAcceleration = empty
	MuscleVelocity = empty
	MuscleLength = np.concatenate((np.array([ContractileElementLength*OptimalLength/100]),np.zeros(len(Time)-1)))
	SeriesElasticElementForce = 0.105

	# filter parameters for Noise
	BButtersCoefficients,AButtersCoefficients = signal.butter(4,100/(SamplingFrequency/2),'low')

	## Activation dynamics parameters
	Y_dot = 0
	Y = 0
	Saf_dot = 0
	Saf = 0
	fint_dot = 0
	fint = 0
	feff_dot = 0
	feff = 0

	ActivationFrequency = 0
	EffectiveNeuralDrive = 0

	DynamicSpindleForce = 0
	StaticSpindleForce = 0
	Bag1TensionSecondDeriv = 0
	Bag1TensionFirstDeriv = 0
	Bag1Tension = 0
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
	TransferFunction = signal.TransferFunction(Num,Den)
	#Hd = signal.TransferFunction(Num,Den,dt = 1/FeedbackSamplingFrequency) #c2d(TransferFunction,1/FeedbackSamplingFrequency)

	## Transcortical loop
	TransCorticalLoopConstant = 0.01

	## Gains
	GammaDynamicGain = gain_parameters['Gamma Dynamic Gain']
	GammaStaticGain = gain_parameters['Gamma Static Gain']
	IaSpindleGain = gain_parameters['Ia Gain']
	IISpindleGain = gain_parameters['II Gain']
	IbGTOGain = gain_parameters['Ib Gain']

	# Convert delays in ms to samples
	IaAfferentDelayTimeStep = IaAfferentDelay*SamplingFrequency/FeedbackSamplingFrequency   #Ia + II 30
	IIAfferentDelayTimeStep = IIAfferentDelay*SamplingFrequency/FeedbackSamplingFrequency
	IbAfferentDelayTimeStep = IbAfferentDelay*SamplingFrequency/FeedbackSamplingFrequency   #Ib 40
	CorticalDelayTimeStep = CorticalDelay*SamplingFrequency/FeedbackSamplingFrequency   #cortical 50

	FeedbackInput = 0
	Count = 0

	# Convert force trajectory to unit of newton
	TargetForceTrajectory = TargetTrajectory*MaximumContractileElementForce

	def bag1_model(Length,LengthFirstDeriv,LengthSecondDeriv,DynamicSpindleForce,Bag1TensionFirstDeriv,Bag1Tension):
		## Feedback system parameters
		## Spindle Model
		p = 2

		Bag1TimeConstant = 0.149
		Bag1Frequency = 60

		Bag1BetaNot = 0.0605
		Bag1Beta = 0.2592
		Bag1Gamma = 0.0289

		G = 20000  #7000

		if LengthFirstDeriv >= 0:
		    C = 1
		else:
		    C = 0.42       
		
		DynamicSpindleForceFirstDeriv = (GammaDynamicGain**p/(GammaDynamicGain**p+Bag1Frequency**p) - DynamicSpindleForce)\
											/ Bag1TimeConstant
		DynamicSpindleForce = (1/SamplingFrequency)*DynamicSpindleForceFirstDeriv + DynamicSpindleForce

		#nan_test(DynamicSpindleForceFirstDeriv,'DynamicSpindleForceFirstDeriv')
		#nan_test(DynamicSpindleForce,'DynamicSpindleForce')

		Bag1Beta = Bag1BetaNot + Bag1Beta * DynamicSpindleForce
		Bag1Gamma = Bag1Gamma * DynamicSpindleForce

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
		Bag1TensionFirstDeriv = Bag1TensionSecondDeriv*(1/SamplingFrequency) + Bag1TensionFirstDeriv
		Bag1Tension = Bag1TensionFirstDeriv*(1/SamplingFrequency) + Bag1Tension

		#nan_test(Bag1TensionSecondDeriv,'Bag1TensionSecondDeriv')
		#nan_test(Bag1TensionFirstDeriv,'Bag1TensionFirstDeriv')
		#nan_test(Bag1Tension,'Bag1Tension')

		Bag1AP = G*(Bag1Tension/K_SR-(LN_SR-L0_SR))

		#nan_test(Bag1AP,'Bag1AP')

		return(Bag1AP,DynamicSpindleForce,Bag1TensionFirstDeriv,Bag1Tension)
	def bag2_model(Length,LengthFirstDeriv,LengthSecondDeriv,StaticSpindleForce,Bag2TensionFirstDeriv,Bag2Tension):
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

		StaticSpindleForceFirstDeriv = (GammaStaticGain**p/(GammaStaticGain**p+Bag2Frequency**p)-StaticSpindleForce) \
											/ Bag2TimeConstant
		StaticSpindleForce = (1/SamplingFrequency)*StaticSpindleForceFirstDeriv + StaticSpindleForce

		#nan_test(StaticSpindleForceFirstDeriv,'StaticSpindleForceFirstDeriv')
		#nan_test(StaticSpindleForce,'StaticSpindleForce')

		Bag2Beta = Bag2BetaNot + Bag2Beta * StaticSpindleForce
		Bag2Gamma = Bag2Gamma * StaticSpindleForce

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
		Bag2TensionFirstDeriv = Bag2TensionSecondDeriv*(1/SamplingFrequency) + Bag2TensionFirstDeriv
		Bag2Tension = Bag2TensionFirstDeriv*(1/SamplingFrequency) + Bag2Tension

		#nan_test(Bag2TensionSecondDeriv,'Bag2TensionSecondDeriv')
		#nan_test(Bag2TensionFirstDeriv,'Bag2TensionFirstDeriv')
		#nan_test(Bag2Tension,'Bag2Tension')

		Bag2APPrimary = G*(Bag2Tension/K_SR-(LN_SR-L0_SR))
		Bag2APSecondary = G*(	X*(L_secondary/L0_SR)*(Bag2Tension/K_SR-(LN_SR-L0_SR)) \
								+(1-X)*(L_secondary/L0_PR)*(Length-Bag2Tension/K_SR-(L0_SR+LN_PR))	)

		#nan_test(Bag2APPrimary,'Bag2APPrimary')
		#nan_test(Bag2APSecondary,'Bag2APSecondary')

		return(Bag2APPrimary,Bag2APSecondary,StaticSpindleForce,Bag2TensionFirstDeriv,Bag2Tension)
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

		ChainStaticSpindleForce = GammaStaticGain**p/(GammaStaticGain**p+ChainFrequency**p)

		#nan_test(ChainStaticSpindleForce,'ChainStaticSpindleForce')

		ChainBeta = ChainBetaNot + ChainBeta * ChainStaticSpindleForce
		ChainGamma = ChainGamma * StaticSpindleForce

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
		ChainTensionFirstDeriv = ChainTensionSecondDeriv*1/SamplingFrequency + ChainTensionFirstDeriv
		ChainTension = ChainTensionFirstDeriv*1/SamplingFrequency + ChainTension

		#nan_test(ChainTensionSecondDeriv,'ChainTensionSecondDeriv')
		#nan_test(ChainTensionFirstDeriv,'ChainTensionFirstDeriv')
		#nan_test(ChainTension,'ChainTension')

		ChainAPPrimary = G*(ChainTension/K_SR-(LN_SR-L0_SR))
		ChainAPSecondary = G*(	X*(L_secondary/L0_SR)*(ChainTension/K_SR-(LN_SR-L0_SR)) \
									+ (1-X)*(L_secondary/L0_PR)*(Length-ChainTension/K_SR-(L0_SR+LN_PR))	)

		#nan_test(ChainAPPrimary,'ChainAPPrimary')
		#nan_test(ChainAPSecondary,'ChainAPSecondary')

		return(ChainAPPrimary,ChainAPSecondary,ChainTensionFirstDeriv,ChainTension)
	def spindle_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,\
		DynamicSpindleForce,Bag1TensionFirstDeriv,Bag1Tension,\
		StaticSpindleForce,Bag2TensionFirstDeriv,Bag2Tension,\
		ChainTensionFirstDeriv,ChainTension):

		S = 0.156

		Bag1AP,DynamicSpindleForce,Bag1TensionFirstDeriv,Bag1Tension = bag1_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,DynamicSpindleForce,Bag1TensionFirstDeriv,Bag1Tension)
		Bag2APPrimary,Bag2APSecondary,StaticSpindleForce,Bag2TensionFirstDeriv,Bag2Tension = bag2_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,StaticSpindleForce,Bag2TensionFirstDeriv,Bag2Tension)
		ChainAPPrimary,ChainAPSecondary,ChainTensionFirstDeriv,ChainTension = chain_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,ChainTensionFirstDeriv,ChainTension)

		if Bag1AP < 0: Bag1AP = 0 
		if Bag2APPrimary < 0: Bag2APPrimary = 0
		if ChainAPPrimary < 0: ChainAPPrimary = 0
		if Bag2APSecondary < 0: Bag2APSecondary = 0
		if ChainAPSecondary < 0: ChainAPSecondary = 0	        

		if Bag1AP > (Bag2APPrimary+ChainAPPrimary):
		    Larger = Bag1AP
		    Smaller = Bag2APPrimary+ChainAPPrimary
		else:
		    Larger = Bag2APPrimary+ChainAPPrimary
		    Smaller = Bag1AP

		PrimaryOutput = Larger + S * Smaller
		SecondaryOutput = Bag2APSecondary + ChainAPSecondary

		if PrimaryOutput < 0:
		    PrimaryOutput = 0
		elif PrimaryOutput > 100000:
		    PrimaryOutput = 100000
		if SecondaryOutput < 0:
		    SecondaryOutput = 0
		elif SecondaryOutput > 100000:
		    SecondaryOutput = 100000

		#nan_test(PrimaryOutput,'PrimaryOutput')
		#nan_test(SecondaryOutput,'SecondaryOutput')

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
		Y = Y_dot*1/SamplingFrequency + Y 

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
		fint = fint_dot*1/SamplingFrequency + fint
		feff_dot = (fint - feff)/Tf
		feff = feff_dot*1/SamplingFrequency + feff

		#nan_test(fint_dot,'fint_dot')
		#nan_test(fint,'fint')
		#nan_test(feff_dot,'feff_dot')
		#nan_test(feff,'feff')

		nf = nf0 + nf1*(1/Length-1)
		ActivationFrequencySlow = 1 - np.exp(-(Y*feff/(af*nf))**nf)

		#nan_test(nf,'nf')
		#nan_test(ActivationFrequencySlow,'ActivationFrequencySlow')

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
		Saf = Saf_dot*1/SamplingFrequency + Saf
		fint_dot = (fenv - fint)/Tf
		fint = fint_dot*1/SamplingFrequency + fint
		feff_dot = (fint - feff)/Tf
		feff = feff_dot*1/SamplingFrequency + feff
		nf = nf0 + nf1*(1/Length-1)
		ActivationFrequencyFast = 1 - np.exp(-(S_af*feff/(af*nf))**nf)

		#nan_test(Saf_dot,'Saf_dot')
		#nan_test(Saf,'Saf')
		#nan_test(fint_dot,'fint_dot')
		#nan_test(fint,'fint')
		#nan_test(feff_dot,'feff_dot')
		#nan_test(feff,'feff')
		#nan_test(nf,'nf')
		#nan_test(ActivationFrequencyFast,'ActivationFrequencyFast')

		return(ActivationFrequencyFast,fint,feff_dot,feff,Saf)
	def force_length(Length):
		beta = 2.3
		omega = 1.12
		rho = 1.62
		ForceLength = np.exp(-abs(float((Length**beta - 1)/omega)**rho))
		#nan_test(ForceLength,'ForceLength')
		return(ForceLength)
	def concentric_force_velocity(Length,V):
		Vmax = -7.88
		cv0 = 5.88
		cv1 = 0
		ConcentricForceVelocity = (Vmax - V)/(Vmax + (cv0 + cv1*Length)*V)
		#nan_test(ConcentricForceVelocity,'ConcentricForceVelocity')
		return(ConcentricForceVelocity)
	def eccentric_force_velocity(Length,V):
		av0 = -4.7
		av1 = 8.41
		av2 = -5.34
		bv = 0.35
		EccentricForceVelocity = (bv - (av0 + av1*Length + av2*Length**2)*V)/(bv+V)
		#nan_test(EccentricForceVelocity,'EccentricForceVelocity')
		return(EccentricForceVelocity)
	def parallel_elastic_element_force_1(Length,MaximumContractileElementLength):
		c1_pe1 = 23.0  #355, 67.1
		k1_pe1 = 0.046 #0.04, 0.056
		Lr1_pe1 = 1.17  #1.35, 1.41
		ParallelElasticElementForce1 = c1_pe1 * k1_pe1 * np.log(np.exp((Length/MaximumContractileElementLength - Lr1_pe1)/k1_pe1)+1)
		#nan_test(ParallelElasticElementForce1,'ParallelElasticElementForce1')
		return(ParallelElasticElementForce1)
	def parallel_elastic_element_force_2(Length):
		c2_pe2 = -0.02 #0.01  -0.1
		k2_pe2 = -21
		Lr2_pe2 = 0.70 #0.79 0.59
		ParallelElasticElementForce2 = c2_pe2*np.exp((k2_pe2*(Length-Lr2_pe2))-1)
		#nan_test(ParallelElasticElementForce2,'ParallelElasticElementForce2')
		return(ParallelElasticElementForce2)
	def normalized_series_elastic_element_force(LT):
		cT_se = 27.8
		kT_se = 0.0047
		LrT_se = 0.964
		NormalizedSeriesElasticElementForce = cT_se * kT_se * np.log(np.exp((LT - LrT_se)/kT_se)+1)
		#nan_test(NormalizedSeriesElasticElementForce,'NormalizedSeriesElasticElementForce')
		return(NormalizedSeriesElasticElementForce)
	
	Input, IaInput, IIInput = empty,empty,empty
	IbInput, x, TemporaryIbInput = empty,empty,empty
	CorticalInput, Noise, FilteredNoise = empty,empty,empty
	LongInput = empty
	OutputForceMuscle,OutputForceTendon,OutputForceLength = empty, empty, empty
	OutputForceVelocity,OutputForcePassive1,OutputForcePassive2 = empty, empty, empty
	OutputSeriesElasticElementLength,OutputContractileElementVelocity = np.concatenate(([(InitialMusculoTendonLength - InitialLength*100)/OptimalTendonLength],np.zeros(len(Time)-1))), empty
	OutputContractileElementLength = np.concatenate(([InitialLength/(OptimalLength/100)],np.zeros(len(Time)-1)))
	OutputContractileElementAcceleration,OutputActivationFrequency,OutputEffectiveNeuralDrive = empty, empty, empty
	StartTime = time.time()
	for i in range(len(Time)): #= 1:length(Time)
		if FeedbackOption == 'ff_only': # Feedforward input only
			FeedforwardInput = TargetForceTrajectory[i]/MaximumContractileElementForce # feedforward input
			Input[i] = FeedforwardInput
			IaInput[i] = 0
			IIInput[i] = 0
			IbInput[i] = 0
		elif FeedbackOption == 'servo_control': # Servo control (feedforward + spindle and GTO)
			PrimaryOutput,SecondaryOutput = spindle_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,\
													DynamicSpindleForce,Bag1TensionFirstDeriv,Bag1Tension,\
													StaticSpindleForce,Bag2TensionFirstDeriv,Bag2Tension,\
													ChainTensionFirstDeriv,ChainTension)
			IaInput[i] = PrimaryOutput
			IIInput[i] = SecondaryOutput
			if i == Count:
				# Get spindle primary and secondary afferent activity

				# Get Ib activity
				x[i] = GTOConstant1*np.log(SeriesElasticElementForce/GTOConstant2+1)
				if i == 0:
					TemporaryIbInput[i] = x[i]
				elif i == 1*SamplingFrequency/FeedbackSamplingFrequency:
					TemporaryIbInput[i] = x[i]
				elif i == 2*SamplingFrequency/FeedbackSamplingFrequency:
					TemporaryIbInput[i] = (Num[0]*x[i])/Den[0]
				elif i == 3*SamplingFrequency/FeedbackSamplingFrequency:
					TemporaryIbInput[i] = (Num[1]*x[i-SamplingFrequency/FeedbackSamplingFrequency] + Num[0]*x[i] \
						- Den[1]*TemporaryIbInput[i-SamplingFrequency/FeedbackSamplingFrequency])/Den[0]
				elif i >= 3*SamplingFrequency/FeedbackSamplingFrequency:
					TemporaryIbInput[i] = (Num[2]*x[i-2*SamplingFrequency/FeedbackSamplingFrequency] + Num[1]*x[i-SamplingFrequency/FeedbackSamplingFrequency] + Num[0]*x[i] \
						- Den[2]*TemporaryIbInput[i-2*SamplingFrequency/FeedbackSamplingFrequency] - Den[1]*TemporaryIbInput[i-SamplingFrequency/FeedbackSamplingFrequency])/Den[0]

				if TemporaryIbInput[i]>0: 
					IbInput[i] = TemporaryIbInput[i]
				else:
					IbInput[i] = 0
			
				if i > IaAfferentDelayTimeStep-1 and i <= IbAfferentDelayTimeStep-1:
					FeedforwardInput = TargetForceTrajectory[i]/MaximumContractileElementForce # feedforward input
			
					Input[i] = IaInput[i-IaAfferentDelayTimeStep-1]/IaSpindleGain \
						+FeedforwardInput #input to the muscle

				elif i > IbAfferentDelayTimeStep-1 and i <= IIAfferentDelayTimeStep-1:
					FeedforwardInput = TargetForceTrajectory[i]/MaximumContractileElementForce # feedforward input

					Input[i] = IaInput[i-IaAfferentDelayTimeStep-1]/IaSpindleGain \
						-IbInput[i-IbAfferentDelayTimeStep-1]/IbGTOGain \
						+FeedforwardInput #input to the muscle
			
				elif i > IIAfferentDelayTimeStep-1 and i <= CorticalDelayTimeStep-1:
					FeedforwardInput = TargetForceTrajectory[i]/MaximumContractileElementForce # feedforward input

					Input[i] = IaInput[i-IaAfferentDelayTimeStep-1]/IaSpindleGain \
						+IIInput[i-IIAfferentDelayTimeStep-1]/IISpindleGain \
						-IbInput[i-IbAfferentDelayTimeStep-1]/IbGTOGain \
						+FeedforwardInput #input to the muscle
			
				else:
					FeedforwardInput = TargetForceTrajectory[i]/MaximumContractileElementForce
					Input[i] = FeedforwardInput
			
				Count = SamplingFrequency/FeedbackSamplingFrequency+Count
			else:
			
				IbInput[i] = IbInput[i-1]
				TemporaryIbInput[i] = TemporaryIbInput[i-1]
				Input[i] = Input[i-1]
	
		elif FeedbackOption == 'fb_control': # Feedback control (proprioceptive systems + supraspinal loop)
			PrimaryOutput,SecondaryOutput = spindle_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,\
													DynamicSpindleForce,Bag1TensionFirstDeriv,Bag1Tension,\
													StaticSpindleForce,Bag2TensionFirstDeriv,Bag2Tension,\
													ChainTensionFirstDeriv,ChainTension)
			IaInput[i] = PrimaryOutput
			IIInput[i] = SecondaryOutput
			if i == Count:
				# Get spindle primary and secondary afferent activity
				
				# Get Ib activity
				x[i] = GTOConstant1*np.log((SeriesElasticElementForce/GTOConstant2+1))
				if i == 0:
					TemporaryIbInput[i] = x[i]
				elif i == 1*SamplingFrequency/FeedbackSamplingFrequency:
					TemporaryIbInput[i] = x[i]
				elif i == 2*SamplingFrequency/FeedbackSamplingFrequency:
					TemporaryIbInput[i] = (Num[0]*x[i])/Den[0]
				elif i == 3*SamplingFrequency/FeedbackSamplingFrequency:
					TemporaryIbInput[i] = (Num[1]*x[i-SamplingFrequency/FeedbackSamplingFrequency] + Num[0]*x[i] \
						- Den[1]*TemporaryIbInput[i-SamplingFrequency/FeedbackSamplingFrequency])/Den[0]
				elif i > 3*SamplingFrequency/FeedbackSamplingFrequency:
					TemporaryIbInput[i] = (Num[2]*x[i-2*SamplingFrequency/FeedbackSamplingFrequency] + Num[1]*x[i-SamplingFrequency/FeedbackSamplingFrequency] + Num[0]*x[i] \
						- Den[2]*TemporaryIbInput[i-2*SamplingFrequency/FeedbackSamplingFrequency] - Den[1]*TemporaryIbInput[i-SamplingFrequency/FeedbackSamplingFrequency])/Den[0]

				if TemporaryIbInput[i]>0: 
					IbInput[i] = TemporaryIbInput[i]
				else:
					IbInput[i] = 0
				
				if i > IaAfferentDelayTimeStep-1 and i <= IbAfferentDelayTimeStep-1:
					FeedforwardInput = TargetForceTrajectory[i]/MaximumContractileElementForce # feedforward input
					
					Input[i] = IaInput[i-IaAfferentDelayTimeStep-1]/IaSpindleGain \
						+FeedforwardInput #input to the muscle
				
				elif i > IbAfferentDelayTimeStep-1 and i <= IIAfferentDelayTimeStep-1:
					FeedforwardInput = TargetForceTrajectory[i]/MaximumContractileElementForce # feedforward input
					
					Input[i] = IaInput[i-IaAfferentDelayTimeStep-1]/IaSpindleGain \
						-IbInput[i-IbAfferentDelayTimeStep-1]/IbGTOGain \
						+FeedforwardInput #input to the muscle
				elif i > IIAfferentDelayTimeStep-1 and i <= CorticalDelayTimeStep-1:
					FeedforwardInput = TargetForceTrajectory[i]/MaximumContractileElementForce # feedforward input
					
					Input[i] = IaInput[i-IaAfferentDelayTimeStep-1]/IaSpindleGain \
						+IIInput[i-IIAfferentDelayTimeStep-1]/IISpindleGain \
						-IbInput[i-IbAfferentDelayTimeStep-1]/IbGTOGain \
						+FeedforwardInput #input to the muscle
				
				elif i > CorticalDelayTimeStep-1:
					#FeedforwardInput = TargetForceTrajectory[i]/MaximumContractileElementForce # feedforward input
					FeedbackInput = TransCorticalLoopConstant*(TargetForceTrajectory[i]-OutputForceTendon[i-CorticalDelayTimeStep-1])/MaximumContractileElementForce + FeedbackInput  # feedback input through cortical pathway
					
					Input[i] = IaInput[i-IaAfferentDelayTimeStep]/IaSpindleGain \
						+IIInput[i-IIAfferentDelayTimeStep]/IISpindleGain \
						-IbInput[i-IbAfferentDelayTimeStep]/IbGTOGain \
						+FeedbackInput #input to the muscle
				
				else:
					FeedforwardInput = TargetForceTrajectory[i]/MaximumContractileElementForce
					Input[i] = FeedforwardInput
				
				Count = SamplingFrequency/FeedbackSamplingFrequency+Count
			else:
				IbInput[i] = IbInput[i-1]
				TemporaryIbInput[i] = TemporaryIbInput[i-1]
				Input[i] = Input[i-1]

		elif FeedbackOption == 'option_4':
			if i > CorticalDelayTimeStep-1:
				#FeedforwardInput = TargetForceTrajectory[i]/MaximumContractileElementForce # feedforward input
				FeedbackInput = TransCorticalLoopConstant*(TargetForceTrajectory[i]-OutputForceTendon[i-CorticalDelayTimeStep-1])/MaximumContractileElementForce + FeedbackInput  # feedback input through cortical pathway
				Input[i] = FeedbackInput #input to the muscle
			else:
				FeedforwardInput = TargetForceTrajectory[i]/MaximumContractileElementForce
				Input[i] = FeedforwardInput

			IaInput[i] = 0
			IIInput[i] = 0
			IbInput[i] = 0

		## Feedforward only
		#Input[i] = input[i] #TargetForceTrajectory[i]/MaximumContractileElementForce
		##
		## Noise + additional input
		if Input[i] < 0:
			Input[i] = 0
		elif Input[i] > 1:
			Input[i] = 1

		#nan_test(Input[i],'Input')

		CorticalInput[i] = CorticalInput[i]

		#nan_test(CorticalInput[i],'CorticalInput')

		if i > 4:
			Noise[i] = 2*(random.random()-0.5)*(np.sqrt(0.01*Input[i])*np.sqrt(3))
			FilteredNoise[i] = (BButtersCoefficients[4]*Noise[i-4] + BButtersCoefficients[3]*Noise[i-3] + BButtersCoefficients[2]*Noise[i-2] + BButtersCoefficients[1]*Noise[i-1] + BButtersCoefficients[0]*Noise[i] \
				- AButtersCoefficients[4]*FilteredNoise[i-4] - AButtersCoefficients[3]*FilteredNoise[i-3] - AButtersCoefficients[2]*FilteredNoise[i-2] - AButtersCoefficients[1]*FilteredNoise[i-1])/AButtersCoefficients[0]
			#nan_test(Noise[i],'Noise')
			#nan_test(FilteredNoise[i],'FilteredNoise')
		else:
			Noise[i] =0
			FilteredNoise[i] = 0
			#nan_test(Noise[i],'Noise')
			#nan_test(FilteredNoise[i],'FilteredNoise')

		Input[i] = Input[i] + FilteredNoise[i] + CorticalInput[i]

		if Input[i] < 0:
			Input[i] = 0
		elif Input[i] > 1:
			Input[i] = 1

		# add delay along efferent pathway
		if i > EfferentDelay*SamplingFrequency/FeedbackSamplingFrequency:
			LongInput[i] = Input[i-EfferentDelay*SamplingFrequency/FeedbackSamplingFrequency]
		else:
			LongInput[i] = 0 

		#nan_test(LongInput[i],'LongInput')

		# add activation filter
		if LongInput[i] >= EffectiveNeuralDrive:
			TU = 0.03
		elif LongInput[i] < EffectiveNeuralDrive:
			TU = 0.15

		EffectiveNeuralDriveFirstDeriv = (LongInput[i] - EffectiveNeuralDrive)/TU
		EffectiveNeuralDrive = EffectiveNeuralDriveFirstDeriv*(1/SamplingFrequency) + EffectiveNeuralDrive # effective neural drive

		#nan_test(EffectiveNeuralDriveFirstDeriv,'EffectiveNeuralDriveFirstDeriv')
		#nan_test(EffectiveNeuralDrive,'EffectiveNeuralDrive')

		ActivationFrequency,Y,fint,feff_dot,feff = activation_frequency_slow(EffectiveNeuralDrive,ContractileElementLength,ContractileElementVelocity,Y,fint,feff_dot,feff,SamplingFrequency,ActivationFrequency) # not used

		# force-velocity relationship
		if ContractileElementVelocity <= 0:
			ForceVelocity = concentric_force_velocity(ContractileElementLength,ContractileElementVelocity)
		elif ContractileElementVelocity > 0:
			ForceVelocity = eccentric_force_velocity(ContractileElementLength,ContractileElementVelocity)

		# force-length relationship
		ForceLength = force_length(ContractileElementLength)
		ContractileElementForce = ForceLength*ForceVelocity
		# viscous property
		ForceViscocity = MuscleViscosity * ContractileElementVelocity
		# passive element 1
		ForcePassive1 = parallel_elastic_element_force_1(ContractileElementLength/MaximumContractileElementLength,MaximumContractileElementLength)
		# passive element 2
		ForcePassive2 = parallel_elastic_element_force_2(ContractileElementLength)
		if ForcePassive2 > 0:
			ForcePassive2 = 0

		# total force from contractile element
		ForceTotal = (EffectiveNeuralDrive*(ContractileElementForce + ForcePassive2) + ForcePassive1 + ForceViscocity)*MaximumContractileElementForce
		if ForceTotal < 0.0:
			ForceTotal = 0.0

		#nan_test(ForceTotal,'ForceTotal')

		# force from series elastic element
		SeriesElasticElementForce = normalized_series_elastic_element_force(SeriesElasticElementLength) * MaximumContractileElementForce

		#nan_test(SeriesElasticElementForce,'SeriesElasticElementForce')

		# calculate muscle excursion acceleration based on the difference
		# between muscle force and tendon force
		if i < len(Time)-1:
			MuscleAcceleration[i+1] = (SeriesElasticElementForce*np.cos(PennationAngle) - ForceTotal*(np.cos(PennationAngle))**2)/(MuscleMass) \
				+ (MuscleVelocity[i])**2*np.tan(PennationAngle)**2/(MuscleLength[i])
			# integrate acceleration to get velocity
			MuscleVelocity[i+1] = (MuscleAcceleration[i+1]+ \
				MuscleAcceleration[i])/2*1/SamplingFrequency+MuscleVelocity[i]
			# integrate velocity to get length
			MuscleLength[i+1] = (MuscleVelocity[i+1]+ \
				MuscleVelocity[i])/2*1/SamplingFrequency+MuscleLength[i]

			#nan_test(MuscleAcceleration[i+1],'MuscleAcceleration')
			#nan_test(MuscleVelocity[i+1],'MuscleVelocity')
			#nan_test(MuscleLength[i+1],'MuscleLength')

			# normalize each variable to optimal muscle length or tendon legnth
			ContractileElementAcceleration = MuscleAcceleration[i+1]/(OptimalLength/100)
			ContractileElementVelocity = MuscleVelocity[i+1]/(OptimalLength/100)
			ContractileElementLength = MuscleLength[i+1]/(OptimalLength/100)
			SeriesElasticElementLength = (InitialMusculoTendonLength - ContractileElementLength*OptimalLength*np.cos(PennationAngle))/OptimalTendonLength

			#nan_test(ContractileElementAcceleration,'ContractileElementAcceleration')
			#nan_test(ContractileElementVelocity,'ContractileElementVelocity')
			#nan_test(ContractileElementLength,'ContractileElementLength')
			#nan_test(SeriesElasticElementLength,'SeriesElasticElementLength')

			# store data
			OutputSeriesElasticElementLength[i+1] = SeriesElasticElementLength
			OutputContractileElementLength[i+1] = ContractileElementLength
			OutputContractileElementVelocity[i+1] = ContractileElementVelocity
			OutputContractileElementAcceleration[i+1] = ContractileElementAcceleration
			OutputActivationFrequency[i+1] = ActivationFrequency
			OutputEffectiveNeuralDrive[i+1] = EffectiveNeuralDrive

		OutputForceMuscle[i] = ForceTotal
		OutputForceTendon[i] = SeriesElasticElementForce
		OutputForceLength[i] = ForceLength
		OutputForceVelocity[i] = ForceVelocity
		OutputForcePassive1[i] = ForcePassive1
		OutputForcePassive2[i] = ForcePassive2
		#OptimalLength = OptimalLength_initial*(0.15*(1-act)+1)
		statusbar(i,len(Time),StartTime=StartTime,Title = 'afferented_muscle_model')

	plt.figure()
	plt.plot(Time,OutputForceTendon)
	plt.plot(Time,TargetForceTrajectory,'r')
	plt.legend(['Output Force','Target Force'])
	plt.xlabel('Time (sec)')
	plt.ylabel('Force (N)')

	# save data as output in structure format
	output = {	'Force' : OutputForceMuscle, 	'ForceTendon' : OutputForceTendon, \
				'FL' : OutputForceLength, 		'FV' : OutputForceVelocity, \
				'PE Force 1' : OutputForcePassive1, 	'PE Force 2' : OutputForcePassive2, \
				'Target' : TargetForceTrajectory, 			'ContractileElementLength' : OutputContractileElementLength, \
				'ContractileElementVelocity' : OutputContractileElementVelocity, 				'ContractileElementAcceleration' : OutputContractileElementAcceleration, \
				'SeriesElasticElementLength' : OutputSeriesElasticElementLength, 				'Activation Frequency' : OutputActivationFrequency, \
				'Input' : Input, \
				'Noise' : Noise, 				'FilteredNoise' : FilteredNoise, \
				'U' : OutputEffectiveNeuralDrive, 				'Ia Input' : IaInput, \
				'II Input' : IIInput, 					'Ib Input' : IbInput, \
				'Ib_temp' : TemporaryIbInput } #, 'CorticalInput' : CorticalInput }

	return(output)

muscle_parameters = {"Pennation Angle":5*np.pi/180, "Muscle Mass":0.075,\
						"Optimal Length":10.1, "Tendon Length":23.5,\
						 "Initial Muscle Length":10.1+0.41, "Initial Tendon Length":23.5+0.09}
# Define delays based on limb length and conduction velocity of each
# pathway
distanceMuscleSpinalCord = 0.8 # cm
conductionVelocity_efferent = 48.5 #m/s (ref. Elias et al. 2014) S-type = 44-51, FR = 51-52, FF 52-53
conductionVelocity_Ia = 64.5 #m/s (ref. Elias et al. 2014)
conductionVelocity_II = 32.5 #m/s (ref. Elias et al. 2014)
conductionVelocity_Ib = 59 #m/s (ref. Elias et al. 2014)
synaptic_delay = 2 #ms (ref. Kandel)

delay_parameters = {"Efferent Delay": round(distanceMuscleSpinalCord/conductionVelocity_efferent*1000), \
					"Ia Delay" : round(distanceMuscleSpinalCord/conductionVelocity_Ia*1000) + synaptic_delay, \
					"II Delay" : round(distanceMuscleSpinalCord/conductionVelocity_II*1000) + 2*synaptic_delay, \
					"Ib Delay" : round(distanceMuscleSpinalCord/conductionVelocity_Ib*1000) + 2*synaptic_delay, \
					"Cortical Delay" : 50 }

# Define gain parameters for each neural pathway
gain_parameters = {"Gamma Dynamic Gain" : 70, \
					"Gamma Static Gain" : 70, \
					"Ia Gain" : 800, \
					"II Gain" : 3000, \
					"Ib Gain" : 1000 }

# Define force trajectory that model needs to track
# You can create any target trajectory as you want, but it has to be
# sampled at 10000 Hz to be consistent with the sampling frequency of
# muscle model
SamplingFrequency = 10000
Time = np.arange(0,15,1/SamplingFrequency) #0:1/SamplingFrequency:15
targetType = 'constant' #'triangle','sinwave','trapezoid' (Still need to do Trap)
targetTrajectory = generate_target_force_trajectory(targetType,Time,0.3,0,0) # in unit of #MVC

# Define additional input to muscle (e.g., 20 Hz oscillatory input)
# amplitude = 1/10*targetTrajectory' # in unit of #MVC
# frequency = 20 # in unit of Hz
# CorticalInput = amplitude.*sin(2*pi*requency*Time)
CorticalInput = np.zeros(len(targetTrajectory))

FeedbackOption = 'fb_control'
if FeedbackOption == 'ff_only':
    FeedbackOptionString = 'FF'
    targetTrajectory = 1.028*targetTrajectory
elif FeedbackOption == 'servo_control':
    FeedbackOptionString = 'SC'
    targetTrajectory = 0.9285*targetTrajectory
elif FeedbackOption == 'fb_control':
    FeedbackOptionString = 'FB'
elif FeedbackOption == 'option_4':
    FeedbackOptionString = 'CC'

# Output from the muscle model
NumberOfTrials = 1
Output = []
PXX = []
Range = np.zeros(NumberOfTrials)
for i in range(NumberOfTrials): #trialN = 1:10
	#rng('shufle')
	output = afferented_muscle_model(muscle_parameters,delay_parameters,gain_parameters,targetTrajectory,CorticalInput,FeedbackOption = FeedbackOption)    
	Output = np.concatenate((Output,[output]))
	#f,pxx = welch(output['ForceTendon'][-10*SamplingFrequency-1:]-np.average(output['ForceTendon'][-10*SamplingFrequency-1:]),window = gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5)),\
	#							noverlap = SamplingFrequency,nperseg = len(gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5))), fs = SamplingFrequency) 
	#pxx = smooth(pxx,5)
	#PXX = np.concatenate((PXX,[pxx]))
	#Range[i] = max(output['Ia Gain'][-10*SamplingFrequency:-1])-min(output['Ia Gain'][-10*SamplingFrequency-1:])

#plt.figure()
#plt.plot(f[:131],np.average(pxx[:131,:],2))

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

