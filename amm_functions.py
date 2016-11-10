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
def is_within_bounds(y,ymin,ymax):
	"""
	This takes in a 1D array/list or scalar value and tests to see if the values are within the allowed bounds [ymin, ymax].
	Returns TRUE if all values are within bounds.

	"""
	import numpy as np 
	y = np.array(y,ndmin=2)
	assert ymin<=ymax, "ymin must be less than or equal to ymax"
	if np.shape(y) != (1,1): assert np.shape(y)[0]==1, "y must either be a scalar or a 1D list/array"
	assert sum([type(elem)!= str and type(elem)!= np.str_ for elem in y[0]])==len(y[0]), "y cannot have string values"
	result = max(y[0])<=ymax and min(y[0])>=ymin
	return(result)
def compare(Variable1,Variable2,epsilon = 1e-3):
	import numpy as np 
	assert np.shape(Variable1)==np.shape(Variable2), "Variables have different shapes. \nVariable 1: "+str(np.shape(Variable1))+"\nVariable 2: "+str(np.shape(Variable2))
	assert type(Variable1)==type(Variable2), "Variables have different types. \nVariable 1: %s\nVariable 2: %s" %(type(Variable1),type(Variable2))
	if type(Variable1)==list: Variable1,Variable2 = np.array(Variable1),np.array(Variable2)
	if type(Variable1)==np.ndarray: error = sum((Variable1-Variable2)**2)
	if type(Variable1)!=np.ndarray: error = abs(Variable1-Variable2)
	assert error < epsilon, "Error threshold of "+str(epsilon)+" has been reached. \nVariable 1: "+str(Variable1)+"\nVariable 2: "+str(Variable2)
def nan_test(x,xname = ""):
	"""
	Takes in a value x and asserts that it is real. Can be either a list, array, or scalar value.
	xname is defaulted to be empty, but can be changed to a string with the variable name.
	"""
	import numpy as np
	from math import isnan
	if xname != None: assert type(xname)==str,'xname must be a string'
	if type(x)==np.ndarray or type(x)==list: 
		assert np.shape(x)[0]==1, 'x must be a 1 x N array/list'
		assert any([isnan(elem) for elem in x])==False, xname + " contains NaN Values!"
	else:
		if xname == "": 
			assert isnan(x)==False, "Value is/contains NaN!"
		else:
			assert isnan(x)==False, xname + " is/contains NaN!"
def smooth(x, window_size):
	"""
	Takes in an array/list and a scalar window_size and returns a smoother
	array by means of a running average convolution.
	"""
	import numpy as np
	assert type(window_size)!=str, "window_size cannot be a string."
	x = np.array(x,ndmin=2)
	assert np.shape(x)[0]==1, "x must be a 1xN list or array"
	x = x[0]
	window= np.ones(int(window_size))/float(window_size)
	return(np.convolve(x, window, 'same'))
def generate_target_force_trajectory(TrajectoryType,Time,Amplitude,AmplitudeModulation,FrequencyModulation,Phase = 0):
	"""
	TrajectoryType must either be 'triangle','sinewave','constant','trapezoid'. Time must be a 1D array. Amplitude,
	AmplitudeModulation, and FrequencyModulation are scalar values. Phase is the delay variable defaulted to zero.

	"""
	import numpy as np
	from scipy.signal import sawtooth

	assert TrajectoryType in ['triangle','sinewave','constant','trapezoid'], "Trajectory type must either be 'triangle'," + \
		" 'sinewave', 'constant', or 'trapezoid'."
	assert Amplitude != str, "Amplitude must be a scalar"
	assert AmplitudeModulation != str, "AmplitudeModulation must be a scalar"
	assert FrequencyModulation != str, "FrequencyModulation must be a scalar"
	assert Phase != str, "Phase must be a scalar"

	if TrajectoryType == 'triangle':
		Delay = (1/FrequencyModulation)/4
		triangle = 2*sawtooth(2*np.pi*(Time-Delay),0.5)-1
		TargetForceTrajectory = AmplitudeModulation*triangle+Amplitude
	elif TrajectoryType == 'sinewave':
		TargetForceTrajectory = Amplitude*np.sin(FrequencyModulation*(2*np.pi*Time-Phase))+AmplitudeModulation
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
def current_pennation_angle(NormalizedMuscleLength,muscle_parameters):
	import numpy as np
	InitialPennationAngle = muscle_parameters['Pennation Angle']
	CurrentPennationAngle = np.arcsin(np.sin(InitialPennationAngle)/NormalizedMuscleLength)
	return(CurrentPennationAngle)
def return_initial_values(muscle_parameters,gain_parameters,TargetTrajectory,CorticalInput):
	import numpy as np
	def contractile_element_parameters(Length,Velocity,Acceleration,MaximumContractileElementLength,FL,FV,Force,MaximumContractileElementForce):
		CE = initialize_dictionary(['Length','Velocity','Acceleration','Maximum Length','FL','FV','Force','Maximum Force'],\
									[[Length],[Velocity],[Acceleration],MaximumContractileElementLength,FL,FV,Force,MaximumContractileElementForce])
		return(CE)
	def bag_1_parameters(GammaDynamicGain,DynamicSpindleFrequency,Bag1Tension,Bag1TensionFirstDeriv):
		Bag1 = initialize_dictionary(['GammaDynamicGain','DynamicSpindleFrequency','Tension','TensionFirstDeriv'],\
										[GammaDynamicGain,[DynamicSpindleFrequency],[Bag1Tension],[Bag1TensionFirstDeriv]])
		return(Bag1)
	def bag_2_parameters(GammaStaticGain,StaticSpindleFrequency,Bag2Tension,Bag2TensionFirstDeriv):
			Bag2 = initialize_dictionary(['GammaStaticGain','StaticSpindleFrequency','Tension','TensionFirstDeriv'],\
											[GammaStaticGain,[StaticSpindleFrequency],[Bag2Tension],[Bag2TensionFirstDeriv]])
			return(Bag2)
	def chain_parameters(ChainTension,ChainTensionFirstDeriv):
		Chain = initialize_dictionary(['Tension','TensionFirstDeriv'],\
										[[ChainTension], [ChainTensionFirstDeriv]])
		return(Chain)
	def slow_twitch_parameters(Y,fint,feff_dot,feff,ActivationFrequencySlow):
		SlowTwitch = initialize_dictionary(['Y','fint','feff_dot','feff','ActivationFrequency'],\
											[Y,fint,feff_dot,feff,[ActivationFrequencySlow]])
		return(SlowTwitch)
	def fast_twitch_parameters(Saf,fint,feff_dot,feff,ActivationFrequencyFast):
		FastTwitch = initialize_dictionary(['Saf','fint','feff_dot','feff','ActivationFrequency'],\
											[Saf,fint,feff_dot,feff,[ActivationFrequencyFast]])
		return(FastTwitch)
	def series_elastic_element_parameters(Length,Force,OptimalTendonLength):
		SEE = initialize_dictionary(['Tendon Length','Tendon Force','Optimal Length'],\
									[[Length], [Force], OptimalTendonLength])
		return(SEE)
	def parallel_elastic_element_parameters(ForcePassive1,ForcePassive2):
		PEE = initialize_dictionary(['Passive Force 1','Passive Force 2'],\
									[ForcePassive1,ForcePassive2])
		return(PEE)
	def muscle_values(EffectiveMuscleActivation,MuscleLength,MuscleVelocity,MuscleAcceleration,MuscleMass,OptimalLength,InitialPennationAngle):
		Muscle = initialize_dictionary(['Effective Activation','Force','Length','Velocity','Acceleration','Mass','Optimal Length','Pennation Angle'],\
									[[EffectiveMuscleActivation],[],[MuscleLength],[MuscleVelocity],[MuscleAcceleration],MuscleMass,OptimalLength,[InitialPennationAngle]])
		return(Muscle)

	PennationAngle = muscle_parameters['Pennation Angle']
	MuscleMass = muscle_parameters['Muscle Mass']
	OptimalLength = muscle_parameters['Optimal Length']
	TendonLength = muscle_parameters['Tendon Length']
	OptimalTendonLength = TendonLength*1.05
	InitialMuscleLength = muscle_parameters['Initial Muscle Length']
	InitialTendonLength = muscle_parameters['Initial Tendon Length']
	# adapted from Elias et al 2014 eq. 4
	# see eq 5 for PennationAngle(NormalizedMuscleLength)
	InitialMusculoTendonLength = InitialMuscleLength*np.cos(PennationAngle)+InitialTendonLength

	GammaDynamicGain = gain_parameters['Gamma Dynamic Gain']
	GammaStaticGain = gain_parameters['Gamma Static Gain']
	IaSpindleGain = gain_parameters['Ia Gain']
	IISpindleGain = gain_parameters['II Gain']
	IbGTOGain = gain_parameters['Ib Gain']

	# calculate maximal force output of the muscle based on PCSA
	MuscleDensity = 1.06 #g/cm**3
	PCSA = (MuscleMass*1000)/MuscleDensity/OptimalLength #physionp.logical cross-sectional area
	SpecificTension = 31.4 
	MaximumContractileElementForce = PCSA * SpecificTension

	# calculate initial length based on balance between stiffness of muscle and
	# tendon. Parameter values are taken from Elias 2014 -- Soleus mm
	# Additional Ref Elias 2014 [75-77,82]
	TendonStiffness = 27.8
	TendonCurvatureConstant = 0.0047
	TendonLengthAtStartOfLinearRegion = 0.964

	MuscleStiffness = 23.0  #355, 67.1
	MuscleCurvatureConstant = 0.046 #0.04, 0.056
	RestingMuscleLength = 1.17  #1.35, 1.41
	MaximumPassiveMuscleForce = float(MuscleStiffness * MuscleCurvatureConstant * np.log(np.exp((1 - RestingMuscleLength)/MuscleCurvatureConstant)+1))
	NormalizedSeriesElasticLength = float(TendonCurvatureConstant*np.log(np.exp(MaximumPassiveMuscleForce/(TendonStiffness*TendonCurvatureConstant))-1)+TendonLengthAtStartOfLinearRegion)
	# Note: this expression was originally float(TendonCurvatureConstant*np.log(np.exp(MaximumPassiveMuscleForce*MaximumContractileElementForce/TendonStiffness/TendonCurvatureConstant)-1)+TendonLengthAtStartOfLinearRegion)
	# BUT Stiffness (27.8 from Elias et al 2014) is in units of MaximumContractileElementForce/OptimalTendonLength's (unitless/unitless). Therefore, 
	# multiplying Stiffness by MaximumContractileElementForce and dividing by OptimalTendonLength will bring back the proper units,
	# and the equation simplifies to the version seen above. (Dan Hagen - 11/9/16)
	SeriesElasticLength = float(TendonLength * NormalizedSeriesElasticLength)

	MaximumMusculoTendonLength = float(OptimalLength * np.cos(PennationAngle) + TendonLength + 1)
	NormalizedMaximumFascicleLength = float((MaximumMusculoTendonLength - SeriesElasticLength)/OptimalLength)
	MaximumContractileElementLength = float(NormalizedMaximumFascicleLength/np.cos(PennationAngle)) # Contractile Element = Muscle

	InitialLength = float((InitialMusculoTendonLength+OptimalTendonLength\
							*((TendonCurvatureConstant/MuscleCurvatureConstant)*RestingMuscleLength-TendonLengthAtStartOfLinearRegion\
							- TendonCurvatureConstant*np.log((MuscleStiffness/TendonStiffness)*(MuscleCurvatureConstant/TendonCurvatureConstant))))\
							/(100*(1+TendonCurvatureConstant/MuscleCurvatureConstant*OptimalTendonLength/MaximumContractileElementLength*(1/OptimalLength))*np.cos(PennationAngle)))
	ContractileElementLength = float(InitialLength/(OptimalLength/100)) # initializing CE Length
	SeriesElasticElementLength = float((InitialMusculoTendonLength - InitialLength*100)/OptimalTendonLength) # MTU - initial muscle = initial SEE length

	# parameter initialization
	ContractileElementVelocity = 0
	ContractileElementAcceleration = 0
	MuscleAcceleration = 0
	MuscleVelocity = 0
	MuscleLength = ContractileElementLength*OptimalLength/100
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
	FL, FV, Force = [],[],0
	
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
	FeedbackInput = 0
	TargetForceTrajectory = TargetTrajectory*MaximumContractileElementForce
	FeedforwardInput = TargetForceTrajectory/MaximumContractileElementForce

	CE = contractile_element_parameters(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,\
											MaximumContractileElementLength,FL,FV,Force,MaximumContractileElementForce)
	Bag1 = bag_1_parameters(GammaDynamicGain,DynamicSpindleFrequency,Bag1Tension,Bag1TensionFirstDeriv)
	Bag2 = bag_2_parameters(GammaStaticGain,StaticSpindleFrequency,Bag2Tension,Bag2TensionFirstDeriv)
	Chain = chain_parameters(ChainTension,ChainTensionFirstDeriv)
	SlowTwitch = slow_twitch_parameters(Y,fint,feff_dot,feff,ActivationFrequency)
	FastTwitch = slow_twitch_parameters(Saf,fint,feff_dot,feff,ActivationFrequency)
	SEE = series_elastic_element_parameters(SeriesElasticElementLength,SeriesElasticElementForce,OptimalTendonLength)
	PEE = parallel_elastic_element_parameters([],[])
	Input = initialize_dictionary(['Ia','II','Ib','IbTemp1','IbTemp2','Total','Noise','FilteredNoise',\
									'Long','Feedback','Target Force Trajectory','Feedforward','Cortical',\
									'Ia Interneuron', 'Ib Interneuron'],\
									[[],[],[],[],[],[],[],[],\
									[],FeedbackInput,TargetForceTrajectory,FeedforwardInput,CorticalInput,\
									[],[]])
	# Removed 'Feedforward','TargetTrajectory','CorticalInput','Feedback',
	Muscle = muscle_values(EffectiveMuscleActivation,MuscleLength,MuscleVelocity,MuscleAcceleration,MuscleMass,OptimalLength,PennationAngle)
	return(Bag1,Bag2,Chain,SlowTwitch,FastTwitch,CE,SEE,PEE,Muscle,Input)
def bag1_model(CE,Bag1,SamplingPeriod):
	# Elias et al 2014 Ref # 84
	import numpy as np
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
def bag2_model(CE,Bag2,SamplingPeriod):
	# Elias et al 2014 Ref # 84
	import numpy as np
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
def chain_model(CE,Chain,Bag2,SamplingPeriod):
	# Elias et al 2014 Ref # 84
	import numpy as np
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
def spindle_model(CE,Bag1,Bag2,Chain,SamplingPeriod):
	# Elias et al 2014 Ref # 84
	S = 0.156
	Bag1AfferentPotential = bag1_model(CE,Bag1,SamplingPeriod)
	Bag2PrimaryAfferentPotential,Bag2SecondaryAfferentPotential = bag2_model(CE,Bag2,SamplingPeriod)
	ChainPrimaryAfferentPotential,ChainSecondaryAfferentPotential = chain_model(CE,Chain,Bag2,SamplingPeriod)

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
	PrimaryOutput = bound(PrimaryOutput,0,100000)
	SecondaryOutput = bound(SecondaryOutput,0,100000)

	return(PrimaryOutput,SecondaryOutput)
def activation_frequency_slow(CE,SlowTwitch,Activation,SamplingPeriod):
	import numpy as np
	Length, LengthFirstDeriv = CE['Length'][-1], CE['Velocity'][-1]
	Y,fint,feff_dot,feff = SlowTwitch['Y'],SlowTwitch['fint'],SlowTwitch['feff_dot'],SlowTwitch['feff']
	ActivationFrequencySlow = SlowTwitch['ActivationFrequency'][-1]

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
	SlowTwitch['ActivationFrequency'].append(1 - np.exp(-(Y*feff/(af*nf))**nf))
	SlowTwitch['Y'],SlowTwitch['fint'],SlowTwitch['feff_dot'],SlowTwitch['feff'] = Y,fint,feff_dot,feff 
	return(SlowTwitch)
def activation_frequency_fast(CE,FastTwitch,Activation,SamplingPeriod):
	import numpy as np
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
	FastTwitch['ActivationFrequency'].append(1 - np.exp(-(S_af*feff/(af*nf))**nf))
	FastTwitch['Saf'],FastTwitch['fint'],FastTwitch['feff_dot'],FastTwitch['feff'] = Saf,fint,feff_dot,feff
	return(FastTwitch)
def force_length(CE):
	import numpy as np
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
	import numpy as np
	Length, MaximumContractileElementLength = CE['Length'][-1], CE['Maximum Length']
	c1_pe1 = 23.0  #355, 67.1
	k1_pe1 = 0.046 #0.04, 0.056
	Lr1_pe1 = 1.17  #1.35, 1.41
	ParallelElasticElementForce1 = c1_pe1 * k1_pe1 * np.log(np.exp((Length/MaximumContractileElementLength - Lr1_pe1)/k1_pe1)+1)
	return(ParallelElasticElementForce1)
def parallel_elastic_element_force_2(CE):
	import numpy as np
	Length = CE['Length'][-1]
	c2_pe2 = -0.02 #0.01  -0.1
	k2_pe2 = -21
	Lr2_pe2 = 0.70 #0.79 0.59
	ParallelElasticElementForce2 = c2_pe2*np.exp((k2_pe2*(Length-Lr2_pe2))-1)
	return(ParallelElasticElementForce2)
def normalized_series_elastic_element_force(SEE):
	# adapted from Elias et al 2014 Equation 3
	# parameter values are for Soleus
	# Additional Ref Elias 2014 [75-77,82]
	import numpy as np
	TendonLength = SEE['Tendon Length'][-1]
	TendonStiffness = 27.8
	TendonCurvatureConstant = 0.0047
	TendonLengthAtStartOfLinearRegion = 0.964
	NormalizedSeriesElasticElementForce = TendonStiffness * TendonCurvatureConstant * np.log(np.exp((TendonLength - TendonLengthAtStartOfLinearRegion)/TendonCurvatureConstant)+1)
	return(NormalizedSeriesElasticElementForce)
def update_ib_input(i,Input,SEE,Num,Den):
	import numpy as np
	## GTO model
	# adjustable parameters chosen to reproduce individual Ib firing rates compatable with literature
	# Ref 59 in Elias et al 2014
	GTOConstant1 = 60 # Hz
	GTOConstant2 = 4 # N
	# Get Ib activity
	# Ref Elias et al 2014 eq. 10
	# 'IbTemp1' is in time domain
	Input['IbTemp1'].append(GTOConstant1*np.log(SEE['Tendon Force'][-1]/GTOConstant2+1))
	if i == 0 or i == 1:
		Input['IbTemp2'].append(Input['IbTemp1'][-1])
	elif i == 2:
		Input['IbTemp2'].append((Num[0]*Input['IbTemp1'][-1])/Den[0])
	elif i == 3:
		Input['IbTemp2'].append((Num[1]*Input['IbTemp1'][-2] + Num[0]*Input['IbTemp1'][-1]- Den[1]*Input['IbTemp2'][-1])/Den[0])
	else:
		Input['IbTemp2'].append((Num[2]*Input['IbTemp1'][-3] + Num[1]*Input['IbTemp1'][-2] + Num[0]*Input['IbTemp1'][-1] \
								- Den[2]*Input['IbTemp2'][-2] - Den[1]*Input['IbTemp2'][-1])/Den[0])
	Input['Ib'].append(Input['IbTemp2'][-1]*(Input['IbTemp2'][-1]>0))
def update_input_with_delayed_spindle_afferents_and_feedback(i,Input,SEE,CE,delay_parameters,gain_parameters,SamplingFrequency,**kwargs):
	ControlStrategy = kwargs.get("ControlStrategy",None)
	SEE_synergist,CE_synergist = kwargs.get("SynergistParameters",[None,None])
	## Transcortical loop
	TransCorticalLoopGain = 0.0005
	IaSpindleGain = gain_parameters['Ia Gain']
	IISpindleGain = gain_parameters['II Gain']

	IaAfferentDelay = delay_parameters['Ia Delay']
	IIAfferentDelay = delay_parameters['II Delay']
	CorticalDelay = delay_parameters['Cortical Delay']

	IaAfferentDelayTimeStep = IaAfferentDelay*SamplingFrequency/1000   #Ia + II 30
	IIAfferentDelayTimeStep = IIAfferentDelay*SamplingFrequency/1000
	CorticalDelayTimeStep = CorticalDelay*SamplingFrequency/1000   #cortical 50

	IaAfferentDelayTimeStep = int(IaAfferentDelayTimeStep)   #Ia + II 30
	IIAfferentDelayTimeStep = int(IIAfferentDelayTimeStep)
	CorticalDelayTimeStep = int(CorticalDelayTimeStep)   #cortical 50

	if i >=IaAfferentDelayTimeStep: DelayedIaInput = Input['Ia'][i-IaAfferentDelayTimeStep]
	if i >=IIAfferentDelayTimeStep: DelayedIIInput = Input['II'][i-IIAfferentDelayTimeStep]
	if i >=CorticalDelayTimeStep:
		if ControlStrategy == "synergist": 
			DelayedTendonForce = SEE['Tendon Force'][i-CorticalDelayTimeStep] + SEE_synergist['Tendon Force'][i-CorticalDelayTimeStep]
			MaximumForce = CE['Maximum Force'] + CE_synergist['Maximum Force']
		else:
			DelayedTendonForce = SEE['Tendon Force'][i-CorticalDelayTimeStep]
			MaximumForce = CE['Maximum Force']
	if i in range(IaAfferentDelayTimeStep,IIAfferentDelayTimeStep):
		Input['Total'].append(DelayedIaInput/IaSpindleGain\
							+ Input['Feedforward'][i]	) #input to the muscle
	elif i in range(IIAfferentDelayTimeStep, CorticalDelayTimeStep):
		Input['Total'].append(DelayedIaInput/IaSpindleGain \
							+ DelayedIIInput/IISpindleGain \
							+ Input['Feedforward'][i]	) #input to the muscle
	elif i >= CorticalDelayTimeStep:
		Input['Feedback'] = TransCorticalLoopGain*(Input['Target Force Trajectory'][i]-DelayedTendonForce)/MaximumForce + Input['Feedback']  # feedback input through cortical pathway
		Input['Total'].append(DelayedIaInput/IaSpindleGain \
							+DelayedIIInput/IISpindleGain \
							+Input['Feedback']	) #input to the muscle
	else:	
		Input['Total'].append(Input['Feedforward'][i])
def update_input_with_delayed_spindle_afferents_and_feedback_multiple_muscles(i,Input,SEE,CE,delay_parameters,gain_parameters,SamplingFrequency,**kwargs):
	ControlStrategy = kwargs.get("ControlStrategy",None)
	SEE_synergist,CE_synergist = kwargs.get("SynergistParameters",[None,None])
	## Transcortical loop
	TransCorticalLoopGain = 0.01
	IaSpindleGain = gain_parameters['Ia Gain']

	IaAfferentDelay = delay_parameters['Ia Delay']
	CorticalDelay = delay_parameters['Cortical Delay']

	IaAfferentDelayTimeStep = IaAfferentDelay*SamplingFrequency/1000   #Ia + II 30
	CorticalDelayTimeStep = CorticalDelay*SamplingFrequency/1000   #cortical 50

	IaAfferentDelayTimeStep = int(IaAfferentDelayTimeStep)   #Ia + II 30
	CorticalDelayTimeStep = int(CorticalDelayTimeStep)   #cortical 50

	if i >=IaAfferentDelayTimeStep: DelayedIaInput = Input['Ia'][i-IaAfferentDelayTimeStep]
	if i >=CorticalDelayTimeStep:
		if ControlStrategy == "synergist": 
			DelayedTendonForce = SEE['Tendon Force'][i-CorticalDelayTimeStep] + SEE_synergist['Tendon Force'][i-CorticalDelayTimeStep]
			MaximumForce = CE['Maximum Force'] + CE_synergist['Maximum Force']
		else:
			DelayedTendonForce = SEE['Tendon Force'][i-CorticalDelayTimeStep]
			MaximumForce = CE['Maximum Force']
	if i in range(IaAfferentDelayTimeStep,CorticalDelayTimeStep):
		Input['Total'].append(DelayedIaInput/IaSpindleGain\
							+ Input['Feedforward'][i]	) #input to the muscle
	elif i >= CorticalDelayTimeStep:
		Input['Feedback'] = TransCorticalLoopGain*(Input['Target Force Trajectory'][i]-DelayedTendonForce)/MaximumForce + Input['Feedback']  # feedback input through cortical pathway
		Input['Total'].append(DelayedIaInput/IaSpindleGain \
							+Input['Feedback']	) #input to the muscle
	else:	
		Input['Total'].append(Input['Feedforward'][i])
def update_input_with_delayed_GTO_afferent(i,Input,delay_parameters,gain_parameters,SamplingFrequency):
	IbGTOGain = gain_parameters['Ib Gain']	
	IbAfferentDelay = delay_parameters['Ib Delay']
	IbAfferentDelayTimeStep = IbAfferentDelay*SamplingFrequency/1000
	IbAfferentDelayTimeStep = int(IbAfferentDelayTimeStep)   #Ib 40
	if i >=IbAfferentDelayTimeStep: 
		DelayedIbInput = Input['Ib'][i-IbAfferentDelayTimeStep]
		Input['Total'][-1] -= DelayedIbInput/IbGTOGain
def update_input_with_delayed_spindle_afferents(i,Input,delay_parameters,gain_parameters,SamplingFrequency):
	## Transcortical loop
	IaSpindleGain = gain_parameters['Ia Gain']
	IISpindleGain = gain_parameters['II Gain']

	IaAfferentDelay = delay_parameters['Ia Delay']
	IIAfferentDelay = delay_parameters['II Delay']
			
	IaAfferentDelayTimeStep = IaAfferentDelay*SamplingFrequency/1000   #Ia + II 30
	IIAfferentDelayTimeStep = IIAfferentDelay*SamplingFrequency/1000
		
	IaAfferentDelayTimeStep = int(IaAfferentDelayTimeStep)   #Ia + II 30
	IIAfferentDelayTimeStep = int(IIAfferentDelayTimeStep)
	
	if i >=IaAfferentDelayTimeStep: DelayedIaInput = Input['Ia'][i-IaAfferentDelayTimeStep]
	if i >=IIAfferentDelayTimeStep: DelayedIIInput = Input['II'][i-IIAfferentDelayTimeStep]
	if i in range(IaAfferentDelayTimeStep,IIAfferentDelayTimeStep):
		Input['Total'].append(DelayedIaInput/IaSpindleGain\
							+ Input['Feedforward'][i]	) #input to the muscle
	elif i >= IIAfferentDelayTimeStep:
		Input['Total'].append(DelayedIaInput/IaSpindleGain \
							+ DelayedIIInput/IISpindleGain \
							+ Input['Feedforward'][i]	) #input to the muscle
	else:	
		Input['Total'].append(Input['Feedforward'][i])
def update_input_with_only_cortical_feedback(Input,SEE,CE,delay_parameters,SamplingFrequency):
	TransCorticalLoopGain = 0.0005
	CorticalDelay = delay_parameters['Cortical Delay']
	CorticalDelayTimeStep = CorticalDelay*SamplingFrequency/1000   #cortical 50
	CorticalDelayTimeStep = int(CorticalDelayTimeStep)   #cortical 50
	if i >=CorticalDelayTimeStep: DelayedTendonForce = SEE['Tendon Force'][i-CorticalDelayTimeStep+1]
	elif i >= CorticalDelayTimeStep:
		Input['Feedback'] = TransCorticalLoopGain*(Input['Target Force Trajectory'][i]-DelayedTendonForce)/CE['Maximum Force'] + Input['Feedback']  # feedback input through cortical pathway
		append_dictionary(Input,['Total','Ia','II','Ib'],[Input['Feedback'],0,0,0])
	else:	
		append_dictionary(Input,['Total','Ia','II','Ib'],[Input['Feedforward'][i],0,0,0])
def update_total_input_at_step_i_single_muscle(i,Input,CE,SEE,Bag1,Bag2,Chain,Num,Den,delay_parameters,gain_parameters,SamplingFrequency,FeedbackOption,**kwargs):
	ControlStrategy = kwargs.get("ControlStrategy",None)
	SynergistParameters = kwargs.get("SynergistParameters",[None,None])
	SamplingPeriod = 1/SamplingFrequency
	if FeedbackOption == 'ff_only': # Feedforward input only
		append_dictionary(Input,['Total','Ia','II','Ib'],[Input['Feedforward'][i],0,0,0])
	elif FeedbackOption == 'servo_control': # Servo control (feedforward + spindle and GTO)
		append_dictionary(Input,['Ia','II'],spindle_model(CE,Bag1,Bag2,Chain,SamplingPeriod))
		update_ib_input(i,Input,SEE,Num,Den)
		update_input_with_delayed_spindle_afferents(i,Input,delay_parameters,gain_parameters,SamplingFrequency)
		update_input_with_delayed_GTO_afferent(i,Input,delay_parameters,gain_parameters,SamplingFrequency)
	elif FeedbackOption == 'fb_control': # Feedback control (proprioceptive systems + supraspinal loop)
		append_dictionary(Input,['Ia','II'],spindle_model(CE,Bag1,Bag2,Chain,SamplingPeriod))
		update_ib_input(i,Input,SEE,Num,Den)
		update_input_with_delayed_spindle_afferents_and_feedback(i,Input,SEE,CE,delay_parameters,gain_parameters,\
			SamplingFrequency,ControlStrategy=ControlStrategy,SynergistParameters=SynergistParameters)
		update_input_with_delayed_GTO_afferent(i,Input,delay_parameters,gain_parameters,SamplingFrequency)
	elif FeedbackOption == 'cortical_fb_only':
		update_input_with_only_cortical_feedback(Input,SEE,CE,delay_parameters,SamplingFrequency)
def update_spindle_afferent_inputs_at_step_i_multiple_muscles(i,Input,CE,SEE,Bag1,Bag2,Chain,Num,Den,delay_parameters,gain_parameters,SamplingFrequency,FeedbackOption,**kwargs):
	ControlStrategy = kwargs.get("ControlStrategy",None)
	SynergistParameters = kwargs.get("SynergistParameters",[None,None])
	SamplingPeriod = 1/SamplingFrequency
	if FeedbackOption == 'ff_only': # Feedforward input only
		append_dictionary(Input,['Total','Ia','II','Ib'],[Input['Feedforward'][i],0,0,0])
	elif FeedbackOption == 'servo_control': # Servo control (feedforward + spindle and GTO)
		append_dictionary(Input,['Ia','II'],spindle_model(CE,Bag1,Bag2,Chain,SamplingPeriod))
		update_ib_input(i,Input,SEE,Num,Den)
		update_input_with_delayed_spindle_afferents(i,Input,delay_parameters,gain_parameters,SamplingFrequency)
	elif FeedbackOption == 'fb_control': # Feedback control (proprioceptive systems + supraspinal loop)
		append_dictionary(Input,['Ia','II'],spindle_model(CE,Bag1,Bag2,Chain,SamplingPeriod))
		update_ib_input(i,Input,SEE,Num,Den)
		update_input_with_delayed_spindle_afferents_and_feedback_multiple_muscles(i,Input,SEE,CE,delay_parameters,gain_parameters,\
			SamplingFrequency,ControlStrategy=ControlStrategy,SynergistParameters=SynergistParameters)
	elif FeedbackOption == 'cortical_fb_only':
		update_input_with_only_cortical_feedback(Input,SEE,CE,delay_parameters,SamplingFrequency)
def update_ia_interneurons_at_step_i(i,Input_1,Input_2,InterneuronType):
	"""
	InterneuronType must be either "excitatory" or "inhibitory"
	"""
	assert InterneuronType in ['excitatory', 'inhibitory'], "InterneuronType must be either 'excitatory' or 'inhibitory'"
	InterneuronDelay = 20
	ExcitatoryWeight = 0.5
	InhibitoryWeight = 0.5
	InterneuronWeight = (InterneuronType=='excitatory')*ExcitatoryWeight - (InterneuronType=='inhibitory')*InhibitoryWeight
	if i <= InterneuronDelay:
		Input_1['Ia Interneuron'].append(0)
		Input_2['Ia Interneuron'].append(0)
	else:
		Input_1['Ia Interneuron'].append(Input_1['Ia'][-1] + InterneuronWeight*Input_2['Ia Interneuron'][i-InterneuronDelay] )
		# Input_1['Ia Interneuron'][-1] = bound(Input_1['Ia Interneuron'][-1],0,1)
		Input_2['Ia Interneuron'].append(Input_2['Ia'][-1] + InterneuronWeight*Input_1['Ia Interneuron'][i-InterneuronDelay] )
		# Input_2['Ia Interneuron'][-1] = bound(Input_2['Ia Interneuron'][-1],0,1)
def update_ib_interneurons_at_step_i(i,Input_1,Input_2,InterneuronType):
	"""
	InterneuronType must be either "excitatory" or "inhibitory"
	"""
	assert InterneuronType in ['excitatory', 'inhibitory'], "InterneuronType must be either 'excitatory' or 'inhibitory'"
	InterneuronDelay = 20
	ExcitatoryWeight = 0.5
	InhibitoryWeight = 0.5
	InterneuronWeight = (InterneuronType=='excitatory')*ExcitatoryWeight - (InterneuronType=='inhibitory')*InhibitoryWeight
	if i <= InterneuronDelay:
		Input_1['Ib Interneuron'].append(0)
		Input_2['Ib Interneuron'].append(0)
	else:
		Input_1['Ib Interneuron'].append(Input_1['Ib'][-1] + InterneuronWeight*Input_2['Ib Interneuron'][i-InterneuronDelay] )
		# Input_1['Ib Interneuron'][-1] = bound(Input_1['Ib Interneuron'][-1],0,1)
		Input_2['Ib Interneuron'].append(Input_2['Ib'][-1] + InterneuronWeight*Input_1['Ib Interneuron'][i-InterneuronDelay] )
		# Input_2['Ib Interneuron'][-1] = bound(Input_2['Ib Interneuron'][-1],0,1)
def add_interneuron_inputs_to_total_at_step_i(i,Input_1,Input_2,gain_1_parameters,delay_1_parameters,gain_2_parameters,delay_2_parameters,IaInterneuronType,IbInterneuronType,SamplingFrequency):
	IaSpindleGain_1, IaSpindleGain_2 = gain_1_parameters['Ia Gain'], gain_2_parameters['Ia Gain']
	IaReciprocalSpindleGain_1, IaReciprocalSpindleGain_2 = gain_1_parameters['Ia Reciprocal Gain'], gain_2_parameters['Ia Reciprocal Gain']
	IbGTOGain_1, IbGTOGain_2 = gain_1_parameters['Ib Gain'], gain_2_parameters['Ib Gain']
	IaAfferentDelay_1, IaAfferentDelay_2 = delay_1_parameters['Ia Delay'], delay_2_parameters['Ia Delay']
	IbAfferentDelay_1, IbAfferentDelay_2 = delay_1_parameters['Ib Delay'], delay_2_parameters['Ib Delay']
	IaAfferentDelayTimeStep_1, IaAfferentDelayTimeStep_2 = IaAfferentDelay_1*SamplingFrequency/1000, IaAfferentDelay_2*SamplingFrequency/1000
	IbAfferentDelayTimeStep_1, IbAfferentDelayTimeStep_2 = IbAfferentDelay_1*SamplingFrequency/1000, IbAfferentDelay_2*SamplingFrequency/1000
	IaAfferentDelayTimeStep_1, IaAfferentDelayTimeStep_2 = int(IaAfferentDelayTimeStep_1), int(IaAfferentDelayTimeStep_2)
	IbAfferentDelayTimeStep_1, IbAfferentDelayTimeStep_2 = int(IbAfferentDelayTimeStep_1), int(IbAfferentDelayTimeStep_2)

	# update the input values for both interneurons (Ia and Ib)
	update_ia_interneurons_at_step_i(i,Input_1,Input_2,IaInterneuronType)
	update_ib_interneurons_at_step_i(i,Input_1,Input_2,IbInterneuronType)

	# create delay time steps
	if i >=IaAfferentDelayTimeStep_1: 
		DelayedIaInput_1 = Input_1['Ia'][i-IaAfferentDelayTimeStep_1]
		DelayedIaInterneuronInput_1 = Input_1['Ia Interneuron'][i-IaAfferentDelayTimeStep_1]
	else:
		DelayedIaInput_1, DelayedIaInterneuronInput_1 = 0,0
	if i >=IaAfferentDelayTimeStep_2: 
		DelayedIaInput_2 = Input_2['Ia'][i-IaAfferentDelayTimeStep_2]
		DelayedIaInterneuronInput_2 = Input_2['Ia Interneuron'][i-IaAfferentDelayTimeStep_2]
	else:
		DelayedIaInput_2, DelayedIaInterneuronInput_2 = 0,0
	if i >=IbAfferentDelayTimeStep_1: 
		DelayedIbInterneuronInput_1 = Input_1['Ib Interneuron'][i-IbAfferentDelayTimeStep_1]
		CorrectiveIbInput_1 = Input_1['Ib'][i-IbAfferentDelayTimeStep_1]
	else:
		DelayedIbInterneuronInput_1, CorrectiveIbInput_1 = 0, 0
	if i >=IbAfferentDelayTimeStep_2: 
		DelayedIbInterneuronInput_2 = Input_2['Ib Interneuron'][i-IbAfferentDelayTimeStep_2]
		CorrectiveIbInput_2 = Input_2['Ib'][i-IbAfferentDelayTimeStep_2]
	else:
		DelayedIbInterneuronInput_2, CorrectiveIbInput_2 = 0, 0

	# update Total input
	# NOTE: WILL THE SIGNS CHANGE AS A FUNCTION OF INTERNEURON TYPE?	
	Input_1['Total'][-1] = Input_1['Total'][-1] \
							- DelayedIaInterneuronInput_2/IaReciprocalSpindleGain_2 \
							- DelayedIbInterneuronInput_1/IbGTOGain_1 \
							+ DelayedIbInterneuronInput_2/IbGTOGain_2 \
							+ DelayedIaInput_2/(IaSpindleGain_2)
	Input_2['Total'][-1] = Input_2['Total'][-1] \
							- DelayedIaInterneuronInput_1/IaReciprocalSpindleGain_1 \
							- DelayedIbInterneuronInput_2/IbGTOGain_2 \
							+ DelayedIbInterneuronInput_1/IbGTOGain_1 \
							+ DelayedIaInput_1/(IaSpindleGain_1)
def bound(value,min_value,max_value):
	if value > max_value:
		result = max_value
	elif value < min_value:
		result = min_value
	else:
		result = value
	return(result)
def add_noise_to_input(i,Input,AButtersCoefficients,BButtersCoefficients):
	"""
	must import and seed random before running this function
	"""
	import numpy as np
	import random
	## Input['Noise'] + additional input
	Input['Total'][-1] = bound(Input['Total'][-1],0,1)
	if i > 4:
		Input['Noise'].append(2*(random.random()-0.5)*(np.sqrt(0.5*Input['Total'][i])*np.sqrt(3)))
		Input['FilteredNoise'].append((BButtersCoefficients[4]*Input['Noise'][i-4] + BButtersCoefficients[3]*Input['Noise'][i-3] \
										+ BButtersCoefficients[2]*Input['Noise'][i-2] + BButtersCoefficients[1]*Input['Noise'][i-1] \
										+ BButtersCoefficients[0]*Input['Noise'][i] \
										- AButtersCoefficients[4]*Input['FilteredNoise'][i-4] - AButtersCoefficients[3]*Input['FilteredNoise'][i-3] \
										- AButtersCoefficients[2]*Input['FilteredNoise'][i-2] - AButtersCoefficients[1]*Input['FilteredNoise'][i-1]) \
										/ AButtersCoefficients[0])
	else:
		append_dictionary(Input,['Noise','FilteredNoise'],[0,0])

	Input['Total'][-1] = Input['Total'][-1] + Input['FilteredNoise'][-1] + Input['Cortical'][i]
	Input['Total'][-1] = bound(Input['Total'][-1],0,1)
def apply_activation_filter(Input,Muscle,SamplingPeriod):
	# add activation filter
	TU = 0.12*(Input['Long'][-1] < Muscle['Effective Activation'][-1]) + 0.03 # When true TU = 0.15, else TU = 0.03

	EffectiveMuscleActivationFirstDeriv = (Input['Long'][-1] - Muscle['Effective Activation'][-1])/TU
	EffectiveMuscleActivation = EffectiveMuscleActivationFirstDeriv*(SamplingPeriod) + Muscle['Effective Activation'][-1] # effective neural drive
	return(EffectiveMuscleActivation)
def update_CE_FLV_relationship(CE):
	# force-velocity relationship
	if CE['Velocity'][-1] <= 0:
		CE['FV'].append(concentric_force_velocity(CE))
	else:
		CE['FV'].append(eccentric_force_velocity(CE))
	# force-length relationship
	CE['FL'].append(force_length(CE))
	CE['Force'] = CE['FL'][-1]*CE['FV'][-1]
def update_PEE_with_current_CE_length(PEE,CE):
	# passive element 1
	PEE['Passive Force 1'].append(parallel_elastic_element_force_1(CE))
	# passive element 2
	ForcePassive2 = parallel_elastic_element_force_2(CE)
	PEE['Passive Force 2'].append((ForcePassive2 <= 0)*ForcePassive2)
def update_muscle_kinematics_with_current_muscle_and_tendon_forces(Muscle,SEE,SamplingPeriod):
	import numpy as np
	# calculate muscle excursion acceleration based on the difference
	# between muscle force and tendon force
	# Ref Elias et al 2014 eq 7
	Muscle['Acceleration'].append((SEE['Tendon Force'][-1]*np.cos(Muscle['Pennation Angle'][-1]) - Muscle['Force'][-1]*(np.cos(Muscle['Pennation Angle'][-1]))**2)/(Muscle['Mass']) \
		+ (Muscle['Velocity'][-1])**2*np.tan(Muscle['Pennation Angle'][-1])**2/(Muscle['Length'][-1]))
	# integrate acceleration to get velocity
	Muscle['Velocity'].append((Muscle['Acceleration'][-1]+ \
		Muscle['Acceleration'][-2])/2*(SamplingPeriod)+Muscle['Velocity'][-1])
	# integrate velocity to get length
	Muscle['Length'].append((Muscle['Velocity'][-1]+ \
		Muscle['Velocity'][-2])/2*(SamplingPeriod)+Muscle['Length'][-1])
def update_CE_kinematics(CE,Muscle):
	# normalize each variable to optimal muscle length or tendon legnth
	CE['Acceleration'].append(Muscle['Acceleration'][-1]/(Muscle['Optimal Length']/100))
	CE['Velocity'].append(Muscle['Velocity'][-1]/(Muscle['Optimal Length']/100))
	CE['Length'].append(Muscle['Length'][-1]/(Muscle['Optimal Length']/100))
def update_kinematics_and_kinetics(CE,PEE,SEE,Muscle,MuscleViscosity,InitialMusculoTendonLength,SamplingPeriod,muscle_parameters):
	import numpy as np
	# update FLV relationship
	update_CE_FLV_relationship(CE)
	# update force due to viscous property
	ForceViscosity = MuscleViscosity * CE['Velocity'][-1]
	# update passive forces with current CE length
	update_PEE_with_current_CE_length(PEE,CE)
	# update total muscle force from contractile element, passive forces, and viscous forces (ONLY POSITIVE VALUES)
	ForceTotal = (Muscle['Effective Activation'][-1]*(CE['Force'] + PEE['Passive Force 2'][-1]) + PEE['Passive Force 1'][-1] + ForceViscosity)*CE['Maximum Force']
	Muscle['Force'].append(ForceTotal*(ForceTotal>=0.0))
	# update force from series elastic element
	SEE['Tendon Force'].append(normalized_series_elastic_element_force(SEE) * CE['Maximum Force'])
	# update muscle kinematics based on most recent muscle/tendon forces
	update_muscle_kinematics_with_current_muscle_and_tendon_forces(Muscle,SEE,SamplingPeriod)
	# update CE kinematics based on most recent muscle kinematics
	update_CE_kinematics(CE,Muscle)
	# update SEE length based on most recent CE length
	SEE['Tendon Length'].append((InitialMusculoTendonLength - CE['Length'][-1]*Muscle['Optimal Length']*np.cos(Muscle['Pennation Angle'][-1]))/SEE['Optimal Length'])
	# update Muscle['Pennation Angle']
	Muscle['Pennation Angle'].append(current_pennation_angle(CE['Length'][-1],muscle_parameters))
def initialize_dictionary(keys,values):
	assert len(keys)==len(values), "keys and values must have the same length."
	result = {}
	for i in range(len(keys)): result[keys[i]]=values[i]
	return(result)
def append_dictionary(dictionary,keys,values):
	assert len(keys)==len(values), "keys and values must have the same length."
	assert type(dictionary)==dict, "Input must be a dictionary"
	for i in range(len(keys)): 
		if keys[i] not in dictionary.keys(): 
			dictionary[keys[i]]=values[i]
		else:
			dictionary[keys[i]].append(values[i])
def test_input_values(muscle_parameters, delay_parameters, gain_parameters, **kwargs):
	FeedbackOption = kwargs.get("FeedbackOption",'ff_only')
	ControlStrategy = kwargs.get("ControlStrategy",'synergist')
	assert FeedbackOption in ['ff_only','servo_control','fb_control','cortical_fb_only'],\
	"FeedbackOption must be either 'ff_only', 'servo_control', 'fb_control', or 'cortical_fb_only'"
	assert ControlStrategy in ['synergist', 'antagonist'], " ControlStrategy must be either 'synergist' or 'antagonist'"
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
	assert len(gain_parameters)==5 or len(gain_parameters)==6, "dict gain_parameters can only have 5 or 6 entries"
	assert 'Gamma Dynamic Gain' in gain_parameters, "'Gamma Dynamic Gain' missing in gain_parameters"
	assert 'Gamma Static Gain' in gain_parameters, "'Gamma Static Gain' missing in gain_parameters"
	assert 'Ia Gain' in gain_parameters, "'Ia Gain' missing in gain_parameters"
	assert 'II Gain' in gain_parameters, "'II Gain' missing in gain_parameters"
	assert 'Ib Gain' in gain_parameters, "'Ib Gain' missing in gain_parameters"
	if len(gain_parameters)==6: assert 'Ia Reciprocal Gain' in gain_parameters, "'Ia Reciprocal Gain' missing when incorporating interneurons"
def compare_output(TheoreticalOutput):
	import pickle
	import numpy as np 
	ActualOutput = pickle.load(open('OutputWithoutError.pkl','rb'))
	Keys = ActualOutput.keys()
	for key in Keys:
		assert type(ActualOutput[key])==type(TheoreticalOutput[key]), "output["+key+"] has changed!"
		if type(ActualOutput[key])==list or type(ActualOutput[key])== int or type(ActualOutput[key])== float:
			assert ActualOutput[key]==TheoreticalOutput[key], "output["+key+"] has changed!"
		elif type(ActualOutput[key])==np.ndarray:
			assert np.array(ActualOutput[key]==TheoreticalOutput[key],ndmin=2).all(1)[0], "output["+key+"] has changed!"