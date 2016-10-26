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


def bag1_model(Length,LengthFirstDeriv,LengthSecondDeriv,GammaDynamicGain,DynamicSpindleFrequency,Bag1TensionFirstDeriv,Bag1Tension):
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

	#nan_test(DynamicSpindleFrequencyFirstDeriv,'DynamicSpindleFrequencyFirstDeriv')
	#nan_test(DynamicSpindleFrequency,'DynamicSpindleFrequency')

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

	#nan_test(Bag1TensionSecondDeriv,'Bag1TensionSecondDeriv')
	#nan_test(Bag1TensionFirstDeriv,'Bag1TensionFirstDeriv')
	#nan_test(Bag1Tension,'Bag1Tension')

	Bag1AfferentPotential = G*(Bag1Tension/K_SR-(LN_SR-L0_SR))

	#nan_test(Bag1AfferentPotential,'Bag1AfferentPotential')

	return(Bag1AfferentPotential,DynamicSpindleFrequency,Bag1TensionFirstDeriv,Bag1Tension)
def bag2_model(Length,LengthFirstDeriv,LengthSecondDeriv,GammaStaticGain,StaticSpindleFrequency,Bag2TensionFirstDeriv,Bag2Tension):
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

	#nan_test(Bag2TensionSecondDeriv,'Bag2TensionSecondDeriv')
	#nan_test(Bag2TensionFirstDeriv,'Bag2TensionFirstDeriv')
	#nan_test(Bag2Tension,'Bag2Tension')

	Bag2PrimaryAfferentPotential = G*(Bag2Tension/K_SR-(LN_SR-L0_SR))
	Bag2SecondaryAfferentPotential = G*(	X*(L_secondary/L0_SR)*(Bag2Tension/K_SR-(LN_SR-L0_SR)) \
						+(1-X)*(L_secondary/L0_PR)*(Length-Bag2Tension/K_SR-(L0_SR+LN_PR))	)

	#nan_test(Bag2PrimaryAfferentPotential,'Bag2PrimaryAfferentPotential')
	#nan_test(Bag2SecondaryAfferentPotential,'Bag2SecondaryAfferentPotential')

	return(Bag2PrimaryAfferentPotential,Bag2SecondaryAfferentPotential,StaticSpindleFrequency,Bag2TensionFirstDeriv,Bag2Tension)
def chain_model(Length,LengthFirstDeriv,LengthSecondDeriv,GammaStaticGain,ChainTensionFirstDeriv,ChainTension):
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

	#nan_test(ChainStaticSpindleFrequency,'ChainStaticSpindleFrequency')

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

	#nan_test(ChainTensionSecondDeriv,'ChainTensionSecondDeriv')
	#nan_test(ChainTensionFirstDeriv,'ChainTensionFirstDeriv')
	#nan_test(ChainTension,'ChainTension')

	ChainPrimaryAfferentPotential = G*(ChainTension/K_SR-(LN_SR-L0_SR))
	ChainSecondaryAfferentPotential = G*(	X*(L_secondary/L0_SR)*(ChainTension/K_SR-(LN_SR-L0_SR)) \
						+ (1-X)*(L_secondary/L0_PR)*(Length-ChainTension/K_SR-(L0_SR+LN_PR))	)

	#nan_test(ChainPrimaryAfferentPotential,'ChainPrimaryAfferentPotential')
	#nan_test(ChainSecondaryAfferentPotential,'ChainSecondaryAfferentPotential')

	return(ChainPrimaryAfferentPotential,ChainSecondaryAfferentPotential,ChainTensionFirstDeriv,ChainTension)
def spindle_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,\
	GammaDynamicGain,DynamicSpindleFrequency,Bag1TensionFirstDeriv,Bag1Tension,\
	GammaStaticGain,StaticSpindleFrequency,Bag2TensionFirstDeriv,Bag2Tension,\
	ChainTensionFirstDeriv,ChainTension):

	S = 0.156

	Bag1AfferentPotential,DynamicSpindleFrequency,Bag1TensionFirstDeriv,Bag1Tension = bag1_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,GammaDynamicGain,DynamicSpindleFrequency,Bag1TensionFirstDeriv,Bag1Tension)
	Bag2PrimaryAfferentPotential,Bag2SecondaryAfferentPotential,StaticSpindleFrequency,Bag2TensionFirstDeriv,Bag2Tension = bag2_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,GammaStaticGain,StaticSpindleFrequency,Bag2TensionFirstDeriv,Bag2Tension)
	ChainPrimaryAfferentPotential,ChainSecondaryAfferentPotential,ChainTensionFirstDeriv,ChainTension = chain_model(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,GammaStaticGain,ChainTensionFirstDeriv,ChainTension)

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

	#nan_test(PrimaryOutput,'PrimaryOutput')
	#nan_test(SecondaryOutput,'SecondaryOutput')

	return(PrimaryOutput,SecondaryOutput,\
		DynamicSpindleFrequency,Bag1TensionFirstDeriv,Bag1Tension,\
		StaticSpindleFrequency,Bag2TensionFirstDeriv,Bag2Tension,\
		ChainTensionFirstDeriv,ChainTension)
def activation_frequency_slow(Activation,Length,LengthFirstDeriv,Y,fint,feff_dot,feff,SamplingFrequency,ActivationFrequencySlow):
	import numpy as np

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
	ActivationFrequencySlow = 1 - np.exp(-(Y*feff/(af*nf))**nf)

	#nan_test(nf,'nf')
	#nan_test(ActivationFrequencySlow,'ActivationFrequencySlow')

	return(ActivationFrequencySlow,Y,fint,feff_dot,feff)
def activation_frequency_fast(Activation,Length,Saf,fint,feff_dot,feff,SamplingFrequency,ActivationFrequencyFast):
	import numpy as np

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
	import numpy as np
	beta = 2.3
	omega = 1.12
	rho = 1.62
	ForceLength = np.exp(-abs((Length**beta - 1)/omega)**rho)
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
	import numpy as np
	c1_pe1 = 23.0  #355, 67.1
	k1_pe1 = 0.046 #0.04, 0.056
	Lr1_pe1 = 1.17  #1.35, 1.41
	ParallelElasticElementForce1 = c1_pe1 * k1_pe1 * np.log(np.exp((Length/MaximumContractileElementLength - Lr1_pe1)/k1_pe1)+1)
	#nan_test(ParallelElasticElementForce1,'ParallelElasticElementForce1')
	return(ParallelElasticElementForce1)
def parallel_elastic_element_force_2(Length):
	import numpy as np
	c2_pe2 = -0.02 #0.01  -0.1
	k2_pe2 = -21
	Lr2_pe2 = 0.70 #0.79 0.59
	ParallelElasticElementForce2 = c2_pe2*np.exp((k2_pe2*(Length-Lr2_pe2))-1)
	#nan_test(ParallelElasticElementForce2,'ParallelElasticElementForce2')
	return(ParallelElasticElementForce2)
def normalized_series_elastic_element_force(LT):
	import numpy as np
	cT_se = 27.8
	kT_se = 0.0047
	LrT_se = 0.964
	NormalizedSeriesElasticElementForce = cT_se * kT_se * np.log(np.exp((LT - LrT_se)/kT_se)+1)
	#nan_test(NormalizedSeriesElasticElementForce,'NormalizedSeriesElasticElementForce')
	return(NormalizedSeriesElasticElementForce)
def return_initial_values(muscle_parameters,gain_parameters):
	import numpy as np
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

	PennationAngle = muscle_parameters['Pennation Angle']
	MuscleMass = muscle_parameters['Muscle Mass']
	OptimalLength = muscle_parameters['Optimal Length']
	TendonLength = muscle_parameters['Tendon Length']
	OptimalTendonLength = TendonLength*1.05
	InitialMuscleLength = muscle_parameters['Initial Muscle Length']
	InitialTendonLength = muscle_parameters['Initial Tendon Length']
	InitialMusculoTendonLength = InitialMuscleLength*np.cos(PennationAngle)+InitialTendonLength

	GammaDynamicGain = gain_parameters['Gamma Dynamic Gain']
	GammaStaticGain = gain_parameters['Gamma Static Gain']
	IaSpindleGain = gain_parameters['Ia Gain']
	IISpindleGain = gain_parameters['II Gain']
	IbGTOGain = gain_parameters['Ib Gain']
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

	CE = contractile_element_parameters(ContractileElementLength,ContractileElementVelocity,ContractileElementAcceleration,\
											MaximumContractileElementLength,[],[],MaximumContractileElementForce)
	Bag1 = bag_1_parameters(GammaDynamicGain,DynamicSpindleFrequency,Bag1Tension,Bag1TensionFirstDeriv)
	Bag2 = bag_2_parameters(GammaStaticGain,StaticSpindleFrequency,Bag2Tension,Bag2TensionFirstDeriv)
	Chain = chain_parameters(ChainTension,ChainTensionFirstDeriv)
	SlowTwitch = slow_twitch_parameters(Y,fint,feff_dot,feff,ActivationFrequency)
	FastTwitch = slow_twitch_parameters(Saf,fint,feff_dot,feff,ActivationFrequency)
	SEE = series_elastic_element_parameters(SeriesElasticElementLength,SeriesElasticElementForce)
	PEE = parallel_elastic_element_parameters([],[])
	Input = initialize_dictionary(['EffectiveMuscleActivation','Ia','II','Ib','TempIb1',\
									'TempIb2','Noise','FilteredNoise','Total','LongInput'],\
									[[EffectiveMuscleActivation],[],[],[],[],[],[],[],[],[]])
	# Removed 'Feedforward','TargetTrajectory','CorticalInput','Feedback',
	Muscle = initialize_dictionary(['Length','Velocity','Acceleration','Muscle Force','Tendon Force','Activation Frequency'],\
									[[MuscleLength],[MuscleVelocity],[MuscleAcceleration],[],[],[]])

	# Input, IaInput, IIInput = [],[],[]
	# IbInput, x, TemporaryIbInput = [],[],[]
	# Noise, FilteredNoise, LongInput = [],[],[]
	# OutputForceMuscle,OutputForceTendon,OutputForceLength = [], [], []
	# OutputForceVelocity,OutputForcePassive1,OutputForcePassive2 = [], [], []
	# OutputSeriesElasticElementLength,OutputContractileElementVelocity = [(InitialMusculoTendonLength - InitialLength*100)/OptimalTendonLength], []
	# OutputContractileElementLength = [InitialLength/(OptimalLength/100)]
	# OutputContractileElementAcceleration,OutputActivationFrequency,OutputEffectiveMuscleActivation = [], [], []
	return(Bag1,Bag2,Chain,SlowTwitch,FastTwitch,CE,SEE,PEE,Muscle,Input)
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

def test_input_values(muscle_parameters, delay_parameters, gain_parameters, FeedbackOption):
	assert FeedbackOption in ['ff_only','servo_control','fb_control','cortical_fb_only'],\
	"FeedbackOption must be either 'ff_only', 'servo_control', 'fb_control', or 'cortical_fb_only'"
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
