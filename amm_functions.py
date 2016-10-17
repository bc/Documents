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