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