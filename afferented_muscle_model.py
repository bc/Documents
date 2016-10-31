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
	import numpy as np 
	from scipy import signal
	import control
	import random
	import matplotlib.pyplot as plt

	# test input values
	FeedbackOption = kwargs.get("FeedbackOption",'ff_only')
	test_input_values(muscle_parameters,delay_parameters,gain_parameters,FeedbackOption = FeedbackOption)

	# muscle architectural parameters
	PennationAngle = muscle_parameters['Pennation Angle']
	InitialMusculoTendonLength = muscle_parameters['Initial Muscle Length']*np.cos(PennationAngle)+muscle_parameters['Initial Tendon Length'] 

	# viscosity of muscle (ref: Elias et al. 2014)
	MuscleViscosity = 0.005 #0.001	

	SamplingFrequency = 10000
	FeedbackSamplingFrequency = 1000
	SamplingRatio = SamplingFrequency/FeedbackSamplingFrequency
	SamplingPeriod = 1/SamplingFrequency
	Time = np.arange(0,len(TargetTrajectory)*SamplingPeriod,SamplingPeriod) 

	# filter parameters for Input['Noise']
	BButtersCoefficients,AButtersCoefficients = signal.butter(4,100/(SamplingFrequency/2),'low')

	# discrete transfer function for GTO output
	Num,Den = [1.7,2.58,0.4],[1,2.2,0.4]
	ContinuousTransferFunction = control.tf(Num,Den)
	DiscreteTransferFunction = control.matlab.c2d(ContinuousTransferFunction,1/FeedbackSamplingFrequency)
	Num,Den = control.matlab.tfdata(DiscreteTransferFunction)
	Num,Den = Num[0][0],Den[0][0]

	# Convert delays in ms to samples
	EfferentDelay = delay_parameters['Efferent Delay']
	EfferentDelayTimeStep = int(EfferentDelay*SamplingRatio)

	# define all dictionaries needed for muscle(s)
	Bag1,Bag2,Chain,SlowTwitch,FastTwitch,CE,SEE,PEE,Muscle,Input = return_initial_values(muscle_parameters,gain_parameters,TargetTrajectory,CorticalInput)

	random.seed()
	StartTime = time.time()
	for i in range(len(Time)): 
		update_total_input_at_step_i_single_muscle(i,Input,CE,SEE,Bag1,Bag2,Chain,Num,Den,delay_parameters,gain_parameters,SamplingRatio,SamplingPeriod,FeedbackOption)
		
		#add noise (and cortical input) to input
		# random.seed(1)
		add_noise_to_input(i,Input,AButtersCoefficients,BButtersCoefficients)

		# add delay along efferent pathway
		if i > EfferentDelayTimeStep: #NEED TO MAKE THIS >= after next run!
			Input['Long'].append(Input['Total'][i-EfferentDelayTimeStep])
		else:
			Input['Long'].append(0)

		# apply activation filter
		Muscle['Effective Activation'].append(apply_activation_filter(Input,Muscle,SamplingPeriod))
		
		# SlowTwitch = activation_frequency_slow(CE,SlowTwitch,Muscle['Effective Activation'][-1],SamplingPeriod) # not used

		# update all kinematics and kinetics from Effective Muscle Activation at timestep i
		update_kinematics_and_kinetics(CE,PEE,SEE,Muscle,MuscleViscosity,PennationAngle,InitialMusculoTendonLength,SamplingPeriod)
		
		statusbar(i,len(Time),StartTime=StartTime,Title = 'afferented_muscle_model')

	plt.figure()
	plt.plot(Time,SEE['Tendon Force'][1:])
	plt.plot(Time,Input['Target Force Trajectory'],'r')
	plt.legend(['Output Force','Target Force'])
	plt.xlabel('Time (sec)')
	plt.ylabel('Force (N)')

	# save data as output in structure format
	output = {	'Force' : Muscle['Force'], 											'ForceTendon' : SEE['Tendon Force'][1:], \
				'FL' : CE['FL'], 												'FV' : CE['FV'], \
				'PE Force 1' : PEE['Passive Force 1'], 							'PE Force 2' : PEE['Passive Force 2'], \
				'Target' : Input['Target Force Trajectory'], 									'ContractileElementLength' : CE['Length'][:-1], \
				'ContractileElementVelocity' : CE['Velocity'][1:-1],				 	'ContractileElementAcceleration' : CE['Acceleration'][1:-1], \
				'SeriesElasticElementLength' : SEE['Tendon Length'][:-1],					 	'Activation Frequency' : [0.0]*149999, \
				'Input' : Input['Total'], 											'Noise' : Input['Noise'], \
				'FilteredNoise' : Input['FilteredNoise'], 							'U' : Muscle['Effective Activation'][1:-1], \
				'Ia Input' : Input['Ia'],											'II Input' : Input['II'], \
				'Ib Input' : Input['Ib'] 	} #, 'CorticalInput' : CorticalInput }

	"""
	NOTES:

	'ContractileElementVelocity','ContractileElementAcceleration', and 'ForceTendon' are from 1: because originally the lists were not initialized even though Vce and Ace are initialized to zero
	'EffectiveMuscleActivation' has the same problem

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
NumberOfTrials = 10
Output = []
PXX = []
Range = []
for i in range(NumberOfTrials): #trialN = 1:10
	#rng('shufle')
	output = afferented_muscle_model(muscle_parameters,delay_parameters,gain_parameters,TargetTrajectory,CorticalInput,FeedbackOption = FeedbackOption)    
	# compare_output(output)
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