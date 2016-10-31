import numpy as np
from scipy.signal import sawtooth,square,gaussian,welch
import matplotlib.pyplot as plt
import time
import ipdb
from amm_functions import *

def afferented_2_muscles_model(muscle_1_parameters,muscle_2_parameters,\
	delay_1_parameters,delay_2_parameters,\
	gain_1_parameters,gain_2_parameters,\
	TargetTrajectory_1,TargetTrajectory2,\
	CorticalInput_1,CorticalInput2,\
	**kwargs):
	"""
	muscle_1/2_parameters must be a dictionary with "Pennation Angle", "Muscle Mass",
	"Optimal Length", "Tendon Length", "Initial Muscle Length", and "Initial Tendon Length"
	values.

	delay_1/2_parameters must be a dictionary with "Efferent", "Ia","II", "Ib", and 
	"Cortical" values.

	gain_1/2_parameters must be a dictionary with "Gamma Dynamic Gain", "Gamma Static Gain",
	"Ia", "II", and "Ib" values.

	~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~

	FeedbackOption is a string that must be either 'ff_only' (feedforward only)
	'servo_control' (feedforward, spindle, and GTO), 'fb_control' (feedback
	control - proprioceptive system + supraspinal loop) and 'cortical_fb_only'

	ControlStrategy is a string that must be either "synergist" or "antagonist"
	"""
	import numpy as np 
	from scipy import signal
	import control
	import random
	import matplotlib.pyplot as plt

	# test input values
	FeedbackOption = kwargs.get("FeedbackOption",'ff_only')
	ControlStrategy = kwargs.get("ControlStrategy", 'synergist')
	test_input_values(muscle_1_parameters,delay_1_parameters,gain_1_parameters,FeedbackOption=FeedbackOption, ControlStrategy=ControlStrategy)
	test_input_values(muscle_2_parameters,delay_2_parameters,gain_2_parameters,FeedbackOption=FeedbackOption, ControlStrategy=ControlStrategy)

	# muscle architectural parameters
	PennationAngle_1 = muscle_1_parameters['Pennation Angle']
	InitialMusculoTendonLength_1 = muscle_1_parameters['Initial Muscle Length']*np.cos(PennationAngle_1)+muscle_1_parameters['Initial Tendon Length'] 
	PennationAngle_2 = muscle_2_parameters['Pennation Angle']	
	InitialMusculoTendonLength_2 = muscle_2_parameters['Initial Muscle Length']*np.cos(PennationAngle_2)+muscle_2_parameters['Initial Tendon Length'] 
	
	# viscosity of muscle (ref: Elias et al. 2014)
	MuscleViscosity = 0.005 #0.001	

	SamplingFrequency = 10000
	FeedbackSamplingFrequency = 1000
	SamplingRatio = SamplingFrequency/FeedbackSamplingFrequency
	SamplingPeriod = 1/SamplingFrequency
	Time = np.arange(0,len(TargetTrajectory_1)*SamplingPeriod,SamplingPeriod) 

	# filter parameters for Input['Noise']
	BButtersCoefficients,AButtersCoefficients = signal.butter(4,100/(SamplingFrequency/2),'low')

	# discrete transfer function for GTO output
	Num,Den = [1.7,2.58,0.4],[1,2.2,0.4]
	ContinuousTransferFunction = control.tf(Num,Den)
	DiscreteTransferFunction = control.matlab.c2d(ContinuousTransferFunction,1/FeedbackSamplingFrequency)
	Num,Den = control.matlab.tfdata(DiscreteTransferFunction)
	Num,Den = Num[0][0],Den[0][0]

	# Convert delays in ms to samples
	EfferentDelay_1 = delay_1_parameters['Efferent Delay']
	EfferentDelayTimeStep_1 = int(EfferentDelay_1*SamplingRatio)
	EfferentDelay_2 = delay_2_parameters['Efferent Delay']	
	EfferentDelayTimeStep_2 = int(EfferentDelay_2*SamplingRatio)
	
	# define all dictionaries needed for muscle(s)
	Bag1_1,Bag2_1,Chain_1,SlowTwitch_1,FastTwitch_1,CE_1,SEE_1,PEE_1,Muscle_1,Input_1 = return_initial_values(muscle_1_parameters,gain_1_parameters,TargetTrajectory_1,CorticalInput_1)
	Bag1_2,Bag2_2,Chain_2,SlowTwitch_2,FastTwitch_2,CE_2,SEE_2,PEE_2,Muscle_2,Input_2 = return_initial_values(muscle_2_parameters,gain_2_parameters,TargetTrajectory_2,CorticalInput_2)
	
	# Account for synergist strategy by changing target force trajectories to be X %MVC OF THE SUM of CE['Maximum Force']
	# NOTE: We only take Trajectory_1 in this case to make for uniform strategy shapes for synergists.
	if ControlStrategy == 'synergist':
		Input_1['Target Force Trajectory'],Input_2['Target Force Trajectory'] = [TargetTrajectory_1*(CE_1['Maximum Force'] + CE_2['Maximum Force'])]*2
	
	random.seed()
	StartTime = time.time()
	for i in range(len(Time)): 
		update_spindle_afferent_inputs_at_step_i_multiple_muscles(i,Input_1,CE_1,SEE_1,Bag1_1,Bag2_1,Chain_1,Num,Den,\
																delay_1_parameters,gain_1_parameters,SamplingRatio,SamplingPeriod,\
																FeedbackOption,ControlStrategy=ControlStrategy,SynergistParameters=[SEE_2,CE_2])
		update_spindle_afferent_inputs_at_step_i_multiple_muscles(i,Input_2,CE_2,SEE_2,Bag1_2,Bag2_2,Chain_2,Num,Den,\
																delay_2_parameters,gain_2_parameters,SamplingRatio,SamplingPeriod,\
																FeedbackOption,ControlStrategy=ControlStrategy,SynergistParameters=[SEE_1,CE_1])	
		
		# add the interneuron input to total input
		# add_interneuron_inputs_to_total_at_step_i(i,Input_1,Input_2,\
		# 	gain_1_parameters,delay_1_parameters,gain_2_parameters,delay_2_parameters,\
		# 	'inhibitory','inhibitory', SamplingRatio)
		
		#add noise (and cortical input) to input
		# random.seed(1)
		add_noise_to_input(i,Input_1,AButtersCoefficients,BButtersCoefficients)
		add_noise_to_input(i,Input_2,AButtersCoefficients,BButtersCoefficients)

		# add delay along efferent pathway
		if i > EfferentDelayTimeStep_1: #NEED TO MAKE THIS >= after next run!
			Input_1['Long'].append(Input_1['Total'][i-EfferentDelayTimeStep_1])
		else:
			Input_1['Long'].append(0)
		if i > EfferentDelayTimeStep_2: #NEED TO MAKE THIS >= after next run!				
			Input_2['Long'].append(Input_2['Total'][i-EfferentDelayTimeStep_2])
		else:
			Input_2['Long'].append(0)

		# apply activation filter
		Muscle_1['Effective Activation'].append(apply_activation_filter(Input_1,Muscle_1,SamplingPeriod))
		Muscle_2['Effective Activation'].append(apply_activation_filter(Input_2,Muscle_2,SamplingPeriod))	
		
		# SlowTwitch_1 = activation_frequency_slow(CE_1,SlowTwitch_1,Muscle_1['Effective Activation'][-1],SamplingPeriod) # not used
		# SlowTwitch_2 = activation_frequency_slow(CE_2,SlowTwitch_2,Muscle_2['Effective Activation'][-1],SamplingPeriod) # not used
		
		# update all kinematics and kinetics from Effective Muscle Activation at timestep i
		update_kinematics_and_kinetics(CE_1,PEE_1,SEE_1,Muscle_1,MuscleViscosity,PennationAngle_1,InitialMusculoTendonLength_1,SamplingPeriod)
		update_kinematics_and_kinetics(CE_2,PEE_2,SEE_2,Muscle_2,MuscleViscosity,PennationAngle_2,InitialMusculoTendonLength_2,SamplingPeriod)	
		
		statusbar(i,len(Time),StartTime=StartTime,Title = '2 Muscles')

	plt.figure()
	plt.plot(Time,SEE_1['Tendon Force'][1:])
	plt.plot(Time,Input_1['Target Force Trajectory'],'r')
	plt.legend(['Output Force','Target Force'])
	plt.xlabel('Time (sec)')
	plt.ylabel('Force (N)')

	plt.figure()
	plt.plot(Time,SEE_2['Tendon Force'][1:])
	plt.plot(Time,Input_2['Target Force Trajectory'],'r')
	plt.legend(['Output Force','Target Force'])
	plt.xlabel('Time (sec)')	
	plt.ylabel('Force (N)')

	if ControlStrategy == 'synergist':
		plt.figure()
		plt.plot(Time,np.array(SEE_2['Tendon Force'][1:])+np.array(SEE_1['Tendon Force'][1:]))
		plt.plot(Time,Input_2['Target Force Trajectory'],'r')
		plt.legend(['Output Force','Target Force'])
		plt.xlabel('Time (sec)')	
		plt.ylabel('Force (N)')
		plt.title('Combined Force')

	# save data as output in structure format
	output_1 = {	'Force' : Muscle_1['Force'], 									'ForceTendon' : SEE_1['Tendon Force'][1:], \
					'FL' : CE_1['FL'], 												'FV' : CE_1['FV'], \
					'PE Force 1' : PEE_1['Passive Force 1'], 						'PE Force 2' : PEE_1['Passive Force 2'], \
					'Target' : Input_1['Target Force Trajectory'], 					'ContractileElementLength' : CE_1['Length'][:-1], \
					'ContractileElementVelocity' : CE_1['Velocity'][1:-1],			'ContractileElementAcceleration' : CE_1['Acceleration'][1:-1], \
					'SeriesElasticElementLength' : SEE_1['Tendon Length'][:-1],		'Activation Frequency' : [0.0]*149999, \
					'Input' : Input_1['Total'], 									'Noise' : Input_1['Noise'], \
					'FilteredNoise' : Input_1['FilteredNoise'], 					'U' : Muscle_1['Effective Activation'][1:-1], \
					'Ia Input' : Input_1['Ia'],										'II Input' : Input_1['II'], \
					'Ib Input' : Input_1['Ib'],										'ALL'  :  Input_1 	} #, 'CorticalInput' : CorticalInput_1 }
	output_2 = {	'Force' : Muscle_2['Force'], 									'ForceTendon' : SEE_2['Tendon Force'][1:], \
					'FL' : CE_2['FL'], 												'FV' : CE_2['FV'], \
					'PE Force 1' : PEE_2['Passive Force 1'], 						'PE Force 2' : PEE_2['Passive Force 2'], \
					'Target' : Input_2['Target Force Trajectory'], 					'ContractileElementLength' : CE_2['Length'][:-1], \
					'ContractileElementVelocity' : CE_2['Velocity'][1:-1],			'ContractileElementAcceleration' : CE_2['Acceleration'][1:-1], \
					'SeriesElasticElementLength' : SEE_2['Tendon Length'][:-1],		'Activation Frequency' : [0.0]*149999, \
					'Input' : Input_2['Total'], 									'Noise' : Input_2['Noise'], \
					'FilteredNoise' : Input_2['FilteredNoise'], 					'U' : Muscle_2['Effective Activation'][1:-1], \
					'Ia Input' : Input_2['Ia'],										'II Input' : Input_2['II'], \
					'Ib Input' : Input_2['Ib'],										'ALL'  :  Input_2 	} #, 'CorticalInput' : CorticalInput_2 }
	"""
	NOTES:

	'ContractileElementVelocity','ContractileElementAcceleration', and 'ForceTendon' are from 1: because originally the lists were not initialized even though Vce and Ace are initialized to zero
	'EffectiveMuscleActivation' has the same problem

	"""
	return(output_1,output_2)

muscle_1_parameters = {	"Pennation Angle":5*np.pi/180, 		"Muscle Mass":0.075,\
						"Optimal Length":10.1, 				"Tendon Length":23.5,\
					 	"Initial Muscle Length":10.1, 	"Initial Tendon Length":23.5}
muscle_2_parameters = {	"Pennation Angle":5*np.pi/180, 		"Muscle Mass":0.075,\
						"Optimal Length":10.1, 				"Tendon Length":23.5,\
					 	"Initial Muscle Length":10.1, 	"Initial Tendon Length":23.5}

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
gain_1_parameters = {	"Gamma Dynamic Gain" : 50, \
					"Gamma Static Gain" : 50, \
					"Ia Gain" : 2000, \
					"II Gain" : 3000, \
					"Ib Gain" : 5000, \
					"Ia Reciprocal Gain" : 2000}
gain_2_parameters = {	"Gamma Dynamic Gain" : 50, \
					"Gamma Static Gain" : 50, \
					"Ia Gain" : 2000, \
					"II Gain" : 3000, \
					"Ib Gain" : 5000, \
					"Ia Reciprocal Gain" : 2000}

# Define force trajectory that model needs to track
# You can create any target trajectory as you want, but it has to be
# sampled at 10000 Hz to be consistent with the sampling frequency of
# muscle model
SamplingFrequency = 10000
Time = np.arange(0,15,1/SamplingFrequency) #0:SamplingPeriod:15
TrajectoryType = 'constant' #'triangle','sinwave','trapezoid' (Still need to do Trap)
TargetTrajectory_1 = generate_target_force_trajectory(TrajectoryType,Time,0.3,0,0) # in unit of #MVC
TargetTrajectory_2 = generate_target_force_trajectory(TrajectoryType,Time,0.3,0,0) # in unit of #MVC

# Define additional input to muscle (e.g., 20 Hz oscillatory input)
# CorticalInput_1Amplitude = (1/10)*TargetTrajectory_1 # in unit of #MVC
# CorticalInput_1Frequency = 20 # in unit of Hz
# CorticalInput_1 = CorticalInput_1Amplitude*np.sin(2*np.pi*CorticalInput_1Frequency*Time)
CorticalInput_1 = np.zeros(len(TargetTrajectory_1))
CorticalInput_2 = np.zeros(len(TargetTrajectory_1))

FeedbackOption = 'fb_control'
if FeedbackOption == 'ff_only':
    FeedbackOptionString = 'FF'
    TargetTrajectory_1 = 1.028*TargetTrajectory_1
    TargetTrajectory_2 = 1.028*TargetTrajectory_2
elif FeedbackOption == 'servo_control':
    FeedbackOptionString = 'SC'
    TargetTrajectory_1 = 0.9285*TargetTrajectory_1
    TargetTrajectory_2 = 0.9285*TargetTrajectory_2
elif FeedbackOption == 'fb_control':
    FeedbackOptionString = 'FB'
elif FeedbackOption == 'cortical_fb_only':
    FeedbackOptionString = 'CC'

ControlStrategy = 'synergist'

# Output from the muscle model
NumberOfTrials = 1
Output_1,Output_2 = [],[]
PXX_1,PXX_2, PXX_total = [],[],[]
Range_1,Range_2 = [],[]
for i in range(NumberOfTrials): #trialN = 1:10
	#rng('shufle')
	output_1,output_2 = afferented_2_muscles_model(muscle_1_parameters,muscle_2_parameters,\
													delay_1_parameters,delay_2_parameters,\
													gain_1_parameters,gain_2_parameters,\
													TargetTrajectory_1,TargetTrajectory_2,\
													CorticalInput_1,CorticalInput_2,\
													FeedbackOption = FeedbackOption,\
													ControlStrategy = ControlStrategy)    
	# compare_output(output)
	Output_1.append(output_1)
	Output_2.append(output_1)
	f_1,pxx_1 = welch(output_1['ForceTendon'][-int(10*SamplingFrequency)-1:]-np.average(output_1['ForceTendon'][int(-10*SamplingFrequency-1):]),window = gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5)),\
		noverlap = SamplingFrequency,nperseg = len(gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5))), fs = SamplingFrequency) 
	f_2,pxx_2 = welch(output_2['ForceTendon'][-int(10*SamplingFrequency)-1:]-np.average(output_2['ForceTendon'][int(-10*SamplingFrequency-1):]),window = gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5)),\
		noverlap = SamplingFrequency,nperseg = len(gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5))), fs = SamplingFrequency) 
	output_total = ((np.array(output_1['ForceTendon'][-int(10*SamplingFrequency)-1:])+np.array(output_2['ForceTendon'][-int(10*SamplingFrequency)-1:]))\
						- np.average(np.array(output_1['ForceTendon'][-int(10*SamplingFrequency)-1:])+np.array(output_2['ForceTendon'][-int(10*SamplingFrequency)-1:])))
	f_total,pxx_total = welch(output_total,window = gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5)),\
		noverlap = SamplingFrequency,nperseg = len(gaussian(5*SamplingFrequency,(5*SamplingFrequency-1)/(2*2.5))), fs = SamplingFrequency) 
	pxx_1 = smooth(pxx_1,5)
	pxx_2 = smooth(pxx_2,5)
	pxx_total = smooth(pxx_total,5)
	PXX_1.append(np.array(pxx_1,ndmin=2))
	PXX_2.append(np.array(pxx_2,ndmin=2))
	PXX_total.append(np.array(pxx_total,ndmin=2))
	Range_1.append(max(output_1['Ia Input'][-10*SamplingFrequency:-1])-min(output_1['Ia Input'][-10*SamplingFrequency-1:]))
	Range_2.append(max(output_2['Ia Input'][-10*SamplingFrequency:-1])-min(output_2['Ia Input'][-10*SamplingFrequency-1:]))
	print("\n")

PXX_1 = np.concatenate(PXX_1,axis=0)
plt.figure()
plt.plot(f_1[:131],np.average(PXX_1[:,:131],axis = 0))

PXX_2 = np.concatenate(PXX_2,axis=0)
plt.figure()
plt.plot(f_2[:131],np.average(PXX_2[:,:131],axis = 0))

PXX_total = np.concatenate(PXX_total,axis=0)
plt.figure()
plt.plot(f_total[:131],np.average(PXX_total[:,:131],axis = 0))
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