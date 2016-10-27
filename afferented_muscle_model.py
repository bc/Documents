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
	InitialMusculoTendonLength = muscle_parameters['Initial Muscle Length']*np.cos(PennationAngle)+muscle_parameters['Initial Tendon Length'] 

	# viscosity of muscle (ref: Elias et al. 2014)
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

	# filter parameters for Input['Noise']
	BButtersCoefficients,AButtersCoefficients = signal.butter(4,100/(SamplingFrequency/2),'low')

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

	Bag1,Bag2,Chain,SlowTwitch,FastTwitch,CE,SEE,PEE,Muscle,Input = return_initial_values(muscle_parameters,gain_parameters)

	# Convert force trajectory to unit of newton and add (with FF Input) to Input Dict
	Input['Target Force Trajectory'] = TargetTrajectory*CE['Maximum Force']
	Input['Feedforward'] = Input['Target Force Trajectory']/CE['Maximum Force']

	StartTime = time.time()
	for i in range(len(Time)): #= 1:length(Time)
		if FeedbackOption == 'ff_only': # Feedforward input onl
			Input['Total'].append(Input['Feedforward'][i])
			Input['Ia'].append(0)
			Input['II'].append(0)
			Input['Ib'].append(0)
		elif FeedbackOption == 'servo_control': # Servo control (feedforward + spindle and GTO)
			PrimaryOutput,SecondaryOutput = spindle_model(CE,Bag1,Bag2,Chain,SamplingPeriod)
			Input['Ia'].append(PrimaryOutput)
			Input['II'].append(SecondaryOutput)
			if i == Count:
				# Get spindle primary and secondary afferent activity

				# Get Ib activity
				Input['IbTemp1'].append(GTOConstant1*np.log(SEE['Force'][-1]/GTOConstant2+1))
				if i == 0:
					Input['IbTemp2'].append(Input['IbTemp1'][-1])
				elif i == 1*SamplingRatio:
					Input['IbTemp2'].append(Input['IbTemp1'][-1])
				elif i == 2*SamplingRatio:
					Input['IbTemp2'].append((Num[0]*Input['IbTemp1'][-1])/Den[0])
				elif i == 3*SamplingRatio:
					Input['IbTemp2'].append((Num[1]*Input['IbTemp1'][-2] + Num[0]*Input['IbTemp1'][-1] \
						- Den[1]*Input['IbTemp2'][-2])/Den[0])
				elif i >= 3*SamplingRatio:
					Input['IbTemp2'].append((Num[2]*Input['IbTemp1'][-3] + Num[1]*Input['IbTemp1'][-2] + Num[0]*Input['IbTemp1'][-1] \
						- Den[2]*Input['IbTemp2'][-3] - Den[1]*Input['IbTemp2'][-2])/Den[0])

				Input['Ib'].append(Input['IbTemp2'][-1]*(Input['IbTemp2'][-1]>0))
			
				if i in range(IaAfferentDelayTimeStep-1,IbAfferentDelayTimeStep):
					Input['Total'].append(Input['Ia'][i-IaAfferentDelayTimeStep+1]/IaSpindleGain \
								+Input['Feedforward'][i]) #input to the muscle
				elif i in range(IbAfferentDelayTimeStep, IIAfferentDelayTimeStep):
					Input['Total'].append(Input['Ia'][i-IaAfferentDelayTimeStep+1]/IaSpindleGain \
								-Input['Ib'][i-IbAfferentDelayTimeStep]/IbGTOGain \
								+Input['Feedforward'][i]) #input to the muscle
				elif i in range(IIAfferentDelayTimeStep, CorticalDelayTimeStep):
					Input['Total'].append(Input['Ia'][i-IaAfferentDelayTimeStep+1]/IaSpindleGain \
								+Input['II'][i-IIAfferentDelayTimeStep]/IISpindleGain \
								-Input['Ib'][i-IbAfferentDelayTimeStep]/IbGTOGain \
								+Input['Feedforward'][i]) #input to the muscle
				else:	
					Input['Total'].append(Input['Feedforward'][i])
				
				Count = SamplingRatio+Count
			else:
				Input['Ib'].append(Input['Ib'][-1])
				Input['Total'].append(Input[-1])
		elif FeedbackOption == 'fb_control': # Feedback control (proprioceptive systems + supraspinal loop)
			PrimaryOutput,SecondaryOutput = spindle_model(CE,Bag1,Bag2,Chain,SamplingPeriod)
			Input['Ia'].append(PrimaryOutput)
			Input['II'].append(SecondaryOutput)
			if i == Count:
				# Get spindle primary and secondary afferent activity
				
				# Get Ib activity
				Input['IbTemp1'].append(GTOConstant1*np.log((SEE['Force'][-1]/GTOConstant2+1)))
				if i == 0:
					Input['IbTemp2'].append(Input['IbTemp1'][-1])
				elif i == 1*SamplingRatio:
					Input['IbTemp2'].append(Input['IbTemp1'][-1])
				elif i == 2*SamplingRatio:
					Input['IbTemp2'].append((Num[0]*Input['IbTemp1'][-1])/Den[0])
				elif i == 3*SamplingRatio:
					Input['IbTemp2'].append((Num[1]*Input['IbTemp1'][-2] + Num[0]*Input['IbTemp1'][-1]- Den[1]*Input['IbTemp2'][-1])/Den[0])
				elif i > 3*SamplingRatio:
					Input['IbTemp2'].append((Num[2]*Input['IbTemp1'][-3] + Num[1]*Input['IbTemp1'][-2] + Num[0]*Input['IbTemp1'][-1] \
														- Den[2]*Input['IbTemp2'][-2] - Den[1]*Input['IbTemp2'][-1])/Den[0])

				Input['Ib'].append(Input['IbTemp2'][-1]*(Input['IbTemp2'][-1]>0))

				if i in range(IaAfferentDelayTimeStep-1,IbAfferentDelayTimeStep):
					Input['Total'].append(Input['Ia'][i-IaAfferentDelayTimeStep+1]/IaSpindleGain\
								+ Input['Feedforward'][i]	) #input to the muscle
				elif i in range(IbAfferentDelayTimeStep, IIAfferentDelayTimeStep):
					Input['Total'].append(Input['Ia'][i-IaAfferentDelayTimeStep+1]/IaSpindleGain \
								- Input['Ib'][i-IbAfferentDelayTimeStep]/IbGTOGain \
								+ Input['Feedforward'][i] ) #input to the muscle
				elif i in range(IIAfferentDelayTimeStep, CorticalDelayTimeStep):
					Input['Total'].append(Input['Ia'][i-IaAfferentDelayTimeStep+1]/IaSpindleGain \
								+ Input['II'][i-IIAfferentDelayTimeStep]/IISpindleGain \
								- Input['Ib'][i-IbAfferentDelayTimeStep]/IbGTOGain \
								+ Input['Feedforward'][i]	) #input to the muscle
				elif i >= CorticalDelayTimeStep:
					FeedbackInput = TransCorticalLoopConstant*(Input['Target Force Trajectory'][i]-SEE['Force'][i-CorticalDelayTimeStep+1])/CE['Maximum Force'] + FeedbackInput  # feedback input through cortical pathway
					Input['Total'].append(Input['Ia'][i-IaAfferentDelayTimeStep+1]/IaSpindleGain \
								+Input['II'][i-IIAfferentDelayTimeStep]/IISpindleGain \
								-Input['Ib'][i-IbAfferentDelayTimeStep]/IbGTOGain \
								+FeedbackInput	) #input to the muscle
				else:	
					Input['Total'].append(Input['Feedforward'][i])
				
				Count = SamplingRatio+Count
			else:
				Input['Ib'].append(Input['Ib'][-1])
				Input['Total'].append(Input['Total'][-1])
		elif FeedbackOption == 'cortical_fb_only':
			if i > CorticalDelayTimeStep-1:
				FeedbackInput = TransCorticalLoopConstant*(Input['Target Force Trajectory'][i]-SEE['Force'][i-CorticalDelayTimeStep+1])/CE['Maximum Force'] + FeedbackInput  # feedback input through cortical pathway
				Input['Total'].append(FeedbackInput) #input to the muscle
			else:
				Input['Total'].append(Input['Feedforward'][i])

			Input['Ia'].append(0)
			Input['II'].append(0)
			Input['Ib'].append(0)

		# Feedforward only
		# Input[i] = input[i] #Input['Target Force Trajectory'][i]/CE['Maximum Force']
		## Input['Noise'] + additional input
		Input['Total'][-1] = bound(Input['Total'][-1],0,1)

		random.seed(1)
		if i > 4:
			Input['Noise'].append(2*(random.random()-0.5)*(np.sqrt(0.01*Input['Total'][i])*np.sqrt(3)))
			Input['FilteredNoise'].append((BButtersCoefficients[4]*Input['Noise'][i-4] + BButtersCoefficients[3]*Input['Noise'][i-3] + BButtersCoefficients[2]*Input['Noise'][i-2] + BButtersCoefficients[1]*Input['Noise'][i-1] + BButtersCoefficients[0]*Input['Noise'][i] \
									- AButtersCoefficients[4]*Input['FilteredNoise'][i-4] - AButtersCoefficients[3]*Input['FilteredNoise'][i-3] - AButtersCoefficients[2]*Input['FilteredNoise'][i-2] - AButtersCoefficients[1]*Input['FilteredNoise'][i-1])/AButtersCoefficients[0])
		else:
			Input['Noise'].append(0)
			Input['FilteredNoise'].append(0)

		Input['Total'][-1] = Input['Total'][-1] + Input['FilteredNoise'][-1] + CorticalInput[-1]

		Input['Total'][-1] = bound(Input['Total'][-1],0,1)

		# add delay along efferent pathway
		if i > EfferentDelayTimeStep: #NEED TO MAKE THIS >= after next run!
			Input['Long'].append(Input['Total'][i-EfferentDelayTimeStep])
		else:
			Input['Long'].append(0)

		# add activation filter
		TU = 0.12*(Input['Long'][-1] < Muscle['Effective Activation'][-1]) + 0.03 # When true TU = 0.15, else TU = 0.03

		EffectiveMuscleActivationFirstDeriv = (Input['Long'][-1] - Muscle['Effective Activation'][-1])/TU
		Muscle['Effective Activation'].append(EffectiveMuscleActivationFirstDeriv*(SamplingPeriod) + Muscle['Effective Activation'][-1]) # effective neural drive

		# SlowTwitch = activation_frequency_slow(CE,SlowTwitch,Muscle['Effective Activation'][-1],SamplingPeriod) # not used

		# force-velocity relationship
		CE['FV'].append((CE['Velocity'][-1] <= 0)*concentric_force_velocity(CE) \
							+ (CE['Velocity'][-1] > 0)*eccentric_force_velocity(CE) )
		# force-length relationship
		CE['FL'].append(force_length(CE))
		CE['Force'] = CE['FL'][-1]*CE['FV'][-1]
		# viscous property
		ForceViscosity = MuscleViscosity * CE['Velocity'][-1]
		# passive element 1
		PEE['Passive Force 1'].append(parallel_elastic_element_force_1(CE))
		# passive element 2
		ForcePassive2 = parallel_elastic_element_force_2(CE)
		PEE['Passive Force 2'].append((ForcePassive2 <= 0)*ForcePassive2)

		# total force from contractile element
		ForceTotal = (Muscle['Effective Activation'][-1]*(CE['Force'] + PEE['Passive Force 2'][-1]) + PEE['Passive Force 1'][-1] + ForceViscosity)*CE['Maximum Force']
		Muscle['Force'].append(ForceTotal*(ForceTotal>=0.0))

		# force from series elastic element
		SEE['Force'].append(normalized_series_elastic_element_force(SEE) * CE['Maximum Force'])

		# calculate muscle excursion acceleration based on the difference
		# between muscle force and tendon force
		if i < len(Time)-1:
			Muscle['Acceleration'].append((SEE['Force'][-1]*np.cos(PennationAngle) - Muscle['Force'][-1]*(np.cos(PennationAngle))**2)/(Muscle['Mass']) \
				+ (Muscle['Velocity'][-1])**2*np.tan(PennationAngle)**2/(Muscle['Length'][-1]))
			# integrate acceleration to get velocity
			Muscle['Velocity'].append((Muscle['Acceleration'][-1]+ \
				Muscle['Acceleration'][-2])/2*(SamplingPeriod)+Muscle['Velocity'][-1])
			# integrate velocity to get length
			Muscle['Length'].append((Muscle['Velocity'][-1]+ \
				Muscle['Velocity'][-2])/2*(SamplingPeriod)+Muscle['Length'][-1])

			# normalize each variable to optimal muscle length or tendon legnth
			CE['Acceleration'].append(Muscle['Acceleration'][-1]/(Muscle['Optimal Length']/100))
			CE['Velocity'].append(Muscle['Velocity'][-1]/(Muscle['Optimal Length']/100))
			CE['Length'].append(Muscle['Length'][-1]/(Muscle['Optimal Length']/100))
			SEE['Length'].append((InitialMusculoTendonLength - CE['Length'][-1]*Muscle['Optimal Length']*np.cos(PennationAngle))/SEE['Optimal Length'])

		statusbar(i,len(Time),StartTime=StartTime,Title = 'afferented_muscle_model')

	# plt.figure()
	# plt.plot(Time,SEE['Force'][1:])
	# plt.plot(Time,Input['Target Force Trajectory'],'r')
	# plt.legend(['Output Force','Target Force'])
	# plt.xlabel('Time (sec)')
	# plt.ylabel('Force (N)')

	# save data as output in structure format
	output = {	'Force' : Muscle['Force'], 											'ForceTendon' : SEE['Force'][1:], \
				'FL' : CE['FL'], 												'FV' : CE['FV'], \
				'PE Force 1' : PEE['Passive Force 1'], 							'PE Force 2' : PEE['Passive Force 2'], \
				'Target' : Input['Target Force Trajectory'], 									'ContractileElementLength' : CE['Length'], \
				'ContractileElementVelocity' : CE['Velocity'][1:],				 	'ContractileElementAcceleration' : CE['Acceleration'][1:], \
				'SeriesElasticElementLength' : SEE['Length'],					 	'Activation Frequency' : [0.0]*149999, \
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