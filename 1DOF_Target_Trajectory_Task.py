# 1DOF Target Trajectory Task
# August 2016 - USC
# Daniel A Hagen

import numpy as np 
from scipy.optimize import linprog
import matplotlib.pyplot as plt 
import random
import sympy as sp 

random.seed()
#MomentArm1, MomentArm2, MomentArm3, MomentArm4 = [random.uniform(0.0,2.0) for i in range(4)]
MomentArm1, MomentArm2, MomentArm3, MomentArm4 = [1, 1.0001, 2, 1]
# Define Moment Arm Matrix (R) as 1X4 Matrix with two positive and two negative
# Moment Arm values to reflect a single joint with two muscle on either side
# of the joint.
R = np.matrix([MomentArm1, MomentArm2, -MomentArm3, -MomentArm4])
force_length_curve = lambda x: 1-(x/0.5)**2 if np.abs(x)<=0.5 else 0
def force_velocity_curve(NormalizedMuscleVelocity,b):
		"""
		NormalizedMuscleVelocity is a scalar and b is a coefficient
		that can be changed to adjust the shape of the curve.
		"""
		if NormalizedMuscleVelocity < -5:
			result = 0
		elif NormalizedMuscleVelocity < 0:
			result = (1 + NormalizedMuscleVelocity/5)/(b - NormalizedMuscleVelocity)
		else:
			result = (1.8 - 0.8*np.exp(-(b+5)*NormalizedMuscleVelocity/4))
		return(result)
def maximum_muscle_force(OptimalLengths, Angle, AngularVelocity, R, OptimalForces):
	"""
	OptimalLengths must be a 1X4 Matrix with optimal muscle lengths for 4 muscles.
	Angle should be the current angle. AngularVelocity should be the current 
	angular velocity. R is the 1X4 moment arm matrix. And OptimalForces should be
	a 4X1 matrix of forces produced by each muscle when at optimal lengths and 
	velocities.
	"""
	# Force-Length Considerations
	CurrentMuscleLengths = OptimalLengths.T - R.T*Angle
	NormalizedMuscleLengths = np.matrix([(CurrentMuscleLengths[i,0]/OptimalLengths.T[i,0])-1 for i in range(4)])
	# We subtract one from every Normalized Muscle Length to find the percentage 
	# above and below the optimal length.
	
	MaximumMuscleForce_FL = np.identity(4)*[force_length_curve(NormalizedMuscleLengths[0,i]) for i in range(4)]

	# Force-Velocity Considerations
	CurrentMuscleVelocity = -R.T*AngularVelocity
	NormalizedMuscleVelocity = [CurrentMuscleVelocity[i,0]/OptimalLengths.T[i,0] for i in range(4)]
	
	MaximumMuscleForce_FV = np.identity(4)* \
								[force_velocity_curve(NormalizedMuscleVelocity[i],1) \
									for i in range(4)]
	MaximumMuscleForce = MaximumMuscleForce_FL*MaximumMuscleForce_FV*[elem for elem in OptimalForces]
	return(MaximumMuscleForce,MaximumMuscleForce_FL,MaximumMuscleForce_FV)

def create_trajectory_and_derivatives(Amplitude, FrequencyOfOscillations, Time):
	time = sp.symbols('time',real = True)
	angle = Amplitude*sp.sin(2*np.pi*FrequencyOfOscillations*time)
	angular_velocity = angle.diff(time)
	angular_acceleration = angle.diff(time,2)
	Angle = np.array([angle.subs([(time,elem)]) for elem in Time])
	AngularVelocity = np.array([angular_velocity.subs([(time,elem)]) for elem in Time])
	AngularAcceleration = np.array([angular_acceleration.subs([(time,elem)]) for elem in Time])
	return(Angle,AngularVelocity,AngularAcceleration)

Time = np.arange(0,10,0.01)
Amplitude = np.pi/4
FrequencyOfOscillations = 1
Angle, AngularVelocity, AngularAcceleration = create_trajectory_and_derivatives(Amplitude,FrequencyOfOscillations,Time)
Mass = 10
LinkLength = 1
OptimalLengths = np.matrix([10,10,10,10])
OptimalForces = np.matrix([200,200,200,200])
Inertia = (1/3)*Mass*LinkLength**2
TorqueDueToGravity = [-0.5*Mass*9.8*LinkLength*np.cos(float(theta)) for theta in Angle]
SumOfMuscleTorques = [Inertia*float(AngularAcceleration[i]) - TorqueDueToGravity[i] for i in range(len(Time))]

def construct_A_matrix(OptimalLengths, Angle, AngularVelocity, R, OptimalForces):
	MaximumForceMatrix,_,_ = maximum_muscle_force(OptimalLengths, Angle, AngularVelocity, R, OptimalForces)
	A = np.concatenate((R*MaximumForceMatrix,-R*MaximumForceMatrix),axis=0)
	return(np.array(A))
def construct_b_vector(SumOfMuscleTorques, Epsilon):
	b = np.array([SumOfMuscleTorques+Epsilon,-SumOfMuscleTorques+Epsilon])
	return(b)
Cost = np.array([1,1,1,1])
Epsilon = 0.01

# Find and print optimal solution

a = [linprog(Cost,A_ub = construct_A_matrix(OptimalLengths, float(Angle[i]), float(AngularVelocity[i]), R, OptimalForces),\
					 b_ub = construct_b_vector(float(SumOfMuscleTorques[i]),Epsilon),\
						bounds = ((0,1),(0,1),(0,1),(0,1))) for i in range(len(Time))]
a_star = np.concatenate([[a[i]['x'] for i in range(len(Time))]])

plt.figure()
ax1 = plt.gca()
plt.plot(Time,a_star)
plt.xlabel('Time (s)')
plt.ylabel('Muscle Activation')
plt.legend(['Muscle 1','Muscle 2','Muscle 3','Muscle 4'])
plt.title("R = [{}, {}, {}, {}]"\
		.format(R[0,0],R[0,1],R[0,2],R[0,3]))

plt.figure()
ax2 = plt.gca()
LengthDomain = np.arange(-0.6,0.6,0.001)
F_L_Curve = [force_length_curve(x) for x in LengthDomain]
plt.plot(LengthDomain,F_L_Curve)
plt.title('Force - Normalized Length Curve')
plt.xlabel('Percent Change in Normalized Muscle Length')
plt.ylabel('Muscle Force Percentage')

plt.figure()
ax3 = plt.gca()
VelocityDomain = np.arange(-6,6,0.01)
F_V_Curve = [force_velocity_curve(x,1) for x in VelocityDomain]
plt.plot(VelocityDomain, F_V_Curve)
plt.title('Force - Normalized Velocity Curve')
plt.xlabel('Muscle Lengths per Second')
plt.ylabel('Muscle Force Percentage')

plt.show((ax1,ax2,ax3))



