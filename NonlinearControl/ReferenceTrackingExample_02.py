import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import sympy as sy
from sympy.utilities import lambdify
import time

def return_muscle_settings():
	"""
	Notes:
	Coefficients from observation, Ramsay, FVC, Holtzbaur, Pigeon, Kuechle, or Banks. Optimal Muscle Length given in mm.

	BIC EFE MA for Ramsay has R² = 0.985 whereas Pigeon has R² = 0.9918. Pigeon, however, only takes elbow angle into account, whereas Ramsay takes in variable PS angles. It appears that because Pigeon uses an average of fully pronated and fully supinated MAs, the BIC moment arm is similar but subject to variation as the level of PS is changed. Coefficients and equation number/type are listed below to test either implementation. (NOTE: BIC becomes slightly negative when x1 > 3.021. If trajectory has elbow angles exceding this value, enter a threshold of 3.021 into the model.)

	Additionally, the SFE MA for the BIC is held constant in Pigeon at 29.21 mm while it was estimated as 15 mm.

	src = 'Ramsay', eq = 2, Coefficients = [8.4533,36.6147,2.4777,-19.432,2.0571,0,13.6502,0,0,-5.6172,0,-2.0854,0,0,0,0], threshold = 3.021
	src = 'Pigeon', dof = 'elbow', eq = None, Coefficients = Pigeon_coeff_conversion([14.660,4.5322,1.8047,-2.9883,0,0]), threshold = 2.9326

	TRI EFE MA for Ramsay has R² = 0.997 whereas Pigeon has R² = 0.9904. Pigeon appears to really fail when the elbow angle is greater than 140°. For this reason, Ramsay should be used. However the approach of fixing the MA for values greater than 140° can be adopted for completeness. Coefficients and equation number/type are listed below to test either implementation.

	Additionally, the SFE MA for the TRI is held constant in Pigeon at -25.40 mm while it was estimated as -15 mm.

	src = 'Ramsay', eq = 1, Coefficients = [-24.5454,-8.8691,9.3509,-1.7518,0]
	src = 'Pigeon', dof = 'elbow', eq = None, Coefficients = Pigeon_coeff_conversion([-23.287,-3.0284,12.886,-19.092,13.277,-3.5171])

	"""
	import sympy as sp
	from sympy.utilities import lambdify
	import numpy as np
	from numpy import pi

	global θ_SFE,θ_EFE,θ_PS
	θ_SFE,θ_EFE,θ_PS = sp.symbols('θ_SFE'),sp.symbols('θ_EFE'),sp.symbols('θ_PS')

	# Coefficients from observation, Ramsay, Pigeon, FVC, Holtzbaur, or Banks.
	# Moment arms are in mm. Mass is in grams. threshold is in radians.

	def Pigeon_coeff_conversion(Coefficients):
		"""
		Takes in Coefficient values from Pigeon (1996) -- which take in angles in degrees -- and coverts them into the properly scaled coefficients for radians, additionally scaled by the magnitude listed in the paper.

		Note that the coefficients listed in Pigeon (1996) are given in decending order (i.e., c5,c4,c3,c₂,c₁,c0). However to maintain continuity with the equations given in Ramsay (2009), we list coefficients in order of increasing power (i.e., c0,c1,c2,c3,c4,c5).
		"""
		import numpy as np
		assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
		assert type(Coefficients)==list, 'Coefficients must be a 6 element list.'
		Rad_Conversion = np.multiply(Coefficients,\
				np.array([1,(180/np.pi),(180/np.pi)**2,(180/np.pi)**3,(180/np.pi)**4,(180/np.pi)**5],dtype = 'float64'))
		new_Coefficients =\
			np.multiply(Rad_Conversion,np.array([1,1e-1,1e-3,1e-5,1e-7,1e-9],dtype='float64'))
		return(new_Coefficients)

	BIC_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : [29.21,0,0,0,0,0],\
			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : [8.4533,36.6147,2.4777,-19.432,2.0571,0,13.6502,0,0,-5.6172,0,-2.0854,0,0,0,0],\
			'Source' : 'Ramsay', 'Equation Number' : 2, 'Threshold' : 3.021, \
			'dof' : 'Elbow'}, \
		'Mass' : 163.8,\
		'Actual No' : 320,\
		'Corrected No' : 292.6,\
		'Relative Abundance' : 1.1,\
		'Optimal Muscle Length' : 116,\
		'Group' : 'flexor'}
	TRI_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : [-25.40,0,0,0,0,0], \
			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : [-24.5454,-8.8691,9.3509,-1.7518,0],\
			'Source' : 'Ramsay', 'Equation Number' : 1, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : (94.2+138.4+92.5), \
		'Actual No' : (200+222+98),\
		'Corrected No' : (223.7+269.6+221.8),\
		'Relative Abundance' : (0.89+0.82+0.44)/3,\
		'Optimal Muscle Length' : 134,\
		'Group' : 'extensor'}

	AllMuscleSettings = {'BIC' : BIC_Settings, 'TRI' : TRI_Settings}
	return(AllMuscleSettings)
def return_MA_matrix_functions(AllMuscleSettings):
	import numpy as np
	import sympy as sp
	from sympy.utilities import lambdify
	def MA_function(Parameters):
		"""
		Note:

		Angles should be a number if Coefficients has a length of 5, or a list of length 2 when the Coefficients have lengths 16 or 18. Angles[0] will be the PRIMARY ANGLE for the DOF being considered while Angles[1] will be the secondary angle.

		Notes:

		threshold is only needed for Pigeon or Ramsay MA functions that are invalid outside of a given value. Must be either None (default) or the radian value of the threshold.

		dof only needed for Pigeon (Ramsay only handles EFE for this 2 DOF system). Must be either 'Shoulder' or 'Elbow'.

		eq is only needed for Ramsay (Pigeon has one quintic polynomial). eq must be either 1, 2, or 3, with list length requirements of 5, 16, or 18, respectively.
		"""
		import sympy as sp
		import numpy as np

		src = Parameters['Source']
		Coefficients = Parameters['MA Coefficients']
		eq = Parameters['Equation Number']
		dof = Parameters['dof']
		threshold = Parameters['Threshold']

		global θ_SFE,θ_EFE,θ_PS
		assert type(src) == str, "src must be a str."
		assert src.capitalize() in ['Ramsay','Pigeon','Est'], "src must be either Ramsay, Pigeon or Est (Estimate)."
		if dof != None:
			assert type(dof) == str, "dof must be a str."
			assert dof.capitalize() in ['Shoulder','Elbow'], "dof must be either Shoulder or Elbow."
		if src.capitalize() == 'Pigeon' :
			assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
			assert dof != None, "For Pigeon (1996), dof must be stated."
			eq = None
			if dof.capitalize() == 'Elbow' :
				θ = θ_EFE
			else:
				θ = θ_SFE
			MomentArm = (np.matrix(Coefficients,dtype='float64')\
							*np.matrix([1,θ,θ**2,θ**3,θ**4,θ**5]).T)[0,0]/1000
		elif src.capitalize() == 'Est' :
			MomentArm = np.array(Coefficients,dtype='float64')/1000
		else: #src.capitalize() == 'Ramsay'
			θ = θ_EFE
			assert type(Coefficients) == list, "Coefficients must be a list."
			assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
			assert eq in [1,2,3], "eq must be either 1, 2, or 3 when using Ramsay (2009)."
			if eq == 1:
				assert len(Coefficients) == 5, "For Eq. 1, Coefficients must be 5 elements long."
				MomentArm = (sp.Matrix(Coefficients,dtype='float64').T\
								*sp.Matrix([1,θ,θ**2,θ**3,θ**4]))[0,0]/1000
			elif eq == 2:
				assert len(Coefficients)==16, "For Eq. 2, Coefficients must be 16 elements long."
				MomentArm = (sp.Matrix(Coefficients,dtype='float64').T*\
								sp.Matrix([1, θ, θ_PS, θ*θ_PS, θ**2, \
											θ_PS**2, (θ**2)*θ_PS, θ*(θ_PS**2), \
											(θ**2)*(θ_PS**2), θ**3, θ_PS**3, \
											(θ**3)*θ_PS, θ*(θ_PS**3), \
											(θ**3)*(θ_PS**2), (θ**2)*(θ_PS**3), \
											(θ**3)*(θ_PS**3)]))[0, 0]/1000
			else: # eq == 3
				assert len(Coefficients)==18, "For Eq. 3, Coefficients must be 18 elements long."
				MomentArm = (sp.Matrix(Coefficients,dtype='float64').T*\
								sp.Matrix([1, θ, θ_PS, θ*θ_PS, θ**2, \
									θ_PS**2, (θ**2)*θ_PS, θ*(θ_PS**2), (θ**2)*(θ_PS**2), \
									θ**3, (θ**3)*θ_PS, (θ**3)*(θ_PS**2), \
									θ**4, (θ**4)*θ_PS, (θ**4)*(θ_PS**2),  \
									θ**5, (θ**5)*θ_PS, (θ**5)*(θ_PS**2)]))[0, 0]/1000
		if threshold == None:
			return(MomentArm)
		else:
			assert type(threshold) in [int,float], "threshold must be a number."
			MomentArm = sp.Piecewise((MomentArm,θ<threshold),(MomentArm.subs(θ,threshold),θ>=threshold))
			return(MomentArm)

	MuscleList = AllMuscleSettings.keys()

	RT_symbolic = sp.Matrix([MA_function(AllMuscleSettings[muscle]['Elbow']) for muscle in MuscleList])
	dRT_symbolic = sp.Matrix(sp.diff(RT_symbolic,θ_EFE))
	d2RT_symbolic = sp.Matrix(sp.diff(sp.diff(RT_symbolic,θ_EFE),θ_EFE))
	# RT_func = lambdify([θ_SFE,x1,θ_PS],RT_symbolic)
	# dRT_func = lambdify([θ_SFE,x1,θ_PS],dRT_symbolic)
	# import ipdb; ipdb.set_trace()
	# returns an (n,m) matrix when n is the number of muscles and m is the number of DOFS. We chose to return R.T because this is commonly utilized in muscle velocity calculations.
	return(RT_symbolic,dRT_symbolic,d2RT_symbolic)
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
	from scipy import interpolate
	import numpy as np
	StartTime = kwargs.get("StartTime",False)
	Title = kwargs.get("Title",'')
	global time_array
	global TimeLeft
	assert type(i)==int, "i must be an int"
	assert type(N)==int, "N must be an int"
	assert N>i, "N must be greater than i"
	assert N>0, "N must be a positive integer"
	assert i>=0, "i must not be negative (can be zero)"
	assert type(Title) == str, "Title should be a string"
	assert len(Title) <= 22, "Title should be less than 25 characters"
	if Title != '' : Title = ' '*(22-len(Title)) + Title + ' : '
	statusbar = Title +'[' + '\u25a0'*int((i+1)/(N/50)) + '\u25a1'*(50-int((i+1)/(N/50))) + '] '
	TimeBreak = abs
	if StartTime != False:
		if i==1:
			time_array = []
			TimeLeft = '--'
		elif i==int(0.02*N):
			time_array.append(time.time()-StartTime)
			TimeLeft = '{0:1.1f}'.format(time_array[-1]*(N/(i+1)))
		elif i%int(0.02*N)==0:
			time_array.append(time.time()-StartTime)
			TimeLeft = '{0:1.1f}'.format(float(interpolate.interp1d(np.arange(len(time_array)),time_array,fill_value='extrapolate')(49))-time_array[-1])
		print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete, ' + '{0:1.1f}'.format(time.time() - StartTime) \
			+ 'sec, (est. ' + TimeLeft,' sec left)		\r', end='')
	else:
		print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete           \r',end = '')

"""
################################
######## Tension Driven ########
################################

x_1 &= \theta \\
x_2 &= \dot{\theta} \\
u_1 &= T_{1} \\
u_2 &= T_{2} \\

################################
#### Muscle Velocity Driven ####
################################

x_1 &= \theta \\
x_2 &= \dot{\theta} \\
x_3 &= T_{1} \\
x_4 &= T_{2} \\
u_1 &= l_{m,1} \\
u_2 &= l_{m,2} \\

"""

N = 20001
Time = np.linspace(0,20,N)
dt = Time[1]-Time[0]

AllMuscleSettings = return_muscle_settings()

g,L = 9.80, 0.45 #m/s², m
M = 2 # kg
α1 = 0 # 10*np.pi/180 # rads
α2 = 0 # 10*np.pi/180 # rads
cT = 27.8
kT = 0.0047
LrT = 0.964
lo1 = AllMuscleSettings["BIC"]["Optimal Muscle Length"]/1000
lo2 = AllMuscleSettings["TRI"]["Optimal Muscle Length"]/1000

######### NEED OPTIMAL TENDON LENGTHS FOR BIC/TRI #########
lTo1 = (1)*AllMuscleSettings["BIC"]["Optimal Muscle Length"]/1000
lTo2 = (1)*AllMuscleSettings["TRI"]["Optimal Muscle Length"]/1000
###########################################################

[[r1,r2],[dr1,dr2],_]= return_MA_matrix_functions(AllMuscleSettings)
PCSA1 = 30**2*np.pi # mm²
PCSA2 = 30**2*np.pi # mm²
F_MAX1 = 0.25*PCSA1
F_MAX2 = 0.25*PCSA2

Amp = 7.5*np.pi/180
Base = 90*np.pi/180
Freq = 2*np.pi

k1,k2,k3,k4 = 100,100,10,100

MaxStep_Tension = 0.01 # percentage of positive maximum.
Tension_Bounds = [[0,F_MAX1],[0,0.10*F_MAX2]]

MaxStep_MuscleVelocity = 5 # percentage of positive maximum.
MuscleVelocity_Bounds =[[-2*lo1,2*lo1],[-0.2*lo2,0.2*lo2]]


"""
c_{1} &= -\frac{3g}{2L} \\
c_{2} &= \frac{3}{ML^2} \\

"""

c1 = -(3*g)/(2*L)
c2 = 3/(M*L**2)
c3 = np.cos(α1)
c4 = np.cos(α2)


'''
g_{1} &= r_{1}\left(\theta,\theta_{PS}=\frac{\pi}{2}\right) \\
g_{2} &= r_{2}\left(\theta,\theta_{PS}=\frac{\pi}{2}\right)  \\
g_{3} &= \frac{F_{\text{MAX},1}c^{T}}{l_{T,o,1}}\left(1 - \exp{\left(\frac{-T_1}{F_{\text{MAX},1}c^{T}k^{T}}\right)}\right) \hspace{1em} \text{(Variable stiffness coefficient from Zajac (1989) ODE for tendon force)} \\
g_{4} &= \text{sgn}\left(-r_1(\theta)\right)\cdot\dot{\theta}\cdot\sqrt{\left(\frac{\partial r_1}{\partial\theta}\right)^2 + r_1^2(\theta)} \\
g_{5} &= \frac{F_{\text{MAX},2}c^{T}}{l_{T,o,2}}\left(1 - \exp{\left(\frac{-T_2}{F_{\text{MAX},2}c^{T}k^{T}}\right)}\right) \hspace{1em} \text{(Variable stiffness coefficient from Zajac (1989) ODE for tendon force)} \\
g_{6} &= \text{sgn}\left(-r_2(\theta)\right)\cdot\dot{\theta}\cdot\sqrt{\left(\frac{\partial r_2}{\partial\theta}\right)^2 + r_2^2(\theta)} \\

'''

r1 = lambdify([θ_EFE],r1.subs([(θ_PS,np.pi/2)]))
r2 = lambdify([θ_EFE],r2.subs([(θ_PS,np.pi/2)]))
dr1_dθ = lambdify([θ_EFE],dr1.subs([(θ_PS,np.pi/2)]))
dr2_dθ = lambdify([θ_EFE],dr2.subs([(θ_PS,np.pi/2)]))

def R1(X):
	return(r1(X[0])) #
def dR1_dx1(X):
	return(dr1_dθ(X[0]))
def R2(X):
	return(r2(X[0])) #
def dR2_dx1(X):
	return(dr2_dθ(X[0]))
def KT_1(X):
	return((F_MAX1*cT/lTo1)*(1-np.exp(-X[2]/(F_MAX1*cT*kT)))) # NOT NORMALIZED (in N/m)
def v_MTU1(X):
	return(np.sign(-R1(X))*X[1]*np.sqrt(dR1_dx1(X)**2 + R1(X)**2)) # NOT NORMALIZED (in m/s)
def KT_2(X):
	return((F_MAX2*cT/lTo2)*(1-np.exp(-X[3]/(F_MAX2*cT*kT)))) # NOT NORMALIZED (in N/m)
def v_MTU2(X):
	return(np.sign(-R2(X))*X[1]*np.sqrt(dR2_dx1(X)**2 + R2(X)**2)) # NOT NORMALIZED (in m/s)

"""
################################
######## Tension Driven ########
################################

\dot{x}_1 &= x_{2} \\
\dot{x}_2 &= c_{1}\sin(x_{1}) + c_{2}R_{1}u_{1} - c_{2}R_{2}u_{2} \\
u_1 &= T_{1} \\
u_2 &= T_{2} \\

################################
#### Muscle Velocity Driven ####
################################

\dot{x}_1 &= x_{2} \\
\dot{x}_2 &= c_{1}\sin(x_{1}) + c_{2}R_{1}x_{3} - c_{2}R_{2}x_{4} \\
\dot{x}_3 &= K_{T,1}(v_{MTU,1} - c_{3}u_1) \\
\dot{x}_4 &= K_{T,2}(v_{MTU,2} - c_{4}u_2) \\
u_1 &= \dot{l}_{m,1} \\
u_2 &= \dot{l}_{m,2} \\

"""

def dX1_dt(X):
	return(X[1])
def dX2_dt(X,U=None):
	if U==None:
		return(c1*np.sin(X[0]) + c2*R1(X)*X[2] + c2*R2(X)*X[3])
	else:
		return(c1*np.sin(X[0]) + c2*R1(X)*U[0] + c2*R2(X)*U[1])
def dX3_dt(X,U=None):
	if U == None:
		return(KT_1(X)*(v_MTU1(X) - c3*X[6]))
	else:
		return(KT_1(X)*(v_MTU1(X) - c3*U[0]))
def dX4_dt(X,U=None):
	if U == None:
		return(KT_2(X)*(v_MTU2(X) - c4*X[7]))
	else:
		return(KT_2(X)*(v_MTU2(X) - c4*U[1]))

r = lambda t: Amp*np.sin(Freq*t) + Base
dr = lambda t: Amp*Freq*np.cos(Freq*t)
d2r = lambda t: -Amp*Freq**2*np.sin(Freq*t)
d3r = lambda t: -Amp*Freq**3*np.cos(Freq*t)

def return_constraint_variables_tension_driven(t,X):
	def Z1(t,X):
		return(r(t) - X[0])
	def dZ1(t,X):
		return(dr(t) - dX1_dt(X))
	def A1(t,X):
		return(dr(t) + k1*Z1(t,X))
	def dA1(t,X):
		return(d2r(t) + k1*dZ1(t,X))
	def Z2(t,X):
		return(X[1] - A1(t,X))
	"""
	def dZ2(t,X,U):
		return(c1*np.sin(X[0]) + c2*R1(X)*U[0] + c2*R2(X)*U[1] - dA1(t,X))
	"""
	def A2(t,X):
		return(Z1(t,X) + dA1(t,X) - c1*np.sin(X[0]) - k2*Z2(t,X))
	Coefficient1 = c2*R1(X)
	Coefficient2 = c2*R2(X)
	Constraint = A2(t,X)
	return(Coefficient1,Coefficient2,Constraint)
def return_constraint_variables_muscle_velocity_driven(t,X):
	def Z1(t,X):
		return(r(t) - X[0])
	def dZ1(t,X):
		return(dr(t) - dX1_dt(X))
	def d2Z1(t,X):
		return(d2r(t) - dX2_dt(X))
	def A1(t,X):
		return(dr(t) + k1*Z1(t,X))
	def dA1(t,X):
		return(d2r(t) + k1*dZ1(t,X))
	def d2A1(t,X):
		return(d3r(t) + k1*d2Z1(t,X))
	def Z2(t,X):
		return(X[1] - A1(t,X))
	def dZ2(t,X):
		return(dX2_dt(X) - dA1(t,X))
	def A2(t,X):
		return(Z1(t,X) + dA1(t,X) - c1*np.sin(X[0]) - k2*Z2(t,X))
	def dA2(t,X):
		return(dZ1(t,X) + d2A1(t,X) - c1*np.cos(X[0])*dX1_dt(X) - k2*dZ2(t,X))
	def Z3(t,X):
		return(c2*R1(X)*X[2] + c2*R2(X)*X[3] - A2(t,X))
	def dZ3(t,X):
		g1 = R1(X)
		g2 = R2(X)
		g3 = KT_1(X)
		g5 = KT_2(X)
		return(c2*dR1_dx1(X)*X[1]*X[2] + c2*dR2_dx1(X)*X[1]*X[3] \
						+ c2*g1*g3*v_MTU1(X) - c2*c3*g1*g3*U[0] \
							+ c2*g2*g5*v_MTU2(X) - c2*c4*g2*g5*U[1] \
								- dA2(t,X))
	def A3(t,X):
		return(Z2(t,X) - dA2(t,X) + k3*Z3(t,X) \
		+ c2*dR1_dx1(X)*X[1]*X[2] + 	c2*dR2_dx1(X)*X[1]*X[3] \
				+ c2*R1(X)*KT_1(X)*v_MTU1(X) + c2*R2(X)*KT_2(X)*v_MTU2(X))
	Coefficient1 = c2*c3*R1(X)*KT_1(X)
	Coefficient2 = c2*c4*R2(X)*KT_2(X)
	Constraint = A3(t,X)
	return(Coefficient1,Coefficient2,Constraint)
def return_U_tension_driven(t,X,U,dt,MaxStep,Bounds,Noise):
	import random
	Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_tension_driven(t,X)
	# elif Method == "Muscle Velocity":
	# 		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t,X)
	assert np.shape(Bounds)==(2,2), "Bounds must be (2,2)."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."
	if Constraint1 != 0:
		assert Coefficient1!=0 and Coefficient2!=0, "Error with Coefficients. Shouldn't be zero with nonzero constraint."
	else:
		assert Coefficient1!=0 and Coefficient2!=0, "Error with Constraint. 0 = 0 implies all inputs valid."

	if Coefficient1 == 0:
		LowerBound = Bounds[0][0]
		UpperBound = Bounds[0][1]
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
	elif Coefficient2 == 0:
		LowerBound = Constraint1/Coefficient1
		UpperBound = Constraint1/Coefficient1
		FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
		FeasibleInput2 = (Bounds[1][1]-Bounds[1][0])*np.random.rand(1000) + Bounds[1][0]
	elif np.sign(-Coefficient1) == np.sign(Coefficient2):
		if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
			LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
		else:
			LowerBound = Bounds[0][0]

		if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
			UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
		else:
			UpperBound = Bounds[0][1]
		assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in FeasibleInput1])
	else: # np.sign(-Coefficient1) != np.sign(Coefficient2)
		if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
			LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
		else:
			LowerBound = Bounds[0][0]

		if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
			UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
		else:
			UpperBound = Bounds[0][1]
		assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in FeasibleInput1])
	#
	# if Constraint1 > 0:
	# 	assert Bounds[0] >= Constraint1/Coefficient1, "Tension cannot be made by muscle 1."
	# 	LowerBound = Constraint1/Coefficient1
	# 	if (Constraint1 - Coefficient1*Bounds[0])/Coefficient2 > Bounds[1]:
	# 		UpperBound = (Constraint1 - Coefficient2*Bounds[1])/Coefficient1
	# 	else:
	# 		UpperBound = Bounds[0]
	# elif Constraint1 < 0:
	# 	assert Bounds[1] >= Constraint1/Coefficient2, "Tension cannot be made by muscle 2."
	# 	LowerBound = 0
	# 	if (Constraint1 - Coefficient2*Bounds[1])/Coefficient1 > Bounds[0]:
	# 		UpperBound = Bounds[0]
	# 	else:
	# 		UpperBound = (Constraint1 - Coefficient2*Bounds[1])/Coefficient1
	# else: # Constraint1 == 0
	# 	LowerBound = 0
	# 	if -Coefficient1*Bounds[0]/Coefficient2 > Bounds[1]:
	# 		UpperBound = -Coefficient2*Bounds[1]/Coefficient1
	# 	else:
	# 		UpperBound = Bounds[0]
	"""
	Checking to see which inputs have the appropriate allowable step size.
	"""
	euclid_dist = np.array(list(map(lambda u1,u2: np.sqrt(((U[0]-u1)/Bounds[0][1])**2 + ((U[1]-u2)/Bounds[1][1])**2),\
							FeasibleInput1,FeasibleInput2)))

	if t<10*dt: MaxStep = 10*MaxStep
	feasible_index = np.where(euclid_dist<=MaxStep)
	# elif Method == "Muscle Velocity":
	# 	feasible_index = np.where(np.logical_and(np.logical_and(euclid_dist>=0, euclid_dist<=MaxStep),np.sign(FeasibleInput1)!=np.sign(FeasibleInput2)))
	if len(feasible_index[0]) == 0: import ipdb; ipdb.set_trace()
	next_index = random.choice(feasible_index[0])
	u1 = FeasibleInput1[next_index]
	u2 = FeasibleInput2[next_index]
	return([u1,u2])
def return_U_muscle_velocity_driven(t,X,U,dt,MaxStep,Bounds,Noise):
	"""
	Enforcing a hyperbolic domain constraint to allow for realistic lengthening/shortenting relationships.
	Input2 = (lo1*0.001)*(lo2*0.001)/Input1 = lo1*lo2/(10^6*Input1)
	"""

	import random
	Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t,X)
	assert np.shape(Bounds)==(2,2), "Bounds must be (2,2)."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."
	if Constraint1 != 0:
		assert Coefficient1!=0 or Coefficient2!=0, "Error with Coefficients. Shouldn't be zero with nonzero constraint."
	else:
		assert Coefficient1!=0 or Coefficient2!=0, "Error with Constraint. 0 = 0 implies all inputs valid."

	if abs(Coefficient1) <= 1e-7:
		LowerBound = Bounds[0][0]
		UpperBound = Bounds[0][1]
		if Constraint1/Coefficient2 > 0:
			LowerBound = Bounds[0][0]
			UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
		else:
			LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
			UpperBound = Bounds[0][1]
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
	elif abs(Coefficient2) <= 1e-7:
		LowerBound = Constraint1/Coefficient1
		UpperBound = Constraint1/Coefficient1
		FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
		if Constraint1/Coefficient1 < 0:
			LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
			UpperBound = Bounds[1][1]
		else:
			LowerBound = Bounds[1][0]
			UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
		FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
	elif np.sign(-Coefficient1) == np.sign(Coefficient2): # DIFFERENT SIGNS
		if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
			LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
		else:
			LowerBound = Bounds[0][0]

		if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
			UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
		else:
			UpperBound = Bounds[0][1]
		assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
		HyperbolicBounds = np.sort([(Constraint1 - \
										np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
											/(2*Coefficient1), \
								 	(Constraint1 + \
										np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
											/(2*Coefficient1)])
		LowerBound = max([LowerBound,HyperbolicBounds[0]])
		UpperBound = min([UpperBound,HyperbolicBounds[1]])
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in FeasibleInput1])
	else: # np.sign(-Coefficient1) != np.sign(Coefficient2) SAME SIGN
		if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
			LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
		else:
			LowerBound = Bounds[0][0]

		if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
			UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
		else:
			UpperBound = Bounds[0][1]
		assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
		if Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001) > 0:
			HyperbolicBounds = np.sort([(Constraint1 - \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1), \
									 	(Constraint1 + \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1)])

			assert LowerBound < HyperbolicBounds[0] or HyperbolicBounds[1] < UpperBound, "No feasible solutions."

			FeasibleInput1 = []
			while len(FeasibleInput1)<1000:
				Random1 = (UpperBound-LowerBound)*np.random.rand() + LowerBound
				if Random1<HyperbolicBounds[0] or Random1>HyperbolicBounds[1]: FeasibleInput1.append(Random1)
			FeasinbleInput1 = np.array(FeasibleInput1)
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		else:
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
	def plot_constraints():
		import matplotlib.pyplot as plt
		plt.figure()
		Input1 = np.linspace(LowerBound,UpperBound,1001)
		Input2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in Input1])
		plt.plot(Input1,Input2,'k--')
		plt.plot(Input1, (lo1*0.001)*(lo2*0.001)/Input1,'r')
		plt.scatter(FeasibleInput1,FeasibleInput2,c='g',marker = '.')
		plt.ylim(MuscleVelocity_Bounds[1])
		plt.show()

	"""
	Checking to see which inputs have the appropriate allowable step size. In normalized muscle velocity.
	"""
	euclid_dist = np.array(list(map(lambda u1,u2: np.sqrt(((U[0]-u1)/lo1)**2 + ((U[1]-u2)/lo2)**2),\
							FeasibleInput1,FeasibleInput2)))

	if t<30*dt: MaxStep = 10*MaxStep
	feasible_index = np.where(euclid_dist<=MaxStep)

	if len(feasible_index[0]) == 0: import ipdb; ipdb.set_trace()
	next_index = random.choice(feasible_index[0])
	u1 = FeasibleInput1[next_index]
	u2 = FeasibleInput2[next_index]
	return([u1,u2])
def plot_MA_values(Time,x1):
	import matplotlib.pyplot as plt
	import numpy as np

	fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8,6))
	plt.subplots_adjust(hspace=0.1,bottom=0.1)

	plt.suptitle("Moment arm equations")
	ax1.plot(np.linspace(0,np.pi*(160/180),1001),\
				np.array(list(map(lambda x1: R1([x1]),np.linspace(0,np.pi*(160/180),1001)))),\
				'0.70')
	ax1.plot(np.linspace(min(x1),max(x1),101),\
				np.array(list(map(lambda x1: R1([x1]),np.linspace(min(x1),max(x1),101)))),\
				'g',lw=3)
	ax1.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
	ax1.set_xticklabels([""]*len(ax1.get_xticks()))
	ax1.set_ylabel("Moment Arm for\n Muscle 1 (m)")

	ax2.plot(Time,np.array(list(map(lambda x1: R1([x1]),x1))),'g')
	ax2.set_ylim(ax1.get_ylim())
	ax2.set_yticks(ax1.get_yticks())
	ax2.set_yticklabels([""]*len(ax1.get_yticks()))
	ax2.set_xticklabels([""]*len(ax2.get_xticks()))

	ax3.plot(np.linspace(0,np.pi*(160/180),1001),\
				np.array(list(map(lambda x1: R2([x1]),np.linspace(0,np.pi*(160/180),1001)))),\
				'0.70')
	ax3.plot(np.linspace(min(x1),max(x1),101),\
				np.array(list(map(lambda x1: R2([x1]),np.linspace(min(x1),max(x1),101)))),\
				'r',lw=3)
	ax3.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
	ax3.set_xticklabels([r"$0$",r"$\frac{\pi}{4}$",r"$\frac{\pi}{2}$",r"$\frac{3\pi}{4}$",r"$\pi$"])
	ax3.set_xlabel("Joint Angle (rads)")
	ax3.set_ylabel("Moment Arm for\n Muscle 2 (m)")

	ax4.plot(Time,np.array(list(map(lambda x1: R2([x1]),x1))),'r')
	ax4.set_ylim(ax3.get_ylim())
	ax4.set_yticks(ax3.get_yticks())
	ax4.set_yticklabels([""]*len(ax3.get_yticks()))
	ax4.set_xlabel("Time (s)")
	return(fig,[ax1,ax2,ax3,ax4])
def animate_test_2(response,t,x1,x2,x3,x4,u1,u2,dt,Bounds):
	assert type(response)==bool, "Input must be either True or False."

	if response == True:
		import numpy as np
		import matplotlib.pyplot as plt
		import matplotlib.animation as animation
		import time

		fig = plt.figure(figsize=(10,8))
		ax1 = plt.gca()

		DescriptiveTitle = "Plotting Constraints vs. Time"

		ax1.set_title(DescriptiveTitle,Fontsize=20,y=0.975)

		#Hyperbolic Constraint
		Input1 = list(np.linspace(Bounds[0][0],Bounds[0][1],1000001))
		Input1.remove(0)
		Input1 = np.array(Input1)
		ax1.plot(Input1,lo1*lo2*0.001**2/Input1,'r',lw=2)
		Input2 = list(np.linspace(Bounds[0][0],Bounds[0][1],1000001))
		Input2.remove(0)
		Input2 = np.array(Input2)
		ax1.plot(lo1*lo2*0.001**2/Input2,Input2,'r',lw=2)
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[0],[x1[0],x2[0],x3[0],x4[0]])
		if abs(Coefficient1) <= 1e-7:
			LowerBound = Bounds[0][0]
			UpperBound = Bounds[0][1]
			if Constraint1/Coefficient2 > 0:
				LowerBound = Bounds[0][0]
				UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
			else:
				LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
				UpperBound = Bounds[0][1]
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
		elif abs(Coefficient2) <= 1e-7:
			LowerBound = Constraint1/Coefficient1
			UpperBound = Constraint1/Coefficient1
			FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
			if Constraint1/Coefficient1 < 0:
				LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
				UpperBound = Bounds[1][1]
			else:
				LowerBound = Bounds[1][0]
				UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
			FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		elif np.sign(-Coefficient1) == np.sign(Coefficient2): # DIFFERENT SIGNS
			if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
				LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
			else:
				LowerBound = Bounds[0][0]

			if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
				UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
			else:
				UpperBound = Bounds[0][1]
			assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
			HyperbolicBounds = np.sort([(Constraint1 - \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1), \
									 	(Constraint1 + \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1)])
			LowerBound = max([LowerBound,HyperbolicBounds[0]])
			UpperBound = min([UpperBound,HyperbolicBounds[1]])
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		else: # np.sign(-Coefficient1) != np.sign(Coefficient2) SAME SIGN
			if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
				LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
			else:
				LowerBound = Bounds[0][0]

			if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
				UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
			else:
				UpperBound = Bounds[0][1]
			assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
			if Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001) > 0:
				HyperbolicBounds = np.sort([(Constraint1 - \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1), \
										 	(Constraint1 + \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1)])

				assert LowerBound < HyperbolicBounds[0] or HyperbolicBounds[1] < UpperBound, "No feasible solutions."

				FeasibleInput1 = []
				while len(FeasibleInput1)<1000:
					Random1 = (UpperBound-LowerBound)*np.random.rand() + LowerBound
					if Random1<HyperbolicBounds[0] or Random1>HyperbolicBounds[1]: FeasibleInput1.append(Random1)
				FeasinbleInput1 = np.array(FeasibleInput1)
				FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
										for el in FeasibleInput1])
			else:
				FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
				FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
										for el in FeasibleInput1])
		cline, = plt.plot(FeasibleInput1,FeasibleInput2,'b',lw=2)
		TimeText = plt.text(0.1,0.1,"t = " + str(t[0]),fontsize=16)
		chosenpoint, = plt.plot(u1[0],u2[0],c='k',marker='o')
		ax1.set_xlim(Bounds[0])
		ax1.set_ylim(Bounds[1])
		ax1.spines['right'].set_visible(False)
		ax1.spines['top'].set_visible(False)
		ax1.set_aspect('equal')

		def animate(i):
			Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[i],[x1[i],x2[i],x3[i],x4[i]])
			if abs(Coefficient1) <= 1e-7:
				LowerBound = Bounds[0][0]
				UpperBound = Bounds[0][1]
				if Constraint1/Coefficient2 > 0:
					LowerBound = Bounds[0][0]
					UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
				else:
					LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
					UpperBound = Bounds[0][1]
				FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
				FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
			elif abs(Coefficient2) <= 1e-7:
				LowerBound = Constraint1/Coefficient1
				UpperBound = Constraint1/Coefficient1
				FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
				if Constraint1/Coefficient1 < 0:
					LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
					UpperBound = Bounds[1][1]
				else:
					LowerBound = Bounds[1][0]
					UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
				FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			elif np.sign(-Coefficient1) == np.sign(Coefficient2): # DIFFERENT SIGNS
				if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
					LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
				else:
					LowerBound = Bounds[0][0]

				if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
					UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
				else:
					UpperBound = Bounds[0][1]
				assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
				HyperbolicBounds = np.sort([(Constraint1 - \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1), \
										 	(Constraint1 + \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1)])
				LowerBound = max([LowerBound,HyperbolicBounds[0]])
				UpperBound = min([UpperBound,HyperbolicBounds[1]])
				FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
				FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
										for el in FeasibleInput1])
			else: # np.sign(-Coefficient1) != np.sign(Coefficient2) SAME SIGN
				if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
					LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
				else:
					LowerBound = Bounds[0][0]

				if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
					UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
				else:
					UpperBound = Bounds[0][1]
				assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
				if Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001) > 0:
					HyperbolicBounds = np.sort([(Constraint1 - \
													np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
														/(2*Coefficient1), \
											 	(Constraint1 + \
													np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
														/(2*Coefficient1)])

					assert LowerBound < HyperbolicBounds[0] or HyperbolicBounds[1] < UpperBound, "No feasible solutions."

					FeasibleInput1 = []
					while len(FeasibleInput1)<1000:
						Random1 = (UpperBound-LowerBound)*np.random.rand() + LowerBound
						if Random1<HyperbolicBounds[0] or Random1>HyperbolicBounds[1]: FeasibleInput1.append(Random1)
					FeasinbleInput1 = np.array(FeasibleInput1)
					FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
											for el in FeasibleInput1])
				else:
					FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
					FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
											for el in FeasibleInput1])
			cline.set_xdata(FeasibleInput1)
			cline.set_ydata(FeasibleInput2)
			chosenpoint.set_xdata(u1[i])
			chosenpoint.set_ydata(u2[i])
			TimeText.set_text("t = " + str(t[i]))
			return cline,chosenpoint,TimeText,


		# Init only required for blitting to give a clean slate.
		def init():
			ax1.plot(Input1,lo1*lo2*0.001**2/Input1,'r',lw=2)
			Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[0],[x1[0],x2[0],x3[0],x4[0]])
			if abs(Coefficient1) <= 1e-7:
				LowerBound = Bounds[0][0]
				UpperBound = Bounds[0][1]
				if Constraint1/Coefficient2 > 0:
					LowerBound = Bounds[0][0]
					UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
				else:
					LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
					UpperBound = Bounds[0][1]
				FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
				FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
			elif abs(Coefficient2) <= 1e-7:
				LowerBound = Constraint1/Coefficient1
				UpperBound = Constraint1/Coefficient1
				FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
				if Constraint1/Coefficient1 < 0:
					LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
					UpperBound = Bounds[1][1]
				else:
					LowerBound = Bounds[1][0]
					UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
				FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			elif np.sign(-Coefficient1) == np.sign(Coefficient2): # DIFFERENT SIGNS
				if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
					LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
				else:
					LowerBound = Bounds[0][0]

				if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
					UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
				else:
					UpperBound = Bounds[0][1]
				assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
				HyperbolicBounds = np.sort([(Constraint1 - \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1), \
										 	(Constraint1 + \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1)])
				LowerBound = max([LowerBound,HyperbolicBounds[0]])
				UpperBound = min([UpperBound,HyperbolicBounds[1]])
				FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
				FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
										for el in FeasibleInput1])
			else: # np.sign(-Coefficient1) != np.sign(Coefficient2) SAME SIGN
				if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
					LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
				else:
					LowerBound = Bounds[0][0]

				if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
					UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
				else:
					UpperBound = Bounds[0][1]
				assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
				if Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001) > 0:
					HyperbolicBounds = np.sort([(Constraint1 - \
													np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
														/(2*Coefficient1), \
											 	(Constraint1 + \
													np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
														/(2*Coefficient1)])

					assert LowerBound < HyperbolicBounds[0] or HyperbolicBounds[1] < UpperBound, "No feasible solutions."

					FeasibleInput1 = []
					while len(FeasibleInput1)<1000:
						Random1 = (UpperBound-LowerBound)*np.random.rand() + LowerBound
						if Random1<HyperbolicBounds[0] or Random1>HyperbolicBounds[1]: FeasibleInput1.append(Random1)
					FeasinbleInput1 = np.array(FeasibleInput1)
					FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
											for el in FeasibleInput1])
				else:
					FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
					FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
											for el in FeasibleInput1])
			cline, = plt.plot(FeasibleInput1,FeasibleInput2,'b',lw=2)
			cline.set_visible(False)
			chosenpoint, = plt.plot(u1[0],u2[0],c='k',marker='o')
			chosenpoint.set_visible(False)
			TimeText = plt.text(0.75,0.75,"t = " + str(t[0]),fontsize=16)
			TimeText.set_visible(False)
			return cline,chosenpoint,TimeText,

		ani = animation.FuncAnimation(fig, animate, np.arange(1, len(t),1), init_func=init,interval=1, blit=False)
		plt.show()
def plot_individual_constraint_versus_time_test_2(t,x1,x2,x3,x4,Return = False):
	import numpy as np
	import matplotlib.pyplot as plt

	DescriptiveTitle = "Plotting Coefficients/Constraints vs. Time"
	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.975)

	"""
	A⋅u₁ + B⋅u₂ = C
	"""

	A,B,C = [],[],[]
	for i in range(len(x1)):
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[i],[x1[i],x2[i],x3[i],x4[i]])
		A.append(Coefficient1)
		B.append(Coefficient2)
		C.append(Constraint1)

	ax1.plot(t[:len(x1)],A,'r',lw=2)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$1^{st}$ Coefficient")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:len(x1)],B,'b',lw=2)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$2^{nd}$ Coefficient")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	ax3.plot(t[:len(x1)],C,'k',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel("Constraint")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_individual_coefficient2_versus_time_test_2(t,x1,x2,x3,x4,Return = False):
	import numpy as np
	import matplotlib.pyplot as plt

	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
	plt.subplots_adjust(top=0.9,hspace=0.4,bottom=0.1,left=0.075,right=0.975)
	plt.suptitle(r"Plotting $2^{nd}$ Coefficient vs. Time",Fontsize=20,y=0.975)

	"""
	B = c2⋅c4⋅R2(X)⋅KT_2(X)
	"""

	r2,kt_2,B = [],[],[]
	for i in range(len(x1)):
		_,Coefficient2,_ = return_constraint_variables_muscle_velocity_driven(t[i],[x1[i],x2[i],x3[i],x4[i]])
		B.append(Coefficient2)
		r2.append(R2([x1[i],x2[i],x3[i],x4[i]]))
		kt_2.append(G5([x1[i],x2[i],x3[i],x4[i]]))

	ax1.plot(t[:len(x1)],r2,'b--',lw=2)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$g_{2}(\vec{x}(t))$")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:len(x1)],kt_2,'b:',lw=2)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$g_{5}(\vec{x}(t))$")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	ax3.plot(t[:len(x1)],B,'b',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel(r"$2^{nd}$ Coefficient")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_individual_coefficient1_versus_time_test_2(t,x1,x2,x3,x4,Return = False):
	import numpy as np
	import matplotlib.pyplot as plt

	DescriptiveTitle = "Plotting Coefficients/Constraints vs. Time"
	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
	plt.subplots_adjust(top=0.9,hspace=0.4,bottom=0.1,left=0.075,right=0.975)
	plt.suptitle(r"Plotting $1^{st}$ Coefficient vs. Time",Fontsize=20,y=0.975)

	"""
	A = c2⋅c3⋅R1(X)⋅KT_1(X)
	"""

	r1,kt_1,B = [],[],[]
	for i in range(len(x1)):
		Coefficient1,_,_ = return_constraint_variables_muscle_velocity_driven(t[i],[x1[i],x2[i],x3[i],x4[i]])
		B.append(Coefficient1)
		r1.append(R1([x1[i],x2[i],x3[i],x4[i]]))
		kt_1.append(KT_1([x1[i],x2[i],x3[i],x4[i]]))

	ax1.plot(t[:len(x1)],r1,'r--',lw=2)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$g_{2}(\vec{x}(t))$")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:len(x1)],kt_1,'r:',lw=2)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$g_{5}(\vec{x}(t))$")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	ax3.plot(t[:len(x1)],B,'r',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel(r"$2^{nd}$ Coefficient")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_states(t,X,Return=False):
	import numpy as np
	import matplotlib.pyplot as plt

	NumStates = np.shape(X)[0]
	NumRows = int(np.ceil(NumStates/5))
	if NumStates < 5:
		NumColumns = NumStates
	else:
		NumColumns = 5

	ColumnNumber = [el%5 for el in np.arange(0,NumStates,1)]
	RowNumber = [int(el/5) for el in np.arange(0,NumStates,1)]

	DescriptiveTitle = "Plotting States vs. Time"
	fig, axes = plt.subplots(NumRows,NumColumns,figsize=(3*NumColumns,2*NumRows + 2))
	plt.subplots_adjust(top=0.85,hspace=0.4,bottom=0.15,left=0.075,right=0.975)
	plt.suptitle(r"Plotting $1^{st}$ Coefficient vs. Time",Fontsize=20,y=0.975)
	if NumStates <=5:
		for j in range(NumStates):
			axes[ColumnNumber[j]].spines['right'].set_visible(False)
			axes[ColumnNumber[j]].spines['top'].set_visible(False)
			axes[ColumnNumber[j]].plot(t[:np.shape(X)[1]],X[j])
			if ColumnNumber[j]!=0:
				axes[ColumnNumber[j]].set_xticklabels(\
									[""]*len(axes[ColumnNumber[j]].get_xticks()))
			else:
				axes[ColumnNumber[j]].set_xlabel("Time (s)")
			axes[ColumnNumber[j]].set_title(r"$x_{" + str(j+1) + "}$")

	else:
		for j in range(NumStates):
			axes[RowNumber[j],ColumnNumber[j]].spines['right'].set_visible(False)
			axes[RowNumber[j],ColumnNumber[j]].spines['top'].set_visible(False)
			axes[RowNumber[j],ColumnNumber[j]].plot(t[:np.shape(X)[1]],X[j])
			if not(RowNumber[j] == RowNumber[-1] and ColumnNumber[j]==0):
				axes[RowNumber[j],ColumnNumber[j]].set_xticklabels(\
									[""]*len(axes[RowNumber[j],ColumnNumber[j]].get_xticks()))
			else:
				axes[RowNumber[j],ColumnNumber[j]].set_xlabel("Time (s)")
			axes[RowNumber[j],ColumnNumber[j]].set_title(r"$x_{" + str(j+1) + "}$")
		if NumStates%5!=0:
			[fig.delaxes(axes[RowNumber[-1],el]) for el in range(ColumnNumber[-1]+1,5)]

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_inputs(t,u1,u2,Return=False):
	import numpy as np
	import matplotlib.pyplot as plt

	DescriptiveTitle = "Plotting Inputs vs. Time"
	fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,5))
	plt.subplots_adjust(top=0.9,hspace=0.4,bottom=0.1,left=0.075,right=0.975)
	plt.suptitle("Plotting Inputs vs. Time",Fontsize=20,y=0.975)

	ax1.plot(t[:len(u1)],u1,'g--',lw=2)
	ax1.plot([-1,t[len(u1)]+1],[0,0],'k--',lw=0.5)
	ax1.set_xlim([t[0],t[len(u1)]])
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$u_1$")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:len(u2)],u2,'g',lw=2)
	ax2.plot([-1,t[len(u2)]+1],[0,0],'k--',lw=0.5)
	ax2.set_xlim([t[0],t[len(u1)]])
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$u_2$")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
x1_1,x2_1 = [Base],[Amp*Freq]
u1_1,u2_1 = [100],[10]

x1_2,x2_2,x3_2,x4_2= [Base],[Amp*Freq],[100],[70]
u1_2,u2_2 = [0.01],[0.01]


CocontractionIndex = 2

AddNoise = False
if AddNoise == True:
    np.random.seed(seed=None)
    NoiseArray = np.random.normal(loc=0.0,scale=0.2,size=(2,len(Time)))
else:
    NoiseArray = np.zeros((2,len(Time)))

def update_policy_tension_driven(t,x1_1,x2_1,dt,NoiseArray):
	import numpy as np
	Method = "Tension"
	X = [x1_1[-1],x2_1[-1]]
	U = [u1_1[-1],u2_1[-1]]
	U = return_U_tension_driven(t,X,U,dt,MaxStep_Tension,Tension_Bounds,NoiseArray[:,int(t/dt)])
	u1_1.append(U[0])
	u2_1.append(U[1])
	x2_1.append(x2_1[-1] + dX2_dt(X,U=U)*dt)
	x1_1.append(x1_1[-1] + dX1_dt(X)*dt)
def update_policy_muscle_velocity_driven(t,x1_2,x2_2,x3_2,x4_2,dt,NoiseArray):
	import numpy as np
	Method = "Muscle Velocity"
	X = [x1_2[-1],x2_2[-1],x3_2[-1],x4_2[-1]]
	U = [u1_2[-1],u2_2[-1]]
	U = return_U_muscle_velocity_driven(t,X,U,dt,MaxStep_MuscleVelocity,MuscleVelocity_Bounds,NoiseArray[:,int(t/dt)])
	u1_2.append(U[0])
	u2_2.append(U[1])
	x4_2.append(x4_2[-1] + dX4_dt(X,U=U)*dt)
	x3_2.append(x3_2[-1] + dX3_dt(X,U=U)*dt)
	x2_2.append(x2_2[-1] + dX2_dt(X)*dt)
	x1_2.append(x1_2[-1] + dX1_dt(X)*dt)

StartTime = time.time()
for t in Time[1:]:
	update_policy_tension_driven(t,x1_1,x2_1,dt,NoiseArray)
	update_policy_muscle_velocity_driven(t,x1_2,x2_2,x3_2,x4_2,dt,NoiseArray)
	statusbar(int(t/dt),len(Time),StartTime=StartTime,Title="Forced-Pendulum")

fig1,[ax1_1,ax2_1,ax3_1,ax4_1] = plot_MA_values(Time,x1_1)
fig2,[ax1_2,ax2_2,ax3_2,ax4_2] = plot_MA_values(Time,x1_2)

plt.figure()
plt.title("Underdetermined Tendon-Tension-Driven\nForced Pendulum Example",\
                fontsize=16,color='gray')
plt.plot(Time,x1_1,'b',lw=2)
plt.plot(Time,x1_2,'g',lw=2)
plt.plot(Time,r(Time),'r--')
plt.xlabel("Time (s)")
plt.ylabel("Desired Measure")
plt.legend([r"Output $y = x_{1}$ (Tension)",r"Output $y = x_{1}$ (mm Velocity)",r"Reference $r(t) = \frac{\pi}{24}\sin(2\pi t) + \frac{\pi}{2}$"],loc='best')

plt.figure()
plt.title('Error vs. Time')
plt.plot(Time, r(Time)-x1_1,color='b')
plt.plot(Time, r(Time)-x1_2,color='g')
plt.legend(["Tension Driven","Muscle Velocity Driven"],loc='best')
plt.xlabel("Time (s)")
plt.ylabel("Error")

plt.figure()
plt.plot(Time,u1_1,'g',Time,u2_1,'r')
plt.title('Tendon Tensions vs. Time')
plt.xlabel("Time (s)")
plt.ylabel("Tendon Tensions (N)")
plt.legend(["Muscle 1","Muscle 2"])

plt.figure()
plt.plot(Time,u1_2,'g',Time,u2_2,'r')
plt.title('Muscle Velocities vs. Time')
plt.xlabel("Time (s)")
plt.ylabel("Muscle Velocities (m/s)")
plt.legend(["Muscle 1","Muscle 2"])

plt.show()
