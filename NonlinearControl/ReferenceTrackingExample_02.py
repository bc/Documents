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
lTo1 = (0.2)*AllMuscleSettings["BIC"]["Optimal Muscle Length"]/1000
lTo2 = (0.2)*AllMuscleSettings["TRI"]["Optimal Muscle Length"]/1000
###########################################################

[[r1,r2],[dr1,dr2],_]= return_MA_matrix_functions(AllMuscleSettings)
PCSA1 = 30**2*np.pi # mm²
PCSA2 = 30**2*np.pi # mm²
F_MAX1 = 0.25*PCSA1
F_MAX2 = 0.25*PCSA2

Amp = 7.5*np.pi/180
Base = 90*np.pi/180
Freq = 2*np.pi

k1,k2,k3,k4 = 100,100,100,100

MaxStep_Tension = np.average([F_MAX1,F_MAX2])*0.01
Tension_Bounds = [[0,F_MAX1],[0,0.10*F_MAX2]]

MaxStep_MuscleVelocity = max([lo1,lo2])
MuscleVelocity_Bounds =[[-1*lo1,1*lo1],[-1*lo2,1*lo2]]


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

R1 = lambdify([θ_EFE],r1.subs([(θ_PS,np.pi/2)]))
R2 = lambdify([θ_EFE],r2.subs([(θ_PS,np.pi/2)]))
dR1_dθ = lambdify([θ_EFE],dr1.subs([(θ_PS,np.pi/2)]))
dR2_dθ = lambdify([θ_EFE],dr2.subs([(θ_PS,np.pi/2)]))

def G1(X):
	return(R1(X[0])) #
def dG1(X):
	return(dR1_dθ(X[0]))
def G2(X):
	return(R2(X[0])) #
def dG2(X):
	return(dR2_dθ(X[0]))
def G3(X):
	return((F_MAX1*cT/lTo1)*(1-np.exp(-X[2]/(F_MAX1*cT*kT)))) # NOT NORMALIZED (in N/m)
def G4(X):
	return(np.sign(-G1(X))*X[1]*np.sqrt(dG1(X)**2 + G1(X)**2)) # NOT NORMALIZED (in m/s)
def G5(X):
	return((F_MAX2*cT/lTo2)*(1-np.exp(-X[3]/(F_MAX2*cT*kT)))) # NOT NORMALIZED (in N/m)
def G6(X):
	return(np.sign(-G2(X))*X[1]*np.sqrt(dG2(X)**2 + G2(X)**2)) # NOT NORMALIZED (in m/s)

"""
################################
######## Tension Driven ########
################################

\dot{x}_1 &= x_{2} \\
\dot{x}_2 &= c_{1}\sin(x_{1}) + c_{2}g_{1}u_{1} - c_{2}g_{2}u_{2} \\
u_1 &= T_{1} \\
u_2 &= T_{2} \\

################################
#### Muscle Velocity Driven ####
################################

\dot{x}_1 &= x_{2} \\
\dot{x}_2 &= c_{1}\sin(x_{1}) + c_{2}g_{1}x_{3} - c_{2}g_{2}x_{4} \\
\dot{x}_3 &= g_{3}(g_{4} - c_{3}u_1) \\
\dot{x}_4 &= g_{5}(g_{6} - c_{4}u_2) \\
u_1 &= \dot{l}_{m,1} \\
u_2 &= \dot{l}_{m,2} \\

"""

def dX1_dt(X):
	return(X[1])
def dX2_dt(X,U=None):
	if U==None:
		return(c1*np.sin(X[0]) + c2*G1(X)*X[2] + c2*G2(X)*X[3])
	else:
		return(c1*np.sin(X[0]) + c2*G1(X)*U[0] + c2*G2(X)*U[1])
def dX3_dt(X,U=None):
	if U == None:
		return(G3(X)*(G4(X) - c3*X[6]))
	else:
		return(G3(X)*(G4(X) - c3*U[0]))
def dX4_dt(X,U=None):
	if U == None:
		return(G5(X)*(G6(X) - c4*X[7]))
	else:
		return(G5(X)*(G6(X) - c4*U[1]))

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
		return(c1*np.sin(X[0]) + c2*G1(X)*U[0] + c2*G2(X)*U[1] - dA1(t,X))
	"""
	def A2(t,X):
		return(Z1(t,X) + dA1(t,X) - c1*np.sin(X[0]) - k2*Z2(t,X))
	Coefficient1 = c2*G1(X)
	Coefficient2 = c2*G2(X)
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
		return(c2*G1(X)*X[2] + c2*G2(X)*X[3] - A2(t,X))
	def dZ3(t,X):
		g1 = G1(X)
		g2 = G2(X)
		g3 = G3(X)
		g5 = G5(X)
		return(c2*dG1(X)*X[1]*X[2] + c2*dG2(X)*X[1]*X[3] \
						+ c2*g1*g3*G4(X) - c2*c3*g1*g3*U[0] \
							+ c2*g2*g5*G6(X) - c2*c4*g2*g5*U[1] \
								- dA2(t,X))
	def A3(t,X):
		return(Z2(t,X) - dA2(t,X) + k3*Z3(t,X) \
		+ c2*dG1(X)*X[1]*X[2] + 	c2*dG2(X)*X[1]*X[3] \
				+ c2*G1(X)*G3(X)*G4(X) + c2*G2(X)*G5(X)*G6(X))
	Coefficient1 = c2*c3*G1(X)*G3(X)
	Coefficient2 = c2*c4*G2(X)*G5(X)
	Constraint = A3(t,X)
	return(Coefficient1,Coefficient2,Constraint)

def return_U(t,X,U,dt,MaxStep,Bounds,Noise,Method):
	import random
	assert Method in ["Tension","Muscle Velocity"]
	if Method == "Tension":
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_tension_driven(t,X)
	elif Method == "Muscle Velocity":
			Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t,X)
	assert np.shape(Bounds)==(2,2), "Bounds must be (2,2)."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."
	assert Coefficient1!=0, "Error with 1st Coefficient. Shouldn't be zero."
	# assert Coefficient2!=0, "Error with 2nd Coefficient. Shouldn't be zero."

	if np.sign(-Coefficient1) == np.sign(Coefficient2):
		if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
			LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
		else:
			LowerBound = Bounds[0][0]

		if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
			UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
		else:
			UpperBound = Bounds[0][1]
	else: # np.sign(-Coefficient1) != np.sign(Coefficient2)
		if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
			LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
		else:
			LowerBound = Bounds[0][0]

		if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
			UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
		else:
			UpperBound = Bounds[0][1]
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
	if t<10*dt: MaxStep = 10*MaxStep
	FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
	FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
							for el in FeasibleInput1])
	euclid_dist = np.array(list(map(lambda u1,u2: np.sqrt((U[0]-u1)**2 + (U[1]-u2)**2),\
							FeasibleInput1,FeasibleInput2)))
	if Method == "Tension":
		feasible_index = np.where(np.logical_and(euclid_dist>=0, euclid_dist<=MaxStep))
	elif Method == "Muscle Velocity":
		feasible_index = np.where(np.logical_and(np.logical_and(euclid_dist>=0, euclid_dist<=MaxStep),np.sign(FeasibleInput1)!=np.sign(FeasibleInput2)))
	if len(feasible_index[0]) == 0: import ipdb; ipdb.set_trace()
	next_index = random.choice(feasible_index[0])
	u1 = FeasibleInput1[next_index]
	u2 = FeasibleInput2[next_index]
	return([u1,u2])


x1_1,x2_1 = [Base],[Amp*Freq]
u1_1,u2_1 = [100],[10]

x1_2,x2_2,x3_2,x4_2= [Base],[Amp*Freq],[170],[70]
u1_2,u2_2 = [-1],[1]


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
	U = return_U(t,X,U,dt,MaxStep_Tension,Tension_Bounds,NoiseArray[:,int(t/dt)],Method)
	u1_1.append(U[0])
	u2_1.append(U[1])
	x2_1.append(x2_1[-1] + dX2_dt(X,U=U)*dt)
	x1_1.append(x1_1[-1] + dX1_dt(X)*dt)
def update_policy_muscle_velocity_driven(t,x1_2,x2_2,x3_2,x4_2,dt,NoiseArray):
	import numpy as np
	Method = "Muscle Velocity"
	X = [x1_2[-1],x2_2[-1],x3_2[-1],x4_2[-1]]
	U = [u1_2[-1],u2_2[-1]]
	U = return_U(t,X,U,dt,MaxStep_Tension,MuscleVelocity_Bounds,NoiseArray[:,int(t/dt)],Method)
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
	statusbar(int(t/dt),len(Time),StartTime=StartTime,Title="T Driven Pendulum")

def plot_MA_values(Time,x1):
	import matplotlib.pyplot as plt
	import numpy as np

	fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8,6))
	plt.subplots_adjust(hspace=0.1,bottom=0.1)

	plt.suptitle("Moment arm equations")
	ax1.plot(np.linspace(0,np.pi*(160/180),1001),\
				np.array(list(map(lambda x1: G1([x1]),np.linspace(0,np.pi*(160/180),1001)))),\
				'0.70')
	ax1.plot(np.linspace(min(x1),max(x1),101),\
				np.array(list(map(lambda x1: G1([x1]),np.linspace(min(x1),max(x1),101)))),\
				'g',lw=3)
	ax1.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
	ax1.set_xticklabels([""]*len(ax1.get_xticks()))
	ax1.set_ylabel("Moment Arm for\n Muscle 1 (m)")

	ax2.plot(Time,np.array(list(map(lambda x1: G1([x1]),x1))),'g')
	ax2.set_ylim(ax1.get_ylim())
	ax2.set_yticks(ax1.get_yticks())
	ax2.set_yticklabels([""]*len(ax1.get_yticks()))
	ax2.set_xticklabels([""]*len(ax2.get_xticks()))

	ax3.plot(np.linspace(0,np.pi*(160/180),1001),\
				np.array(list(map(lambda x1: G2([x1]),np.linspace(0,np.pi*(160/180),1001)))),\
				'0.70')
	ax3.plot(np.linspace(min(x1),max(x1),101),\
				np.array(list(map(lambda x1: G2([x1]),np.linspace(min(x1),max(x1),101)))),\
				'r',lw=3)
	ax3.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
	ax3.set_xticklabels([r"$0$",r"$\frac{\pi}{4}$",r"$\frac{\pi}{2}$",r"$\frac{3\pi}{4}$",r"$\pi$"])
	ax3.set_xlabel("Joint Angle (rads)")
	ax3.set_ylabel("Moment Arm for\n Muscle 2 (m)")

	ax4.plot(Time,np.array(list(map(lambda x1: G2([x1]),x1))),'r')
	ax4.set_ylim(ax3.get_ylim())
	ax4.set_yticks(ax3.get_yticks())
	ax4.set_yticklabels([""]*len(ax3.get_yticks()))
	ax4.set_xlabel("Time (s)")
	return(fig,[ax1,ax2,ax3,ax4])

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
