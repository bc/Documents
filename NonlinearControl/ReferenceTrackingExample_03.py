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

t = sy.Symbol('t')

"""
x_1 &= \theta \\
x_2 &= \dot{\theta} \\
x_3 &= f_{T,1} \\
x_4 &= f_{T,2} \\
x_5 &= l_{m,1} \\
x_6 &= l_{m,2} \\
x_7 &= v_{m,1} \\
x_8 &= v_{m,2} \\
u_1 &= \alpha_1 \\
u_2 &= \alpha_2 \\

"""

g,l = 9.80, 0.45 #m/s², m
M = 2 # kg
α1 = 0 # 10*np.pi/180 # rads
α2 = 0 # 10*np.pi/180 # rads
m1 = 1 # kg
m2 = 1 # kg
bm1 = 1 # kg/s
bm2 = 1 # kg/s
AllMuscleSettings = return_muscle_settings()
[[r1,r2],[dr1,dr2],[d2r1,d2r2]]= return_MA_matrix_functions(AllMuscleSettings)
cT = 27.8
kT = 0.0047
LrT = 0.964
β = 1.55
ω = 0.75
ρ = 2.12
PCSA1 = 30**2*np.pi # mm²
PCSA2 = 30**2*np.pi # mm²
V_max = -9.15
cv0 = -5.78
cv1 = 9.18
av0 = -1.53
av1 = 0
av2 = 0
bv = 0.69
lo1 = AllMuscleSettings["BIC"]["Optimal Muscle Length"]/1000
lo2 = AllMuscleSettings["TRI"]["Optimal Muscle Length"]/1000
c_1 = 23.0
k_1 = 0.046
Lr1 = 1.17
η = 0.01
F_MAX1 = 0.25*PCSA1
F_MAX2 = 0.25*PCSA2

Amp = 7.5*np.pi/180
Base = 90*np.pi/180
Freq = 2*np.pi

k1,k2,k3,k4 = 100,10,10,100

max_allowable_act_jump = 0.3


"""
c_{1} &= -\frac{3g}{2l_{\text{arm}}} \\
c_{2} &= \frac{3}{m_{\text{arm}}l_{\text{arm}}^2} \\
c_{3} &= \frac{1}{\cos(\alpha_{1})} \\
c_{4} &= \frac{1}{\cos(\alpha_{2})} \\
c_{5} &= \frac{\cos(\alpha_{1})}{m_1} \\
c_{6} &= \frac{\cos^2(\alpha_{1})}{m_1} \\
c_{7} &= \frac{b_{m,1}\cos^2(\alpha_{1})}{m_1} \\
c_{8} &= \tan^2(\alpha_{1}) \\
c_{9} &= \frac{\cos(\alpha_{2})}{m_2} \\
c_{10} &= \frac{\cos^2(\alpha_{2})}{m_2} \\
c_{11} &= \frac{b_{m,2}\cos^2(\alpha_{2})}{m_2} \\
c_{12} &= \tan^2(\alpha_{2}) \\

"""

c1 = -(3*g)/(2*l)
c2 = 3/(M*l**2)
c3 = 1/(np.cos(α1))
c4 = 1/(np.cos(α2))
c5 = np.cos(α1)/m1
c6 = np.cos(α1)*c5
c7 = bm1*c6
c8 = np.tan(α1)**2
c9 = np.cos(α2)/m2
c10 = np.cos(α2)*c9
c11 = bm2*c10
c12 = np.tan(α2)**2

'''
g_{1} &= r_{1}\left(\theta,\theta_{PS}=\frac{\pi}{2}\right) \\
g_{2} &= r_{2}\left(\theta,\theta_{PS}=\frac{\pi}{2}\right)  \\
g_{3} &= \kappa^{T}(f_{T,1}) \hspace{1em} \text{(Variable stiffness coefficient from Zajac (1989) ODE for tendon force)} \\
g_{4} &= \kappa^{T}(f_{T,2}) \hspace{1em} \\
g_{5} &= \text{sgn}(-r_1(\theta)\cdot\dot{\theta})\cdot\dot{\theta}\cdot\sqrt{\left(\frac{\partial r_1}{\partial\theta}\right)^2 + r_1^2(\theta)} \\
g_{6} &= \text{sgn}(-r_2(\theta)\cdot\dot{\theta})\cdot\dot{\theta}\cdot\sqrt{\left(\frac{\partial r_2}{\partial\theta}\right)^2 + r_2^2(\theta)} \\
g_{7} &= F_{\text{MAX},1}\cdot f_{L}(l_{m,1}) \cdot f_{V}(l_{m,1},v_{m,1}) \\
g_{8} &= F_{PE}(l_{m,1},v_{m,1}) \\
g_{9} &= F_{PE}(l_{m,2},v_{m,2}) \\
g_{10} &= F_{\text{MAX},2}\cdot f_{L}(l_{m,2}) \cdot f_{V}(l_{m,2},v_{m,2}) \\
'''
# FL = lambda l,lo: 1
# FV = lambda l,v,lo: 1
FL = lambda l,lo: np.exp(-abs(((l/lo)**β-1)/ω)**ρ)
FV = lambda l,v,lo: np.piecewise(v,[v<=0, v>0],\
	[lambda v: (V_max - v/lo)/(V_max + (cv0 + cv1*(l/lo))*(v/lo)),\
	lambda v: (bv-(av0 + av1*(l/lo) + av2*(l/lo)**2)*(v/lo))/(bv + (v/lo))])

def G1(X):
	return(lambdify([θ_EFE],r1.subs([(θ_PS,np.pi/2)]))(X[0])) #
def dG1(X):
	return(lambdify([θ_EFE],dr1.subs([(θ_PS,np.pi/2)]))(X[0]))
def d2G1(X):
	return(lambdify([θ_EFE],d2r1.subs([(θ_PS,np.pi/2)]))(X[0]))
def G2(X):
	return(lambdify([θ_EFE],r2.subs([(θ_PS,np.pi/2)]))(X[0])) #
def dG2(X):
	return(lambdify([θ_EFE],dr2.subs([(θ_PS,np.pi/2)]))(X[0]))
def d2G2(X):
	return(lambdify([θ_EFE],d2r2.subs([(θ_PS,np.pi/2)]))(X[0]))
def G3(X):
	return((cT*F_MAX1/lo1)*(1-np.exp(-X[2]/(cT*kT)))) # Maybe its negative because of theta sign convention
def dG3(X):
	return(np.exp(-X[2]/(cT*kT))/(kT*lo1))
def G4(X):
	return(np.sign(-G1(X))*X[1]*np.sqrt(dG1(X)**2 + G1(X)**2)/lo1) # Should be G5
def dG4_dt(X):
	return(np.sign(-G1(X))*(dX2_dt(X)*np.sqrt(dG1(X)**2 + G1(X)**2) +\
		  ((X[1]**2)*dG1(X)*(d2G1(X) + G1(X)))/np.sqrt(dG1(X)**2 + G1(X)**2))/lo1)
def G5(X):
	return((cT*F_MAX2/lo2)*(1-np.exp(-X[3]/(cT*kT))))
def dG5(X):
	return(np.exp(-X[3]/(cT*kT))/(kT*lo2))
def G6(X):
	return(np.sign(-G2(X))*X[1]*np.sqrt(dG2(X)**2 + G2(X)**2)/lo2) # Should be G4
def dG6_dt(X):
	return(np.sign(-G2(X))*(dX2_dt(X)*np.sqrt(dG2(X)**2+ G2(X)**2) +\
			(X[1]**2*dG2(X)*(d2G2(X) + G2(X)))/np.sqrt(dG2(X)**2 + G2(X)**2))/lo2)
def G7(X):
	return(F_MAX1*(c_1*k_1*np.log(np.exp(((X[4]/(lo1*1.2))-Lr1)/k_1) + 1) + η*(X[6]/lo1)))  # Should be G8
def G8(X):
	return(0.25*PCSA1*FL(X[4],lo1)*FV(X[4],X[6],lo1)) # Should be G6
def G9(X):
	return(F_MAX2*(c_1*k_1*np.log(np.exp(((X[5]/(lo2*1.2))-Lr1)/k_1) + 1) + η*(X[7]/lo2)))
def G10(X):
	return(0.25*PCSA2*FL(X[5],lo2)*FV(X[5],X[7],lo2)) # Should be G7

"""
\dot{x}_1 &= x_{2} \\
\dot{x}_2 &= c_{1}\cos(x_{1}) + c_{2}g_{1}x_{3} - c_{2}g_{2}x_{4} \\
\dot{x}_3 &= g_{3}g_{5} - c_{3}g_{3}x_{7} \\
\dot{x}_4 &= g_{4}g_{6} - c_{4}g_{4}x_{8} \\
\dot{x}_5 &= x_{7} \\
\dot{x}_6 &= x_{8} \\
\dot{x}_7 &= c_{5}x_{3} - c_{6}g_{8} - c_{7}x_{7} + \frac{c_{8}x_{7}^2}{x_{5}} - c_{6}g_{7}u_{1}\\
\dot{x}_8 &= c_{9}x_{4} - c_{10}g_{9} - c_{11}x_{8} + \frac{c_{12}x_{8}^2}{x_{6}} - c_{10}g_{10}u_{2} \\
"""

def dX1_dt(X):
	return(X[1])
def dX2_dt(X):
	return(c1*np.sin(X[0]) + c2*G1(X)*X[2] + c2*G2(X)*X[3])
def d2X2_dt2(X):
	return(c1*np.cos(X[0])*dX1_dt(X) + c2*dG1(X)*dX1_dt(X)*X[2] + c2*G1(X)*dX3_dt(X)\
					+ c2*dG2(X)*dX1_dt(X)*X[3] + c2*G2(X)*dX4_dt(X))
def dX3_dt(X):
	return(G3(X)*G4(X) - c3*G3(X)*X[6])
def dX4_dt(X):
	return(G5(X)*G6(X) - c4*G5(X)*X[7])
def dX5_dt(X):
	return(X[6])
def dX6_dt(X):
	return(X[7])
def dX7_dt(X,U):
	return(c5*X[2] - c6*G8(X)*U[0] - c6*G7(X) - c7*X[6] \
	+ c8*X[6]**2/X[4])
def dX8_dt(X,U):
	return(c9*X[3] - c10*G10(X)*U[1] - c10*G9(X) - c11*X[7] \
	+ c12*X[7]**2/X[5])

r = lambda t: Amp*np.sin(Freq*t) + Base
dr = lambda t: Amp*Freq*np.cos(Freq*t)
d2r = lambda t: -Amp*Freq**2*np.sin(Freq*t)
d3r = lambda t: -Amp*Freq**3*np.cos(Freq*t)
d4r = lambda t: Amp*Freq**4*np.sin(Freq*t)

def Z1(t,X):
	return(r(t) - X[0])
def dZ1(t,X):
	return(dr(t) - X[1])
def d2Z1(t,X):
	return(d2r(t) - dX2_dt(X))
def d3Z1(t,X):
	return(d3r(t) - d2X2_dt2(X))
def A1(t,X):
	return(dr(t) + k1*Z1(t,X))
def dA1(t,X):
	return(d2r(t) + k1*dZ1(t,X))
def d2A1(t,X):
	return(d3r(t) + k1*d2Z1(t,X))
def d3A1(t,X):
	return(d4r(t) + k1*d3Z1(t,X))
def Z2(t,X):
	return(X[1] - A1(t,X))
def dZ2(t,X):
	return(c1*np.sin(X[0]) + c2*G1(X)*X[2] + c2*G2(X)*X[3] - dA1(t,X))
def d2Z2(t,X):
	dx1_dt = dX1_dt(X)
	return(c1*np.cos(X[0])*dx1_dt + c2*dG1(X)*dx1_dt*X[2] + c2*G1(X)*dX3_dt(X) \
					+ c2*dG2(X)*dx1_dt*X[3] + c2*G2(X)*dX4_dt(X) - d2A1(t,X))
def A2(t,X):
	return(Z1(t,X) + dA1(t,X) - c1*np.sin(X[0]) - k2*Z2(t,X))
def dA2(t,X):
	return(dZ1(t,X) + d2A1(t,X) - c1*np.cos(X[0])*X[1] - k2*dZ2(t,X))
def d2A2(t,X):
	return(d2Z1(t,X) + d3A1(t,X) + c1*np.sin(X[0])*X[1]**2 - k2*d2Z2(t,X))
def Z3(t,X):
	return(c2*G1(X)*X[2] + c2*G2(X)*X[3] - A2(t,X))
def dZ3(t,X):
	g1 = G1(X)
	g2 = G2(X)
	g3 = G3(X)
	g5 = G5(X)
	return(c2*dG1(X)*X[1]*X[2] + c2*dG2(X)*X[1]*X[3] \
					+ c2*g1*g3*G4(X) - c2*c3*g1*g3*X[6] \
						+ c2*g2*g5*G6(X) - c2*c4*g2*g5*X[7] \
							- dA2(t,X))
def A3(t,X):
	return(Z2(t,X) + c2*dG1(X)*X[1]*X[2] + c2*dG2(X)*X[1]*X[3] \
					+ c2*G1(X)*G3(X)*G4(X) + c2*G2(X)*G5(X)*G6(X) \
						- dA2(t,X) + k3*Z3(t,X))
def dA3(t,X):
	g1 = G1(X)
	dg1 = dG1(X)
	g2 = G2(X)
	dg2 = dG2(X)
	dx1_dt = dX1_dt(X)
	dx2_dt = dX2_dt(X)
	dx3_dt = dX3_dt(X)
	dx4_dt = dX4_dt(X)
	g3 = G3(X)
	g4 = G4(X)
	g5 = G5(X)
	g6 = G6(X)
	return(dZ2(t,X) + c2*d2G1(X)*(X[1]**2)*X[2] + c2*dg1*dx2_dt*X[2]\
	 		+ c2*dg1*X[1]*dx3_dt + c2*d2G2(X)*(X[1]**2)*X[3] \
				+ c2*dg2*dx2_dt*X[3] + c2*dg2*X[1]*dx4_dt \
					+ c2*dg1*dx1_dt*g3*g4 + c2*g1*dG3(X)*dx3_dt*g4 \
						+ c2*g1*g3*dG4_dt(X) \
							+ c2*dg2*dx1_dt*g5*g6 + c2*g2*dG5(X)*dx4_dt*g6 \
								+ c2*g2*g5*dG6_dt(X) \
									- d2A2(t,X) + k3*dZ3(t,X))
def Z4(t,X):
	return(c2*c3*G1(X)*G3(X)*X[6] + c2*c4*G2(X)*G5(X)*X[7] - A3(t,X))
def dZ4(t,X,U):
	g1 = G1(X)
	g2 = G2(X)
	g3 = G3(X)
	dg3 = dG3(X)
	g5 = G5(X)
	dg5 = dG5(X)
	return(c2*c3*X[6]*(dG1(X)*g3*X[1] + g1*dg3*g3*G4(X) - c7*g1*g3) \
			+ c2*c4*X[7]*(dG2(X)*g5*X[1] + g2*dg5*g5*G6(X) - c11*g2*g5) \
			+ c2*c3*X[6]**2*(c8*g1*g3/X[4] - c3*g1*dg3*g3) \
			+ c2*c4*X[7]**2*(c12*g2*g5/X[5] - c4*g2*dg5*g5) \
			+ c2*c3*c5*g1*g3*X[2] - c2*c3*c6*g1*g3*G7(X) \
			+ c2*c4*c9*g2*g5*X[3] - c2*c4*c10*g2*g5*G9(X) \
			- c2*c3*c6*g1*g3*G8(X)*U[0] - c2*c4*c10*g2*g5*G10(X)*U[1] \
			-dA3(t,X))
Coefficient_11 = lambda t,X: -c2*c3*c6*G1(X)*G3(X)*G8(X)
Coefficient_12 = lambda t,X: c2*c4*c10*G2(X)*G5(X)*G10(X)
def Constraint_1(t,X):
	g1 = G1(X)
	g2 = G2(X)
	g3 = G3(X)
	g5 = G5(X)
	dx1_dt = dX1_dt(X)
	return(	c2*c3*dG1(X)*dx1_dt*g3*X[6] + c2*c3*g1*dG3(X)*dX3_dt(X)*X[6] \
			+ c2*c4*dG2(X)*dx1_dt*g5*X[7] + c2*c4*g2*dG5(X)*dX4_dt(X)*X[7] \
			+ c2*c3*c5*g1*g3*X[2] - c2*c3*c6*g1*g3*G7(X) - c2*c3*c7*g1*g3*X[6] \
			+ c2*c4*c9*g2*g5*X[3] - c2*c4*c10*g2*g5*G9(X) - c2*c4*c11*g2*g5*X[7] \
			+ c2*c3*c8*g1*g3*((X[6]**2)/X[4]) + c2*c4*c12*g2*g5*((X[7]**2)/X[5]) \
			- dA3(t,X) - Z3(t,X) + k4*Z4(t,X))
					# = c2*c3*c6*G1(X)*G3(X)*G8(X)*U[0] + c2*c4*c10*G2(X)*G5(X)*G10(X)*U[1]
# Constraint_1 = lambda t,X: \
# 				c2*c3*X[6]*(dG1(X)*G3(X)*X[1] + G1(X)*dG3(X)*G3(X)*G4(X) - c7*G1(X)*G3(X)) \
# 				+ c2*c4*X[7]*(dG2(X)*G5(X)*X[1] + G2(X)*dG5(X)*G5(X)*G6(X) - c11*G2(X)*G5(X)) \
# 				+ c2*c3*X[6]**2*(c8*G1(X)*G3(X)/X[4] - c3*G1(X)*dG3(X)*G3(X)) \
# 				+ c2*c4*X[7]**2*(c12*G2(X)*G5(X)/X[5] - c4*G2(X)*dG5(X)*G5(X)) \
# 				+ c2*c3*c5*G1(X)*G3(X)*X[2] - c2*c3*c6*G1(X)*G3(X)*G7(X) \
# 				+ c2*c4*c9*G2(X)*G5(X)*X[3] - c2*c4*c10*G2(X)*G5(X)*G9(X) \
# 				-dA3(t,X) - Z3(t,X) + k4*Z4(t,X) \
# 					# = c2*c3*c6*G1(X)*G3(X)*G8(X)*U[0] + c2*c4*c10*G2(X)*G5(X)*G10(X)*U[1]
# Z1 = lambda t,X: r(t) - X[0]
# dZ1 = lambda t,X: dr(t) - X[1]
# d2Z1 = lambda t,X: d2r(t) - dX2_dt(X)
# d3Z1 = lambda t,X: d3r(t) - d2X2_dt2(X)
# a1 = lambda t,X: dr(t) - k1*Z1(t,X)
# da1 = lambda t,X: d2r(t) - k1*dZ1(t,X)
# d2a1 = lambda t,X: d3r(t) - k1*d2Z1(t,X)
# d3a1 = lambda t,X: d4r(t) - k1*d3Z1(t,X)
# z2 = lambda t,X: X[1] - A1(t,X)
# dz2 = lambda t,X: c1*np.cos(X[0]) + c2*G1(X)*X[2] - c2*G2(X)*X[3] - dA1(t,X)
# d2z2 = lambda t,X: -c1*np.sin(X[0])*dX1_dt(X) + c2*dG1(X)*dX1_dt(X)*X[2] + c2*G1(X)*dX3_dt(X) \
# 					- c2*dG2(X)*dX1_dt(X)*X[3] - c2*G2(X)*dX4_dt(X) - d2A1(t,X)
# a2 = lambda t,X: Z1(t,X) + dA1(t,X) - c1*np.cos(X[0]) - k2*Z2(t,X)
# da2 = lambda t,X: dZ1(t,X) + d2A1(t,X) + c1*np.sin(X[0])*X[1] - k2*dZ2(t,X)
# d2a2 = lambda t,X: d2Z1(t,X) + d3A1(t,X) + c1*np.cos(X[0])*X[1]**2 - k2*d2Z2(t,X)
# z3 = lambda t,X: c2*G1(X)*X[2] - c2*G2(X)*X[3] - A2(t,X)
# dz3 = lambda t,X: c2*dG1(X)*X[1]*X[2] - c2*dG2(X)*X[1]*X[3] \
# 					+ c2*G1(X)*G3(X)*G4(X) - c2*c3*G1(X)*G3(X)*X[6] \
# 						-c2*G2(X)*G5(X)*G6(X) + c2*c4*G2(X)*G5(X)*X[7] \
# 							- dA2(t,X)
# a3 = lambda t,X: Z2(t,X) + c2*dG1(X)*X[1]*X[2] - c2*dG2(X)*X[1]*X[3] \
# 					+ c2*G1(X)*G3(X)*G4(X) - c2*G2(X)*G5(X)*G6(X) \
# 						- dA2(t,X) + k3*Z3(t,X)
# da3 = lambda t,X: dZ2(t,X) + c2*d2G1(X)*X[1]**2*X[2] + c2*dG1(X)*dX2_dt(X)*X[2]\
# 			 		+ c2*dG1(X)*X[1]*dX3_dt(X) - c2*d2G2(X)*X[1]**2*X[3] \
# 						- c2*dG2(X)*dX2_dt(X)*X[3] - c2*dG2(X)*X[1]*dX4_dt(X) \
# 							+ c2*dG1(X)*dX1_dt(X)*G3(X)*G4(X) + c2*G1(X)*dG3(X)*dX3_dt(X)*G4(X) \
# 								+ c2*G1(X)*G3(X)*dG4_dt(X) \
# 									- c2*dG2(X)*dX1_dt(X)*G5(X)*G6(X) - c2*G2(X)*dG5(X)*dX4_dt(X)*G6(X) \
# 										- c2*G2(X)*G5(X)*dG6_dt(X) \
# 											- d2A2(t,X) + k3*dZ3(t,X)
# z4 = lambda t,X: c2*c3*G1(X)*G3(X)*X[6] - c2*c4*G2(X)*G5(X)*X[7] - A3(t,X)
# dz4 = lambda X,U,t: \
# 				c2*c3*X[6]*(dG1(X)*G3(X)*X[1] + G1(X)*dG3(X)*G3(X)*G4(X) - c7*G1(X)*G3(X)) \
# 				- c2*c4*X[7]*(dG2(X)*G5(X)*X[1] + G2(X)*dG5(X)*G5(X)*G6(X) - c11*G2(X)*G5(X)) \
# 				+ c2*c3*X[6]**2*(c8*G1(X)*G3(X)/X[4] - c3*G1(X)*dG3(X)*G3(X)) \
# 				- c2*c4*X[7]**2*(c12*G2(X)*G5(X)/X[5] - c4*G2(X)*dG5(X)*G5(X)) \
# 				+ c2*c3*c5*G1(X)*G3(X)*X[2] - c2*c3*c6*G1(X)*G3(X)*G7(X) \
# 				- c2*c4*c9*G2(X)*G5(X)*X[3] + c2*c4*c10*G2(X)*G5(X)*G9(X) \
# 				- c2*c3*c6*G1(X)*G3(X)*G8(X)*U[0] + c2*c4*c10*G2(X)*G5(X)*G10(X)*U[1] \
# 				-dA3(t,X)
# Coefficient_11 = lambda t,X: c2*c3*c6*G1(X)*G3(X)*G8(X)
# Coefficient_12 = lambda t,X: -c2*c4*c10*G2(X)*G5(X)*G10(X)
# Constraint_1 = lambda t,X: \
# 				c2*c3*X[6]*(dG1(X)*G3(X)*X[1] + G1(X)*dG3(X)*G3(X)*G4(X) - c7*G1(X)*G3(X)) \
# 				- c2*c4*X[7]*(dG2(X)*G5(X)*X[1] + G2(X)*dG5(X)*G5(X)*G6(X) - c11*G2(X)*G5(X)) \
# 				+ c2*c3*X[6]**2*(c8*G1(X)*G3(X)/X[4] - c3*G1(X)*dG3(X)*G3(X)) \
# 				- c2*c4*X[7]**2*(c12*G2(X)*G5(X)/X[5] - c4*G2(X)*dG5(X)*G5(X)) \
# 				+ c2*c3*c5*G1(X)*G3(X)*X[2] - c2*c3*c6*G1(X)*G3(X)*G7(X) \
# 				- c2*c4*c9*G2(X)*G5(X)*X[3] + c2*c4*c10*G2(X)*G5(X)*G9(X) \
# 				-dA3(t,X) - Z3(t,X) + k4*Z4(t,X) \
# 					# = c2*c3*c6*G1(X)*G3(X)*G8(X)*U[0] - c2*c4*c10*G2(X)*G5(X)*G10(X)*U[1]

def return_U(t,X,U,e,Noise):
	import random
	Constraint1 = Constraint_1(t,X)
	d1 =Coefficient_11(t,X)
	d2 = Coefficient_12(t,X)
	# assert d1>=0, "Error with 1st Coefficient. Should be positive."
	# assert d2<=0, "Error with 2nd Coefficient. Should be negative."
	# assert Constraint1>d2, "HUGE ERROR... Constraints create infeasible activations (i.e., negative)"
	Feasible = (True,"YAY!")
	if Constraint1 >=0:
		# assert d1>0 or d2>0, "HUGE PROBLEM!! Infeasible activations!"
		if d1 > 0 and d2 > 0:
			lowerbound = 0
			if Constraint1<d1:
				upperbound = Constraint1/d1
			else:
				upperbound = 1
		elif d1 > 0 and d2 < 0:
			# assert Constraint1 < d1, "HUGE ERROR, infeasible activations"
			if Constraint1<d1:
				lowerbound = Constraint1/d1
				if Constraint1 < d1 + d2:
					upperbound = (Constraint1 - d2)/d1
				else:
					upperbound = 1
			else:
				Feasible = (False,1)
		elif d1 < 0 and d2 > 0:
			# assert Constraint1 < d2, "HUGE ERROR, infeasible activations"
			if Constraint1 < d2:
				lowerbound = 0
				if Constraint1 > d1 + d2:
					upperbound = (Constraint1 - d2)/d1
				else:
					upperbound = 1
			else:
				Feasible = (False,2)
		else:
			Feasible = (False,3)
	else: #Constraint1 <0
		# assert d1<0 or d2<0, "HUGE PROBLEM!! Infeasible activations!"
		if d1 < 0 and d2 < 0:
			lowerbound = 0
			if Constraint1>d1:
				upperbound = Constraint1/d1
			else:
				upperbound = 1
		elif d1 < 0 and d2 > 0:
			# assert Constraint1 > d1, "HUGE ERROR, infeasible activations"
			if Constraint1 > d1:
				lowerbound = Constraint1/d1
				if Constraint1 > d1 + d2:
					upperbound = (Constraint1 - d2)/d1
				else:
					upperbound = 1
			else:
				Feasible = (False,4)
		elif d1 > 0 and d2 < 0:
			# assert Constraint1 > d2, "HUGE ERROR, infeasible activations"
			if Constraint1 > d2:
				lowerbound = 0
				if Constraint1 < d1 + d2:
					upperbound = (Constraint1 - d2)/d1
				else:
					upperbound = 1
			else:
				Feasible = (False,5)
		else:
			Feasible = (False,6)
	if Feasible[0] == True:
		poss_u1 = (upperbound-lowerbound)*np.random.rand(1000) + lowerbound
		poss_u2 = np.array([Constraint1/d2 - (d1/d2)*el for el in poss_u1])
		euclid_dist = np.array(list(map(lambda u1,u2: np.sqrt((U[0]-u1)**2 + (U[1]-u2)**2),\
								poss_u1,poss_u2)))
		feasible_index = np.where(np.logical_and(euclid_dist>=0, euclid_dist<=max_allowable_act_jump))
		if len(feasible_index[0]) == 0: import ipdb; ipdb.set_trace()
		next_index = random.choice(feasible_index[0])
		u1 = poss_u1[next_index]
		u2 = poss_u2[next_index]
	else:
		u1 = U[0]
		u2 = U[1]
	# poss_u1 = poss_u1[feasible_index]
	# poss_u2 = poss_u2[feasible_index]
	# assert len(poss_u1)!=0, "Error! No feasible activations. Perhaps a finer sample is needed."
	# import ipdb; ipdb.set_trace()
	# if Constraint1<=0:
	#     u1 = Constraint1/(Coefficient_11(t,X) + e*Coefficient_12(t,X)) + Noise[0]
	#     u2 = e*(u1-Noise[0]) + Noise[1]
	# else:
	#     u2 = Constraint1/(e*Coef,ficient_11(t,X) + Coefficient_12(t,X)) + Noise[1]
	#     u1 = e*(u2-Noise[1]) + Noise[0]
	return([u1,u2],Feasible)

N = 10001
Time = np.linspace(0,10,N)
dt = Time[1]-Time[0]
x1,x2,x3,x4,x5,x6,x7,x8 = [Base],[Amp*Freq],[2000],[2000],[lo1*0.9],[lo2*1.1],[-1],[1]
u1,u2 = [0.1],[0.1]
CocontractionIndex = 2

AddNoise = False
if AddNoise == True:
    np.random.seed(seed=None)
    NoiseArray = np.random.normal(loc=0.0,scale=0.2,size=(2,len(Time)))
else:
    NoiseArray = np.zeros((2,len(Time)))

test = []
def update_policy(t,x1,x2,x3,x4,x5,x6,x7,x8,dt,NoiseArray,e=2):
	import numpy as np
	X = [x1[-1],x2[-1],x3[-1],x4[-1],x5[-1],x6[-1],x7[-1],x8[-1]]
	U = [u1[-1],u2[-1]]
	U,Feasible = return_U(t,X,U,e,NoiseArray[:,int(t/dt)])
	test.append(Feasible)
	u1.append(U[0])
	u2.append(U[1])
	x8.append(x8[-1] + dX8_dt(X,U)*dt)
	x7.append(x7[-1] + dX7_dt(X,U)*dt)
	x6.append(x6[-1] + dX6_dt(X)*dt)
	x5.append(x5[-1] + dX5_dt(X)*dt)
	x4.append(x4[-1] + dX4_dt(X)*dt)
	x3.append(x3[-1] + dX3_dt(X)*dt)
	x2.append(x2[-1] + dX2_dt(X)*dt)
	x1.append(x1[-1] + dX1_dt(X)*dt)


StartTime = time.time()
for t in Time[1:]:
	update_policy(t,x1,x2,x3,x4,x5,x6,x7,x8,dt,NoiseArray,e=CocontractionIndex)
	statusbar(int(t/dt),len(Time[1:]),StartTime=StartTime,Title="mm Driven Pendulum")

plt.figure()
# plt.title(r'$\dot{x}_{1} = x_{1}^{2} + x_{2}; \hspace{1em} \dot{x}_{2} = u$',\
#                 fontsize=16,color='gray')
plt.title("Underdetermined Muscle-Driven\nForced Pendulum Example",\
                fontsize=16,color='gray')
plt.plot(Time,x1,'b',lw=2)
plt.plot(Time,r(Time),'r--')
plt.xlabel("Time (s)")
plt.ylabel("Desired Measure")
plt.legend([r"Output $y = x_{1}$",r"Reference $r(t) = " + str(Amp) + "\sin(" + str(Freq) + "t)$ " + str(Base)])

plt.figure()
plt.title('Error vs. Time')
plt.plot(Time, r(Time)-x1,color='r')
plt.xlabel("Time (s)")
plt.ylabel("Error")

plt.figure()
plt.plot(Time[-1:],u1,'g',Time[-1:],u2,'r')
plt.title('Muscle Activations vs. Time')
plt.xlabel("Time (s)")
plt.ylabel("Muscle Activations")
plt.legend(["Muscle 1","Muscle 2"])

# Constraint_1 = Constraint_1.subs([(sy.Derivative(r,t),dr),(sy.Derivative(dr,t),d2r)\
# 									,(sy.Derivative(d2r,t),d3r),(sy.Derivative(d3r,t),d4r)])
#
# ConstraintFunction_1 = lambdify([x1,x2,x3,x4,x5,x6,x7,x8,r,dr,d2r,d3r,d4r],Constraint_1)
# Coefficient11Function = lambdify([x1,x2,x3,x4,x5,x6,x7,x8],Coefficient_11)
# Coefficient12Function = lambdify([x1,x2,x3,x4,x5,x6,x7,x8],Coefficient_12)
# """
# import numpy as np
# from scipy.integrate import odeint
# import matplotlib.pyplot as plt
# import sympy as sy
# from sympy.utilities import lambdify
# import time
#
# def return_muscle_settings():
# 	"""
# 	# Notes:
# 	# Coefficients from observation, Ramsay, FVC, Holtzbaur, Pigeon, Kuechle, or Banks. Optimal Muscle Length given in mm.
# 	#
# 	# BIC EFE MA for Ramsay has R² = 0.985 whereas Pigeon has R² = 0.9918. Pigeon, however, only takes elbow angle into account, whereas Ramsay takes in variable PS angles. It appears that because Pigeon uses an average of fully pronated and fully supinated MAs, the BIC moment arm is similar but subject to variation as the level of PS is changed. Coefficients and equation number/type are listed below to test either implementation. (NOTE: BIC becomes slightly negative when x1 > 3.021. If trajectory has elbow angles exceding this value, enter a threshold of 3.021 into the model.)
# 	#
# 	# Additionally, the SFE MA for the BIC is held constant in Pigeon at 29.21 mm while it was estimated as 15 mm.
# 	#
# 	# src = 'Ramsay', eq = 2, Coefficients = [8.4533,36.6147,2.4777,-19.432,2.0571,0,13.6502,0,0,-5.6172,0,-2.0854,0,0,0,0], threshold = 3.021
# 	# src = 'Pigeon', dof = 'elbow', eq = None, Coefficients = Pigeon_coeff_conversion([14.660,4.5322,1.8047,-2.9883,0,0]), threshold = 2.9326
# 	#
# 	# TRI EFE MA for Ramsay has R² = 0.997 whereas Pigeon has R² = 0.9904. Pigeon appears to really fail when the elbow angle is greater than 140°. For this reason, Ramsay should be used. However the approach of fixing the MA for values greater than 140° can be adopted for completeness. Coefficients and equation number/type are listed below to test either implementation.
# 	#
# 	# Additionally, the SFE MA for the TRI is held constant in Pigeon at -25.40 mm while it was estimated as -15 mm.
# 	#
# 	# src = 'Ramsay', eq = 1, Coefficients = [-24.5454,-8.8691,9.3509,-1.7518,0]
# 	# src = 'Pigeon', dof = 'elbow', eq = None, Coefficients = Pigeon_coeff_conversion([-23.287,-3.0284,12.886,-19.092,13.277,-3.5171])
#
# 	"""
# 	import sympy as sp
# 	from sympy.utilities import lambdify
# 	import numpy as np
# 	from numpy import pi
#
# 	global θ_SFE,θ_EFE,θ_PS
# 	θ_SFE,θ_EFE,θ_PS = sp.symbols('θ_SFE'),sp.symbols('θ_EFE'),sp.symbols('θ_PS')
#
# 	# Coefficients from observation, Ramsay, Pigeon, FVC, Holtzbaur, or Banks.
# 	# Moment arms are in mm. Mass is in grams. threshold is in radians.
#
# 	def Pigeon_coeff_conversion(Coefficients):
# 		"""
# 		# Takes in Coefficient values from Pigeon (1996) -- which take in angles in degrees -- and coverts them into the properly scaled coefficients for radians, additionally scaled by the magnitude listed in the paper.
# 		#
# 		# Note that the coefficients listed in Pigeon (1996) are given in decending order (i.e., c5,c4,c3,c₂,c₁,c0). However to maintain continuity with the equations given in Ramsay (2009), we list coefficients in order of increasing power (i.e., c0,c1,c2,c3,c4,c5).
# 		"""
# 		import numpy as np
# 		assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
# 		assert type(Coefficients)==list, 'Coefficients must be a 6 element list.'
# 		Rad_Conversion = np.multiply(Coefficients,\
# 				np.array([1,(180/np.pi),(180/np.pi)**2,(180/np.pi)**3,(180/np.pi)**4,(180/np.pi)**5],dtype = 'float64'))
# 		new_Coefficients =\
# 			np.multiply(Rad_Conversion,np.array([1,1e-1,1e-3,1e-5,1e-7,1e-9],dtype='float64'))
# 		return(new_Coefficients)
#
# 	BIC_Settings = {\
# 		'Shoulder' : {\
# 			'MA Coefficients' : [29.21,0,0,0,0,0],\
# 			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : None, \
# 			'dof' : 'Shoulder'}, \
# 		'Elbow' : {\
# 			'MA Coefficients' : [8.4533,36.6147,2.4777,-19.432,2.0571,0,13.6502,0,0,-5.6172,0,-2.0854,0,0,0,0],\
# 			'Source' : 'Ramsay', 'Equation Number' : 2, 'Threshold' : 3.021, \
# 			'dof' : 'Elbow'}, \
# 		'Mass' : 163.8,\
# 		'Actual No' : 320,\
# 		'Corrected No' : 292.6,\
# 		'Relative Abundance' : 1.1,\
# 		'Optimal Muscle Length' : 116,\
# 		'Group' : 'flexor'}
# 	TRI_Settings = {\
# 		'Shoulder' : {\
# 			'MA Coefficients' : [-25.40,0,0,0,0,0], \
# 			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : None, \
# 			'dof' : 'Shoulder'}, \
# 		'Elbow' : {\
# 			'MA Coefficients' : [-24.5454,-8.8691,9.3509,-1.7518,0],\
# 			'Source' : 'Ramsay', 'Equation Number' : 1, 'Threshold' : None, \
# 			'dof' : 'Elbow'}, \
# 		'Mass' : (94.2+138.4+92.5), \
# 		'Actual No' : (200+222+98),\
# 		'Corrected No' : (223.7+269.6+221.8),\
# 		'Relative Abundance' : (0.89+0.82+0.44)/3,\
# 		'Optimal Muscle Length' : 134,\
# 		'Group' : 'extensor'}
#
# 	AllMuscleSettings = {'BIC' : BIC_Settings, 'TRI' : TRI_Settings}
# 	return(AllMuscleSettings)
# def return_MA_matrix_functions(AllMuscleSettings):
# 	import numpy as np
# 	import sympy as sp
# 	from sympy.utilities import lambdify
# 	def MA_function(Parameters):
# 		"""
# 		# Note:
# 		#
# 		# Angles should be a number if Coefficients has a length of 5, or a list of length 2 when the Coefficients have lengths 16 or 18. Angles[0] will be the PRIMARY ANGLE for the DOF being considered while Angles[1] will be the secondary angle.
# 		#
# 		# Notes:
# 		#
# 		# threshold is only needed for Pigeon or Ramsay MA functions that are invalid outside of a given value. Must be either None (default) or the radian value of the threshold.
# 		#
# 		# dof only needed for Pigeon (Ramsay only handles EFE for this 2 DOF system). Must be either 'Shoulder' or 'Elbow'.
# 		#
# 		# eq is only needed for Ramsay (Pigeon has one quintic polynomial). eq must be either 1, 2, or 3, with list length requirements of 5, 16, or 18, respectively.
# 		"""
# 		import sympy as sp
# 		import numpy as np
#
# 		src = Parameters['Source']
# 		Coefficients = Parameters['MA Coefficients']
# 		eq = Parameters['Equation Number']
# 		dof = Parameters['dof']
# 		threshold = Parameters['Threshold']
#
# 		global θ_SFE,θ_EFE,θ_PS
# 		assert type(src) == str, "src must be a str."
# 		assert src.capitalize() in ['Ramsay','Pigeon','Est'], "src must be either Ramsay, Pigeon or Est (Estimate)."
# 		if dof != None:
# 			assert type(dof) == str, "dof must be a str."
# 			assert dof.capitalize() in ['Shoulder','Elbow'], "dof must be either Shoulder or Elbow."
# 		if src.capitalize() == 'Pigeon' :
# 			assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
# 			assert dof != None, "For Pigeon (1996), dof must be stated."
# 			eq = None
# 			if dof.capitalize() == 'Elbow' :
# 				θ = θ_EFE
# 			else:
# 				θ = θ_SFE
# 			MomentArm = (np.matrix(Coefficients,dtype='float64')\
# 							*np.matrix([1,θ,θ**2,θ**3,θ**4,θ**5]).T)[0,0]/1000
# 		elif src.capitalize() == 'Est' :
# 			MomentArm = np.array(Coefficients,dtype='float64')/1000
# 		else: #src.capitalize() == 'Ramsay'
# 			θ = θ_EFE
# 			assert type(Coefficients) == list, "Coefficients must be a list."
# 			assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
# 			assert eq in [1,2,3], "eq must be either 1, 2, or 3 when using Ramsay (2009)."
# 			if eq == 1:
# 				assert len(Coefficients) == 5, "For Eq. 1, Coefficients must be 5 elements long."
# 				MomentArm = (sp.Matrix(Coefficients,dtype='float64').T\
# 								*sp.Matrix([1,θ,θ**2,θ**3,θ**4]))[0,0]/1000
# 			elif eq == 2:
# 				assert len(Coefficients)==16, "For Eq. 2, Coefficients must be 16 elements long."
# 				MomentArm = (sp.Matrix(Coefficients,dtype='float64').T*\
# 								sp.Matrix([1, θ, θ_PS, θ*θ_PS, θ**2, \
# 											θ_PS**2, (θ**2)*θ_PS, θ*(θ_PS**2), \
# 											(θ**2)*(θ_PS**2), θ**3, θ_PS**3, \
# 											(θ**3)*θ_PS, θ*(θ_PS**3), \
# 											(θ**3)*(θ_PS**2), (θ**2)*(θ_PS**3), \
# 											(θ**3)*(θ_PS**3)]))[0, 0]/1000
# 			else: # eq == 3
# 				assert len(Coefficients)==18, "For Eq. 3, Coefficients must be 18 elements long."
# 				MomentArm = (sp.Matrix(Coefficients,dtype='float64').T*\
# 								sp.Matrix([1, θ, θ_PS, θ*θ_PS, θ**2, \
# 									θ_PS**2, (θ**2)*θ_PS, θ*(θ_PS**2), (θ**2)*(θ_PS**2), \
# 									θ**3, (θ**3)*θ_PS, (θ**3)*(θ_PS**2), \
# 									θ**4, (θ**4)*θ_PS, (θ**4)*(θ_PS**2),  \
# 									θ**5, (θ**5)*θ_PS, (θ**5)*(θ_PS**2)]))[0, 0]/1000
# 		if threshold == None:
# 			return(MomentArm)
# 		else:
# 			assert type(threshold) in [int,float], "threshold must be a number."
# 			MomentArm = sp.Piecewise((MomentArm,θ<threshold),(MomentArm.subs(θ,threshold),θ>=threshold))
# 			return(MomentArm)
#
# 	MuscleList = AllMuscleSettings.keys()
#
# 	RT_symbolic = sp.Matrix([MA_function(AllMuscleSettings[muscle]['Elbow']) for muscle in MuscleList])
# 	dRT_symbolic = sp.Matrix(sp.diff(RT_symbolic,θ_EFE))
# 	d2RT_symbolic = sp.Matrix(sp.diff(sp.diff(RT_symbolic,θ_EFE),θ_EFE))
# 	# RT_func = lambdify([θ_SFE,x1,θ_PS],RT_symbolic)
# 	# dRT_func = lambdify([θ_SFE,x1,θ_PS],dRT_symbolic)
# 	# import ipdb; ipdb.set_trace()
# 	# returns an (n,m) matrix when n is the number of muscles and m is the number of DOFS. We chose to return R.T because this is commonly utilized in muscle velocity calculations.
# 	return(RT_symbolic,dRT_symbolic,d2RT_symbolic)
# def statusbar(i,N,**kwargs):
# 	"""
# 	# i is the current iteration (must be an int) and N is the length of
# 	# the range (must be an int). i must also be in [0,N).
# 	#
# 	# ~~~~~~~~~~~~~~
# 	# **kwargs
# 	# ~~~~~~~~~~~~~~
# 	#
# 	# StartTime should equal time.time() and should be defined before your
# 	# loop to ensure that you get an accurate representation of elapsed time.
# 	#
# 	# Title should be a str that will be displayed before the statusbar. Title
# 	# should be no longer than 25 characters.
# 	#
# 	# ~~~~~~~~~~~~~~
# 	#
# 	# NOTE: you should place a print('\n') after the loop to ensure you
# 	# begin printing on the next line.
#
# 	"""
# 	import time
# 	from scipy import interpolate
# 	import numpy as np
# 	StartTime = kwargs.get("StartTime",False)
# 	Title = kwargs.get("Title",'')
# 	global time_array
# 	global TimeLeft
# 	assert type(i)==int, "i must be an int"
# 	assert type(N)==int, "N must be an int"
# 	assert N>i, "N must be greater than i"
# 	assert N>0, "N must be a positive integer"
# 	assert i>=0, "i must not be negative (can be zero)"
# 	assert type(Title) == str, "Title should be a string"
# 	assert len(Title) <= 22, "Title should be less than 25 characters"
# 	if Title != '' : Title = ' '*(22-len(Title)) + Title + ' : '
# 	statusbar = Title +'[' + '\u25a0'*int((i+1)/(N/50)) + '\u25a1'*(50-int((i+1)/(N/50))) + '] '
# 	TimeBreak = abs
# 	if StartTime != False:
# 		if i==1:
# 			time_array = []
# 			TimeLeft = '--'
# 		elif i==int(0.02*N):
# 			time_array.append(time.time()-StartTime)
# 			TimeLeft = '{0:1.1f}'.format(time_array[-1]*(N/(i+1)))
# 		elif i%int(0.02*N)==0:
# 			time_array.append(time.time()-StartTime)
# 			TimeLeft = '{0:1.1f}'.format(float(interpolate.interp1d(np.arange(len(time_array)),time_array,fill_value='extrapolate')(49))-time_array[-1])
# 		print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete, ' + '{0:1.1f}'.format(time.time() - StartTime) \
# 			+ 'sec, (est. ' + TimeLeft,' sec left)		\r', end='')
# 	else:
# 		print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete           \r',end = '')
#
# t = sy.Symbol('t')
#
# """
# # x_1 &= \theta \\
# # x_2 &= \dot{\theta} \\
# # x_3 &= f_{T,1} \\
# # x_4 &= f_{T,2} \\
# # x_5 &= l_{m,1} \\
# # x_6 &= l_{m,2} \\
# # x_7 &= v_{m,1} \\
# # x_8 &= v_{m,2} \\
# # u_1 &= \alpha_1 \\
# # u_2 &= \alpha_2 \\
#
# """
#
# g,l = 9.80, 0.45 #m/s², m
# M = 2 # kg
# α1 = 0 # 10*np.pi/180 # rads
# α2 = 0 # 10*np.pi/180 # rads
# m1 = 1 # kg
# m2 = 1 # kg
# bm1 = 1 # kg/s
# bm2 = 1 # kg/s
# AllMuscleSettings = return_muscle_settings()
# [[r1,r2],[dr1,dr2],[d2r1,d2r2]]= return_MA_matrix_functions(AllMuscleSettings)
# cT = 27.8
# kT = 0.0047
# LrT = 0.964
# β = 1.55
# ω = 0.75
# ρ = 2.12
# PCSA1 = 100**2*np.pi # mm²
# PCSA2 = 100**2*np.pi # mm²
# V_max = -9.15
# cv0 = -5.78
# cv1 = 9.18
# av0 = -1.53
# av1 = 0
# av2 = 0
# bv = 0.69
# lo1 = AllMuscleSettings["BIC"]["Optimal Muscle Length"]/100
# lo2 = AllMuscleSettings["TRI"]["Optimal Muscle Length"]/100
# c_1 = 23.0
# k_1 = 0.046
# Lr1 = 1.17
# η = 0.01
#
# Amp = 7.5*np.pi/180
# Base = 90*np.pi/180
# Freq = 2*np.pi
#
# k1,k2,k3,k4 = 100,10,10,100
#
# max_allowable_act_jump = 0.3
#
#
# """
# # c_{1} &= -\frac{3g}{2l_{\text{arm}}} \\
# # c_{2} &= \frac{3}{m_{\text{arm}}l_{\text{arm}}^2} \\
# # c_{3} &= \frac{1}{\cos(\alpha_{1})} \\
# # c_{4} &= \frac{1}{\cos(\alpha_{2})} \\
# # c_{5} &= \frac{\cos(\alpha_{1})}{m_1} \\
# # c_{6} &= \frac{\cos^2(\alpha_{1})}{m_1} \\
# # c_{7} &= \frac{b_{m,1}\cos^2(\alpha_{1})}{m_1} \\
# # c_{8} &= \tan^2(\alpha_{1}) \\
# # c_{9} &= \frac{\cos(\alpha_{2})}{m_2} \\
# # c_{10} &= \frac{\cos^2(\alpha_{2})}{m_2} \\
# # c_{11} &= \frac{b_{m,2}\cos^2(\alpha_{2})}{m_2} \\
# # c_{12} &= \tan^2(\alpha_{2}) \\
#
# """
#
# c1 = -(3*g)/(2*l)
# c2 = 3/(M*l**2)
# c3 = 1/(np.cos(α1)*lo1)
# c4 = 1/(np.cos(α2)*lo2)
# c5 = np.cos(α1)/m1
# c6 = np.cos(α1)*c5
# c7 = bm1*c6
# c8 = np.tan(α1)**2
# c9 = np.cos(α2)/m2
# c10 = np.cos(α2)*c9
# c11 = bm2*c10
# c12 = np.tan(α2)**2
#
# '''
# g_{1} &= r_{1}\left(\theta,\theta_{PS}=\frac{\pi}{2}\right) \\
# g_{2} &= r_{2}\left(\theta,\theta_{PS}=\frac{\pi}{2}\right)  \\
# g_{3} &= \kappa^{T}(f_{T,1}) \hspace{1em} \text{(Variable stiffness coefficient from Zajac (1989) ODE for tendon force)} \\
# g_{4} &= \kappa^{T}(f_{T,2}) \hspace{1em} \\
# g_{5} &= \text{sgn}(-r_1(\theta)\cdot\dot{\theta})\cdot\dot{\theta}\cdot\sqrt{\left(\frac{\partial r_1}{\partial\theta}\right)^2 + r_1^2(\theta)} \\
# g_{6} &= \text{sgn}(-r_2(\theta)\cdot\dot{\theta})\cdot\dot{\theta}\cdot\sqrt{\left(\frac{\partial r_2}{\partial\theta}\right)^2 + r_2^2(\theta)} \\
# g_{7} &= F_{\text{MAX},1}\cdot f_{L}(l_{m,1}) \cdot f_{V}(l_{m,1},v_{m,1}) \\
# g_{8} &= F_{PE}(l_{m,1},v_{m,1}) \\
# g_{9} &= F_{PE}(l_{m,2},v_{m,2}) \\
# g_{10} &= F_{\text{MAX},2}\cdot f_{L}(l_{m,2}) \cdot f_{V}(l_{m,2},v_{m,2}) \\
# '''
# FL = lambda l,lo: 1
# FV = lambda l,v,lo: 1
# # FL = lambda l,lo: np.exp(-abs(((l/lo)**β-1)/ω)**ρ)
# # FV = lambda l,v,lo: np.piecewise(v,[v<=0, v>0],\
# # 	[lambda v: (V_max - v/lo)/(V_max + (cv0 + cv1*(l/lo))*(v/lo)),\
# # 	lambda v: (bv-(av0 + av1*(l/lo) + av2*(l/lo)**2)*(v/lo))/(bv + (v/lo))])
#
# def G1(X):
# 	return(lambdify([θ_EFE],r1.subs([(θ_PS,np.pi/2)]))(X[0])) #
# def dG1(X):
# 	return(lambdify([θ_EFE],dr1.subs([(θ_PS,np.pi/2)]))(X[0]))
# def d2G1(X):
# 	return(lambdify([θ_EFE],d2r1.subs([(θ_PS,np.pi/2)]))(X[0]))
# def G2(X):
# 	return(lambdify([θ_EFE],r2.subs([(θ_PS,np.pi/2)]))(X[0])) #
# def dG2(X):
# 	return(lambdify([θ_EFE],dr2.subs([(θ_PS,np.pi/2)]))(X[0]))
# def d2G2(X):
# 	return(lambdify([θ_EFE],d2r2.subs([(θ_PS,np.pi/2)]))(X[0]))
# def G3(X):
# 	return(-(cT/lo1)*(1-np.exp(-X[2]/(cT*kT)))) # Maybe its negative because of theta sign convention
# def dG3(X):
# 	return(np.exp(-X[2]/(cT*kT))/(kT*lo1))
# def G4(X):
# 	return(np.sign(-G1(X))*X[1]*np.sqrt(dG1(X)**2 + G1(X)**2)/lo1) # Should be G5
# def dG4_dt(X):
# 	return(np.sign(-G1(X))*(dX2_dt(X)*np.sqrt(dG1(X)**2 + G1(X)**2) +\
# 		  ((X[1]**2)*dG1(X)*(d2G1(X) + G1(X)))/np.sqrt(dG1(X)**2 + G1(X)**2))/lo1)
# def G5(X):
# 	return((cT/lo2)*(1-np.exp(-X[3]/(cT*kT))))
# def dG5(X):
# 	return(np.exp(-X[3]/(cT*kT))/(kT*lo2))
# def G6(X):
# 	return(np.sign(-G2(X))*X[1]*np.sqrt(dG2(X)**2 + G2(X)**2)/lo2) # Should be G4
# def dG6_dt(X):
# 	return(np.sign(-G2(X))*(dX2_dt(X)*np.sqrt(dG2(X)**2+ G2(X)**2) +\
# 			(X[1]**2*dG2(X)*(d2G2(X) + G2(X)))/np.sqrt(dG2(X)**2 + G2(X)**2))/lo2)
# def G7(X):
# 	return(c_1*k_1*np.log(np.exp(((X[4]/lo1)-Lr1)/k_1) + 1) + η*(X[6]/lo1))  # Should be G8
# def G8(X):
# 	return(0.25*PCSA1*FL(X[4],lo1)*FV(X[4],X[6],lo1)) # Should be G6
# def G9(X):
# 	return(c_1*k_1*np.log(np.exp(((X[5]/lo2)-Lr1)/k_1) + 1) + η*(X[7]/lo2))
# def G10(X):
# 	return(0.25*PCSA2*FL(X[5],lo2)*FV(X[5],X[7],lo2)) # Should be G7
#
# """
# # \dot{x}_1 &= x_{2} \\
# # \dot{x}_2 &= c_{1}\cos(x_{1}) + c_{2}g_{1}x_{3} - c_{2}g_{2}x_{4} \\
# # \dot{x}_3 &= g_{3}g_{5} - c_{3}g_{3}x_{7} \\
# # \dot{x}_4 &= g_{4}g_{6} - c_{4}g_{4}x_{8} \\
# # \dot{x}_5 &= x_{7} \\
# # \dot{x}_6 &= x_{8} \\
# # \dot{x}_7 &= c_{5}x_{3} - c_{6}g_{8} - c_{7}x_{7} + \frac{c_{8}x_{7}^2}{x_{5}} - c_{6}g_{7}u_{1}\\
# # \dot{x}_8 &= c_{9}x_{4} - c_{10}g_{9} - c_{11}x_{8} + \frac{c_{12}x_{8}^2}{x_{6}} - c_{10}g_{10}u_{2} \\
# """
#
# def dX1_dt(X):
# 	return(X[1])
# def dX2_dt(X):
# 	return(c1*np.sin(X[0]) + c2*G1(X)*X[2] + c2*G2(X)*X[3])
# def d2X2_dt2(X):
# 	return(c1*np.cos(X[0])*dX1_dt(X) + c2*dG1(X)*dX1_dt(X)*X[2] + c2*G1(X)*dX3_dt(X)\
# 					+ c2*dG2(X)*dX1_dt(X)*X[3] + c2*G2(X)*dX4_dt(X))
# def dX3_dt(X):
# 	return(G3(X)*G4(X) - c3*G3(X)*X[6])
# def dX4_dt(X):
# 	return(G5(X)*G6(X) - c4*G5(X)*X[7])
# def dX5_dt(X):
# 	return(X[6])
# def dX6_dt(X):
# 	return(X[7])
# def dX7_dt(X,U):
# 	return(c5*X[2] - c6*G8(X)*U[0] - c6*G7(X) - c7*X[6] \
# 	+ c8*X[6]**2/X[4])
# def dX8_dt(X,U):
# 	return(c9*X[3] - c10*G10(X)*U[1] - c10*G9(X) - c11*X[7] \
# 	+ c12*X[7]**2/X[5])
#
# r = lambda t: Amp*np.sin(Freq*t) + Base
# dr = lambda t: Amp*Freq*np.cos(Freq*t)
# d2r = lambda t: -Amp*Freq**2*np.sin(Freq*t)
# d3r = lambda t: -Amp*Freq**3*np.cos(Freq*t)
# d4r = lambda t: Amp*Freq**4*np.sin(Freq*t)
#
# def Z1(t,X):
# 	return(r(t) - X[0])
# def dZ1(t,X):
# 	return(dr(t) - X[1])
# def d2Z1(t,X):
# 	return(d2r(t) - dX2_dt(X))
# def d3Z1(t,X):
# 	return(d3r(t) - d2X2_dt2(X))
# def A1(t,X):
# 	return(dr(t) + k1*Z1(t,X))
# def dA1(t,X):
# 	return(d2r(t) + k1*dZ1(t,X))
# def d2A1(t,X):
# 	return(d3r(t) + k1*d2Z1(t,X))
# def d3A1(t,X):
# 	return(d4r(t) + k1*d3Z1(t,X))
# def Z2(t,X):
# 	return(X[1] - A1(t,X))
# def dZ2(t,X):
# 	return(c1*np.sin(X[0]) + c2*G1(X)*X[2] + c2*G2(X)*X[3] - dA1(t,X))
# def d2Z2(t,X):
# 	dx1_dt = dX1_dt(X)
# 	return(c1*np.cos(X[0])*dx1_dt + c2*dG1(X)*dx1_dt*X[2] + c2*G1(X)*dX3_dt(X) \
# 					+ c2*dG2(X)*dx1_dt*X[3] + c2*G2(X)*dX4_dt(X) - d2A1(t,X))
# def A2(t,X):
# 	return(Z1(t,X) + dA1(t,X) - c1*np.sin(X[0]) - k2*Z2(t,X))
# def dA2(t,X):
# 	return(dZ1(t,X) + d2A1(t,X) - c1*np.cos(X[0])*X[1] - k2*dZ2(t,X))
# def d2A2(t,X):
# 	return(d2Z1(t,X) + d3A1(t,X) + c1*np.sin(X[0])*X[1]**2 - k2*d2Z2(t,X))
# def Z3(t,X):
# 	return(c2*G1(X)*X[2] + c2*G2(X)*X[3] - A2(t,X))
# def dZ3(t,X):
# 	g1 = G1(X)
# 	g2 = G2(X)
# 	g3 = G3(X)
# 	g5 = G5(X)
# 	return(c2*dG1(X)*X[1]*X[2] + c2*dG2(X)*X[1]*X[3] \
# 					+ c2*g1*g3*G4(X) - c2*c3*g1*g3*X[6] \
# 						+ c2*g2*g5*G6(X) - c2*c4*g2*g5*X[7] \
# 							- dA2(t,X))
# def A3(t,X):
# 	return(Z2(t,X) + c2*dG1(X)*X[1]*X[2] + c2*dG2(X)*X[1]*X[3] \
# 					+ c2*G1(X)*G3(X)*G4(X) + c2*G2(X)*G5(X)*G6(X) \
# 						- dA2(t,X) + k3*Z3(t,X))
# def dA3(t,X):
# 	g1 = G1(X)
# 	dg1 = dG1(X)
# 	g2 = G2(X)
# 	dg2 = dG2(X)
# 	dx1_dt = dX1_dt(X)
# 	dx2_dt = dX2_dt(X)
# 	dx3_dt = dX3_dt(X)
# 	dx4_dt = dX4_dt(X)
# 	g3 = G3(X)
# 	g4 = G4(X)
# 	g5 = G5(X)
# 	g6 = G6(X)
# 	return(dZ2(t,X) + c2*d2G1(X)*(X[1]**2)*X[2] + c2*dg1*dx2_dt*X[2]\
# 	 		+ c2*dg1*X[1]*dx3_dt + c2*d2G2(X)*(X[1]**2)*X[3] \
# 				+ c2*dg2*dx2_dt*X[3] + c2*dg2*X[1]*dx4_dt \
# 					+ c2*dg1*dx1_dt*g3*g4 + c2*g1*dG3(X)*dx3_dt*g4 \
# 						+ c2*g1*g3*dG4_dt(X) \
# 							+ c2*dg2*dx1_dt*g5*g6 + c2*g2*dG5(X)*dx4_dt*g6 \
# 								+ c2*g2*g5*dG6_dt(X) \
# 									- d2A2(t,X) + k3*dZ3(t,X))
# def Z4(t,X):
# 	return(c2*c3*G1(X)*G3(X)*X[6] + c2*c4*G2(X)*G5(X)*X[7] - A3(t,X))
# def dZ4(t,X,U):
# 	g1 = G1(X)
# 	g2 = G2(X)
# 	g3 = G3(X)
# 	dg3 = dG3(X)
# 	g5 = G5(X)
# 	dg5 = dG5(X)
# 	return(c2*c3*X[6]*(dG1(X)*g3*X[1] + g1*dg3*g3*G4(X) - c7*g1*g3) \
# 			+ c2*c4*X[7]*(dG2(X)*g5*X[1] + g2*dg5*g5*G6(X) - c11*g2*g5) \
# 			+ c2*c3*X[6]**2*(c8*g1*g3/X[4] - c3*g1*dg3*g3) \
# 			+ c2*c4*X[7]**2*(c12*g2*g5/X[5] - c4*g2*dg5*g5) \
# 			+ c2*c3*c5*g1*g3*X[2] - c2*c3*c6*g1*g3*G7(X) \
# 			+ c2*c4*c9*g2*g5*X[3] - c2*c4*c10*g2*g5*G9(X) \
# 			- c2*c3*c6*g1*g3*G8(X)*U[0] - c2*c4*c10*g2*g5*G10(X)*U[1] \
# 			-dA3(t,X))
# Coefficient_11 = lambda t,X: c2*c3*c6*G1(X)*G3(X)*G8(X)
# Coefficient_12 = lambda t,X: c2*c4*c10*G2(X)*G5(X)*G10(X)
# def Constraint_1(t,X):
# 	g1 = G1(X)
# 	g2 = G2(X)
# 	g3 = G3(X)
# 	g5 = G5(X)
# 	dx1_dt = dX1_dt(X)
# 	return(	c2*c3*dG1(X)*dx1_dt*g3*X[6] + c2*c3*g1*dG3(X)*dX3_dt(X)*X[6] \
# 			+ c2*c4*dG2(X)*dx1_dt*g5*X[7] + c2*c4*g2*dG5(X)*dX4_dt(X)*X[7] \
# 			+ c2*c3*c5*g1*g3*X[2] - c2*c3*c6*g1*g3*G7(X) - c2*c3*c7*g1*g3*X[6] \
# 			+ c2*c4*c9*g2*g5*X[3] - c2*c4*c10*g2*g5*G9(X) - c2*c4*c11*g2*g5*X[7] \
# 			+ c2*c3*c8*g1*g3*((X[6]**2)/X[4]) + c2*c4*c12*g2*g5*((X[7]**2)/X[5]) \
# 			- dA3(t,X) - Z3(t,X) + k4*Z4(t,X))
# 					# = c2*c3*c6*G1(X)*G3(X)*G8(X)*U[0] + c2*c4*c10*G2(X)*G5(X)*G10(X)*U[1]
# # Constraint_1 = lambda t,X: \
# # 				c2*c3*X[6]*(dG1(X)*G3(X)*X[1] + G1(X)*dG3(X)*G3(X)*G4(X) - c7*G1(X)*G3(X)) \
# # 				+ c2*c4*X[7]*(dG2(X)*G5(X)*X[1] + G2(X)*dG5(X)*G5(X)*G6(X) - c11*G2(X)*G5(X)) \
# # 				+ c2*c3*X[6]**2*(c8*G1(X)*G3(X)/X[4] - c3*G1(X)*dG3(X)*G3(X)) \
# # 				+ c2*c4*X[7]**2*(c12*G2(X)*G5(X)/X[5] - c4*G2(X)*dG5(X)*G5(X)) \
# # 				+ c2*c3*c5*G1(X)*G3(X)*X[2] - c2*c3*c6*G1(X)*G3(X)*G7(X) \
# # 				+ c2*c4*c9*G2(X)*G5(X)*X[3] - c2*c4*c10*G2(X)*G5(X)*G9(X) \
# # 				-dA3(t,X) - Z3(t,X) + k4*Z4(t,X) \
# # 					# = c2*c3*c6*G1(X)*G3(X)*G8(X)*U[0] + c2*c4*c10*G2(X)*G5(X)*G10(X)*U[1]
# # Z1 = lambda t,X: r(t) - X[0]
# # dZ1 = lambda t,X: dr(t) - X[1]
# # d2Z1 = lambda t,X: d2r(t) - dX2_dt(X)
# # d3Z1 = lambda t,X: d3r(t) - d2X2_dt2(X)
# # a1 = lambda t,X: dr(t) - k1*Z1(t,X)
# # da1 = lambda t,X: d2r(t) - k1*dZ1(t,X)
# # d2a1 = lambda t,X: d3r(t) - k1*d2Z1(t,X)
# # d3a1 = lambda t,X: d4r(t) - k1*d3Z1(t,X)
# # z2 = lambda t,X: X[1] - A1(t,X)
# # dz2 = lambda t,X: c1*np.cos(X[0]) + c2*G1(X)*X[2] - c2*G2(X)*X[3] - dA1(t,X)
# # d2z2 = lambda t,X: -c1*np.sin(X[0])*dX1_dt(X) + c2*dG1(X)*dX1_dt(X)*X[2] + c2*G1(X)*dX3_dt(X) \
# # 					- c2*dG2(X)*dX1_dt(X)*X[3] - c2*G2(X)*dX4_dt(X) - d2A1(t,X)
# # a2 = lambda t,X: Z1(t,X) + dA1(t,X) - c1*np.cos(X[0]) - k2*Z2(t,X)
# # da2 = lambda t,X: dZ1(t,X) + d2A1(t,X) + c1*np.sin(X[0])*X[1] - k2*dZ2(t,X)
# # d2a2 = lambda t,X: d2Z1(t,X) + d3A1(t,X) + c1*np.cos(X[0])*X[1]**2 - k2*d2Z2(t,X)
# # z3 = lambda t,X: c2*G1(X)*X[2] - c2*G2(X)*X[3] - A2(t,X)
# # dz3 = lambda t,X: c2*dG1(X)*X[1]*X[2] - c2*dG2(X)*X[1]*X[3] \
# # 					+ c2*G1(X)*G3(X)*G4(X) - c2*c3*G1(X)*G3(X)*X[6] \
# # 						-c2*G2(X)*G5(X)*G6(X) + c2*c4*G2(X)*G5(X)*X[7] \
# # 							- dA2(t,X)
# # a3 = lambda t,X: Z2(t,X) + c2*dG1(X)*X[1]*X[2] - c2*dG2(X)*X[1]*X[3] \
# # 					+ c2*G1(X)*G3(X)*G4(X) - c2*G2(X)*G5(X)*G6(X) \
# # 						- dA2(t,X) + k3*Z3(t,X)
# # da3 = lambda t,X: dZ2(t,X) + c2*d2G1(X)*X[1]**2*X[2] + c2*dG1(X)*dX2_dt(X)*X[2]\
# # 			 		+ c2*dG1(X)*X[1]*dX3_dt(X) - c2*d2G2(X)*X[1]**2*X[3] \
# # 						- c2*dG2(X)*dX2_dt(X)*X[3] - c2*dG2(X)*X[1]*dX4_dt(X) \
# # 							+ c2*dG1(X)*dX1_dt(X)*G3(X)*G4(X) + c2*G1(X)*dG3(X)*dX3_dt(X)*G4(X) \
# # 								+ c2*G1(X)*G3(X)*dG4_dt(X) \
# # 									- c2*dG2(X)*dX1_dt(X)*G5(X)*G6(X) - c2*G2(X)*dG5(X)*dX4_dt(X)*G6(X) \
# # 										- c2*G2(X)*G5(X)*dG6_dt(X) \
# # 											- d2A2(t,X) + k3*dZ3(t,X)
# # z4 = lambda t,X: c2*c3*G1(X)*G3(X)*X[6] - c2*c4*G2(X)*G5(X)*X[7] - A3(t,X)
# # dz4 = lambda X,U,t: \
# # 				c2*c3*X[6]*(dG1(X)*G3(X)*X[1] + G1(X)*dG3(X)*G3(X)*G4(X) - c7*G1(X)*G3(X)) \
# # 				- c2*c4*X[7]*(dG2(X)*G5(X)*X[1] + G2(X)*dG5(X)*G5(X)*G6(X) - c11*G2(X)*G5(X)) \
# # 				+ c2*c3*X[6]**2*(c8*G1(X)*G3(X)/X[4] - c3*G1(X)*dG3(X)*G3(X)) \
# # 				- c2*c4*X[7]**2*(c12*G2(X)*G5(X)/X[5] - c4*G2(X)*dG5(X)*G5(X)) \
# # 				+ c2*c3*c5*G1(X)*G3(X)*X[2] - c2*c3*c6*G1(X)*G3(X)*G7(X) \
# # 				- c2*c4*c9*G2(X)*G5(X)*X[3] + c2*c4*c10*G2(X)*G5(X)*G9(X) \
# # 				- c2*c3*c6*G1(X)*G3(X)*G8(X)*U[0] + c2*c4*c10*G2(X)*G5(X)*G10(X)*U[1] \
# # 				-dA3(t,X)
# # Coefficient_11 = lambda t,X: c2*c3*c6*G1(X)*G3(X)*G8(X)
# # Coefficient_12 = lambda t,X: -c2*c4*c10*G2(X)*G5(X)*G10(X)
# # Constraint_1 = lambda t,X: \
# # 				c2*c3*X[6]*(dG1(X)*G3(X)*X[1] + G1(X)*dG3(X)*G3(X)*G4(X) - c7*G1(X)*G3(X)) \
# # 				- c2*c4*X[7]*(dG2(X)*G5(X)*X[1] + G2(X)*dG5(X)*G5(X)*G6(X) - c11*G2(X)*G5(X)) \
# # 				+ c2*c3*X[6]**2*(c8*G1(X)*G3(X)/X[4] - c3*G1(X)*dG3(X)*G3(X)) \
# # 				- c2*c4*X[7]**2*(c12*G2(X)*G5(X)/X[5] - c4*G2(X)*dG5(X)*G5(X)) \
# # 				+ c2*c3*c5*G1(X)*G3(X)*X[2] - c2*c3*c6*G1(X)*G3(X)*G7(X) \
# # 				- c2*c4*c9*G2(X)*G5(X)*X[3] + c2*c4*c10*G2(X)*G5(X)*G9(X) \
# # 				-dA3(t,X) - Z3(t,X) + k4*Z4(t,X) \
# # 					# = c2*c3*c6*G1(X)*G3(X)*G8(X)*U[0] - c2*c4*c10*G2(X)*G5(X)*G10(X)*U[1]
#
# def return_U(t,X,U,e,Noise):
# 	import random
# 	Constraint1 = Constraint_1(t,X)
# 	d1 =Coefficient_11(t,X)
# 	d2 = Coefficient_12(t,X)
# 	# assert d1>=0, "Error with 1st Coefficient. Should be positive."
# 	# assert d2<=0, "Error with 2nd Coefficient. Should be negative."
# 	# assert Constraint1>d2, "HUGE ERROR... Constraints create infeasible activations (i.e., negative)"
# 	Feasible = True
# 	if Constraint1 >=0:
# 		# assert d1>0 or d2>0, "HUGE PROBLEM!! Infeasible activations!"
# 		if d1 > 0 and d2 > 0:
# 			lowerbound = 0
# 			if Constraint1<d1:
# 				upperbound = Constraint1/d1
# 			else:
# 				upperbound = 1
# 		elif d1 > 0 and d2 < 0:
# 			# assert Constraint1 < d1, "HUGE ERROR, infeasible activations"
# 			if Constraint1<d1:
# 				lowerbound = Constraint1/d1
# 				if Constraint1 < d1 + d2:
# 					upperbound = (Constraint1 - d2)/d1
# 				else:
# 					upperbound = 1
# 			else:
# 				Feasible = False
# 		elif d1 < 0 and d2 > 0:
# 			# assert Constraint1 < d2, "HUGE ERROR, infeasible activations"
# 			if Constraint1 < d2:
# 				lowerbound = 0
# 				if Constraint1 > d1 + d2:
# 					upperbound = (Constraint1 - d2)/d1
# 				else:
# 					upperbound = 1
# 			else:
# 				Feasible = False
# 		else:
# 			Feasible = False
# 	else: #Constraint1 <0
# 		# assert d1<0 or d2<0, "HUGE PROBLEM!! Infeasible activations!"
# 		if d1 < 0 and d2 < 0:
# 			lowerbound = 0
# 			if Constraint1>d1:
# 				upperbound = Constraint1/d1
# 			else:
# 				upperbound = 1
# 		elif d1 < 0 and d2 > 0:
# 			# assert Constraint1 > d1, "HUGE ERROR, infeasible activations"
# 			if Constraint1 > d1:
# 				lowerbound = Constraint1/d1
# 				if Constraint1 > d1 + d2:
# 					upperbound = (Constraint1 - d2)/d1
# 				else:
# 					upperbound = 1
# 			else:
# 				Feasible = False
# 		elif d1 > 0 and d2 < 0:
# 			# assert Constraint1 > d2, "HUGE ERROR, infeasible activations"
# 			if Constraint1 > d2:
# 				lowerbound = 0
# 				if Constraint1 < d1 + d2:
# 					upperbound = (Constraint1 - d2)/d1
# 				else:
# 					upperbound = 1
# 			else:
# 				Feasible = False
# 		else:
# 			Feasible = False
# 	if Feasible == True:
# 		poss_u1 = (upperbound-lowerbound)*np.random.rand(1000) + lowerbound
# 		poss_u2 = np.array([Constraint1/d2 - (d1/d2)*el for el in poss_u1])
# 		euclid_dist = np.array(list(map(lambda u1,u2: np.sqrt((U[0]-u1)**2 + (U[1]-u2)**2),\
# 								poss_u1,poss_u2)))
# 		feasible_index = np.where(np.logical_and(euclid_dist>=0, euclid_dist<=max_allowable_act_jump))
# 		next_index = random.choice(feasible_index[0])
# 		u1 = poss_u1[next_index]
# 		u2 = poss_u2[next_index]
# 	else:
# 		u1 = U[0]
# 		u2 = U[1]
# 	# poss_u1 = poss_u1[feasible_index]
# 	# poss_u2 = poss_u2[feasible_index]
# 	# assert len(poss_u1)!=0, "Error! No feasible activations. Perhaps a finer sample is needed."
# 	# import ipdb; ipdb.set_trace()
# 	# if Constraint1<=0:
# 	#     u1 = Constraint1/(Coefficient_11(t,X) + e*Coefficient_12(t,X)) + Noise[0]
# 	#     u2 = e*(u1-Noise[0]) + Noise[1]
# 	# else:
# 	#     u2 = Constraint1/(e*Coefficient_11(t,X) + Coefficient_12(t,X)) + Noise[1]
# 	#     u1 = e*(u2-Noise[1]) + Noise[0]
# 	return([u1,u2])
#
# N = 10001
# Time = np.linspace(0,10,N)
# dt = Time[1]-Time[0]
# x1,x2,x3,x4,x5,x6,x7,x8 = [Base],[Amp*Freq],[100],[100],[lo1],[lo2],[0],[0]
# u1,u2 = [0.1],[0.1]
# CocontractionIndex = 2
#
# AddNoise = False
# if AddNoise == True:
#     np.random.seed(seed=None)
#     NoiseArray = np.random.normal(loc=0.0,scale=0.2,size=(2,len(Time)))
# else:
#     NoiseArray = np.zeros((2,len(Time)))
#
# def update_policy(t,x1,x2,x3,x4,x5,x6,x7,x8,dt,NoiseArray,e=2):
# 	import numpy as np
# 	X = [x1[-1],x2[-1],x3[-1],x4[-1],x5[-1],x6[-1],x7[-1],x8[-1]]
# 	U = [u1[-1],u2[-1]]
# 	U = return_U(t,X,U,e,NoiseArray[:,int(t/dt)])
# 	u1.append(U[0])
# 	u2.append(U[1])
# 	x8.append(x8[-1] + dX8_dt(X,U)*dt)
# 	x7.append(x7[-1] + dX7_dt(X,U)*dt)
# 	x6.append(x6[-1] + dX6_dt(X)*dt)
# 	x5.append(x5[-1] + dX5_dt(X)*dt)
# 	x4.append(x4[-1] + dX4_dt(X)*dt)
# 	x3.append(x3[-1] + dX3_dt(X)*dt)
# 	x2.append(x2[-1] + dX2_dt(X)*dt)
# 	x1.append(x1[-1] + dX1_dt(X)*dt)
#
#
# StartTime = time.time()
# for t in Time[1:]:
# 	update_policy(t,x1,x2,x3,x4,x5,x6,x7,x8,dt,NoiseArray,e=CocontractionIndex)
# 	statusbar(int(t/dt),len(Time[1:]),StartTime=StartTime,Title="mm Driven Pendulum")
#
# plt.figure()
# # plt.title(r'$\dot{x}_{1} = x_{1}^{2} + x_{2}; \hspace{1em} \dot{x}_{2} = u$',\
# #                 fontsize=16,color='gray')
# plt.title("Underdetermined Muscle-Driven\nForced Pendulum Example",\
#                 fontsize=16,color='gray')
# plt.plot(Time,x1,'b',lw=2)
# plt.plot(Time,r(Time),'r--')
# plt.xlabel("Time (s)")
# plt.ylabel("Desired Measure")
# plt.legend([r"Output $y = x_{1}$",r"Reference $r(t) = " + str(Amp) + "\sin(" + str(Freq) + "t)$ " + str(Base)])
#
# plt.figure()
# plt.title('Error vs. Time')
# plt.plot(Time, r(Time)-x1,color='r')
# plt.xlabel("Time (s)")
# plt.ylabel("Error")
#
# plt.figure()
# plt.plot(Time[-1:],u1,'g',Time[-1:],u2,'r')
# plt.title('Muscle Activations vs. Time')
# plt.xlabel("Time (s)")
# plt.ylabel("Muscle Activations")
# plt.legend(["Muscle 1","Muscle 2"])
#
# # Constraint_1 = Constraint_1.subs([(sy.Derivative(r,t),dr),(sy.Derivative(dr,t),d2r)\
# # 									,(sy.Derivative(d2r,t),d3r),(sy.Derivative(d3r,t),d4r)])
# #
# # ConstraintFunction_1 = lambdify([x1,x2,x3,x4,x5,x6,x7,x8,r,dr,d2r,d3r,d4r],Constraint_1)
# # Coefficient11Function = lambdify([x1,x2,x3,x4,x5,x6,x7,x8],Coefficient_11)
# # Coefficient12Function = lambdify([x1,x2,x3,x4,x5,x6,x7,x8],Coefficient_12)
#
# """
