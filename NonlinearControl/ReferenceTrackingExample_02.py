import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import sympy as sym
def return_muscle_settings():
	"""
	Notes:
	Coefficients from observation, Ramsay, FVC, Holtzbaur, Pigeon, Kuechle, or Banks.

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

		Note that the coefficients listed in Pigeon (1996) are given in decending order (i.e., c₅,c₄,c₃,c₂,c₁,c₀). However to maintain continuity with the equations given in Ramsay (2009), we list coefficients in order of increasing power (i.e., c₀,c₁,c₂,c₃,c₄,c₅).
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
    						*np.matrix([1,θ,θ**2,θ**3,θ**4,θ**5]).T)[0,0]
    	elif src.capitalize() == 'Est' :
    		MomentArm = np.array(Coefficients,dtype='float64')
    	else: #src.capitalize() == 'Ramsay'
    		θ = θ_EFE
    		assert type(Coefficients) == list, "Coefficients must be a list."
    		assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
    		assert eq in [1,2,3], "eq must be either 1, 2, or 3 when using Ramsay (2009)."
    		if eq == 1:
    			assert len(Coefficients) == 5, "For Eq. 1, Coefficients must be 5 elements long."
    			MomentArm = (sp.Matrix(Coefficients,dtype='float64').T\
    							*sp.Matrix([1,θ,θ**2,θ**3,θ**4]))[0,0]
    		elif eq == 2:
    			assert len(Coefficients)==16, "For Eq. 2, Coefficients must be 16 elements long."
    			MomentArm = (sp.Matrix(Coefficients,dtype='float64').T*\
    							sp.Matrix([1, θ, θ_PS, θ*θ_PS, θ**2, \
    										θ_PS**2, (θ**2)*θ_PS, θ*(θ_PS**2), \
    										(θ**2)*(θ_PS**2), θ**3, θ_PS**3, \
    										(θ**3)*θ_PS, θ*(θ_PS**3), \
    										(θ**3)*(θ_PS**2), (θ**2)*(θ_PS**3), \
    										(θ**3)*(θ_PS**3)]))[0, 0]
    		else: # eq == 3
    			assert len(Coefficients)==18, "For Eq. 3, Coefficients must be 18 elements long."
    			MomentArm = (sp.Matrix(Coefficients,dtype='float64').T*\
    							sp.Matrix([1, θ, θ_PS, θ*θ_PS, θ**2, \
    								θ_PS**2, (θ**2)*θ_PS, θ*(θ_PS**2), (θ**2)*(θ_PS**2), \
    								θ**3, (θ**3)*θ_PS, (θ**3)*(θ_PS**2), \
    								θ**4, (θ**4)*θ_PS, (θ**4)*(θ_PS**2),  \
    								θ**5, (θ**5)*θ_PS, (θ**5)*(θ_PS**2)]))[0, 0]
    	if threshold == None:
    		return(MomentArm)
    	else:
    		assert type(threshold) in [int,float], "threshold must be a number."
    		MomentArm = sp.Piecewise((MomentArm,θ<threshold),(MomentArm.subs(θ,threshold),θ>=threshold))
    		return(MomentArm)

    MuscleList = AllMuscleSettings.keys()

    Rᵀ_symbolic = sp.Matrix([MA_function(AllMuscleSettings[muscle]['Elbow']) for muscle in MuscleList])
    Ṙᵀ_symbolic = sp.Matrix(sp.diff(Rᵀ_symbolic,θ_EFE))
    # Rᵀ_func = lambdify([θ_SFE,x1,θ_PS],Rᵀ_symbolic)
    # Ṙᵀ_func = lambdify([θ_SFE,x1,θ_PS],Ṙᵀ_symbolic)
    # import ipdb; ipdb.set_trace()
    # returns an (n,m) matrix when n is the number of muscles and m is the number of DOFS. We chose to return R.T because this is commonly utilized in muscle velocity calculations.
    return(Rᵀ_symbolic,Ṙᵀ_symbolic)


# G,Ġ
t = sym.Symbol('t')
x1 = sym.Function('x1')(t)
x2 = sym.Function('x2')(t)
x3 = sym.Function('x3')(t)
x4 = sym.Function('x4')(t)
x5 = sym.Function('x5')(t)
x6 = sym.Function('x6')(t)
x7 = sym.Function('x7')(t)
x8 = sym.Function('x8')(t)

AllMuscleSettings = return_muscle_settings()
[[r1,r2],_]= return_MA_matrix_functions(AllMuscleSettings)
g1 = r1.subs([(θ_EFE,x1),(θ_PS,np.pi/2)])
g2 = r2.subs([(θ_EFE,x1),(θ_PS,np.pi/2)])
# g1 =
# g2
# g3
# g4
# g5
# g6
# g7
# g8


ẋ1 = x2
ẋ2 = x2
ẋ3 = x2
ẋ4 = x2
ẋ5 = x2
ẋ6 = x2
ẋ7 = x2
ẋ8 = x2

# g,l = 9.80, 0.40
# k1,k2 = 50,50
# m1,m2,M = 1,1,1
# A,w = 0.10,0.5*np.pi
# b1,b2,b3,b4 = 20,20,20,20
# CocontractionIndex = 2
#
# def dx1(t,X):
#     return(X[1])
# def dx2(t,X):
#     return(-(k1+k2)/M*X[0] + (k1/M)*X[2] + (k2/M)*X[3])
# def dx3(t,X):
#     return(X[4])
# def dx4(t,X):
#     return(X[5])
# def dx5(t,X,U):
#     return(k1/m1*X[0] -k1/m1*X[2] + U[0]/m1)
# def dx6(t,X,U):
#     return(k2/m2*X[0] -k2/m2*X[3] - U[1]/m2)
#
# r = lambda t: A*np.sin(w*t)
# dr = lambda t: A*w*np.cos(w*t)
# d2r = lambda t: -A*w**2*np.sin(w*t)
# d3r = lambda t: -A*w**3*np.cos(w*t)
# d4r = lambda t: A*w**4*np.sin(w*t)
# def z1(t,X):
#     return(r(t) - X[0])
# def dz1(t,X):
#     return(dr(t) - X[1])
# def d2z1(t,X):
#     return(d2r(t) - dx2(t,X))
# def d3z1(t,X):
#     return(d3r(t) + (k1+k2)/M*X[1] - k1/M*X[4] - k2/M*X[5])
# def a1(t,X):
#     return(dr(t) - b1*z1(t,X))
# def da1(t,X):
#     return(d2r(t) - b1*dz1(t,X))
# def z2(t,X):
#     return(X[1] - a1(t,X))
# def dz2(t,X):
#     return(dx2(t,X) - da1(t,X))
# def a2(t,X):
#     return((k1+k2)/M*X[0] + (1 + b1*b2)*z1(t,X) + (b1+b2)*dz1(t,X) + d2r(t))
# def da2(t,X):
#     return((k1+k2)/M*X[1] + (1 + b1*b2)*dz1(t,X) + (b1+b2)*d2z1(t,X) + d3r(t))
# def d2a2(t,X):
#     return((k1+k2)/M*dx2(t,X) + (1 + b1*b2)*d2z1(t,X) + (b1+b2)*d3z1(t,X) + d4r(t))
# def z3(t,X):
#     return(k1/M*X[2] + k2/M*X[3] - a2(t,X))
# def dz3(t,X):
#     return(k1/M*X[4] + k2/M*X[5] - da2(t,X))
# def a3(t,X):
#     return(da2(t,X) - z2(t,X) -b3*z3(t,X))
# def da3(t,X):
#     return(d2a2(t,X) - dz2(t,X) -b3*dz3(t,X))
# def z4(t,X):
#     return(k1/M*X[4] + k2/M*X[5] - a3(t,X))
# # def dz4(t,X):
# #     return(k1/M*dx5(t,X,U) + k2/M*dx6(t,X,U) - da3(t,X))
# def c1(t,X):
#     return(-(k1**2*m2 + k2**2*m1)/(m1*m2*M)*X[0] + k1**2/(m1*M)*X[2] + k2**2/(m2*M)*X[3] + \
#             da3(t,X) - z3(t,X) - b4*z4(t,X))
# def return_U(t,X,e,Noise):
#     if c1(t,X)<=0:
#         u1 = ((m1*m2*M)/(k1*m2-e*k2*m1))*c1(t,X) + Noise[0]
#         u2 = e*(u1-Noise[0]) + Noise[1]
#     else:
#         u2 = ((m1*m2*M)/(e*k1*m2-k2*m1))*c1(t,X) + Noise[1]
#         u1 = e*(u2-Noise[1]) + Noise[0]
#     return([u1,u2])
# def animate_trajectory(response,Time,x1,x3,x4,u1,u2):
#     assert type(response)==bool, "Input must be either True or False."
#
#     if response == True:
#         import numpy as np
#         import matplotlib.pyplot as plt
#         from matplotlib.patches import Ellipse
#         import matplotlib.animation as animation
#         import matplotlib.patches as patches
#         from scipy import signal
#
#         fig = plt.figure(figsize=(10,8))
#         ax1 = plt.subplot2grid((3,4),(0,0),colspan=4)
#         ax2 = plt.subplot2grid((3,4),(1,0),colspan=2)
#         ax3 = plt.subplot2grid((3,4),(1,2),colspan=2)
#         ax4 = plt.subplot2grid((3,4),(2,0),colspan=3)
#         ax5 = plt.subplot2grid((3,4),(2,3))
#
#         plt.suptitle("Underdetermined Mass-Spring System",Fontsize=28,y=0.95)
#
#         # Model Drawing
#         IdealBoxScalingFactor = 0.78533496170320571 # Calculated from w = np.pi
#         CurrentTrialScalingFactor = max([max(x3)-min(x1),max(x1)-min(x4)])
#         StraightLength = 0.05*CurrentTrialScalingFactor/IdealBoxScalingFactor
#         RestingLength = max([max(x1)-min(x3),max(x4)-min(x1)])+2*StraightLength\
#                         +0.30*CurrentTrialScalingFactor/IdealBoxScalingFactor
#         CenterBoxHalfWidth = 0.15*CurrentTrialScalingFactor/IdealBoxScalingFactor
#         CenterBoxHalfHeight = 0.2*CurrentTrialScalingFactor/IdealBoxScalingFactor
#         SideBoxHalfWidth = 0.1*CurrentTrialScalingFactor/IdealBoxScalingFactor
#         SideBoxHalfHeight = 0.075*CurrentTrialScalingFactor/IdealBoxScalingFactor
#         ForceScaling = 1*CurrentTrialScalingFactor/IdealBoxScalingFactor
#
#         Spring_array =\
#          SideBoxHalfWidth\
#             *np.abs(signal.sawtooth(5*2*np.pi*np.linspace(0,1,1001)-np.pi/2))\
#                 -(1/2)*SideBoxHalfWidth
#
#         Spring1, =\
#             ax1.plot(np.linspace(x1[0]+CenterBoxHalfWidth+StraightLength,\
#                                     RestingLength+x3[0]-SideBoxHalfWidth-StraightLength,1001),\
#                                         Spring_array,'k')
#         Spring1_left, = ax1.plot([x1[0]+CenterBoxHalfWidth,x1[0]+CenterBoxHalfWidth+StraightLength],[0,0],'k')
#         Spring1_right, = \
#             ax1.plot([RestingLength+x3[0]-SideBoxHalfWidth-StraightLength,\
#                         RestingLength+x3[0]-SideBoxHalfWidth],\
#                             [0,0],'k')
#
#         Spring2, =\
#             ax1.plot(np.linspace(-RestingLength+x4[0]+SideBoxHalfWidth+StraightLength,\
#                                     x1[0]-CenterBoxHalfWidth-StraightLength,1001),\
#                                         Spring_array,'k')
#         Spring2_left, = ax1.plot([x1[0]-CenterBoxHalfWidth-StraightLength,x1[0]-CenterBoxHalfWidth],[0,0],'k')
#         Spring2_right, = \
#             ax1.plot([-RestingLength+x4[0]+SideBoxHalfWidth,\
#                         -RestingLength+x4[0]+SideBoxHalfWidth+StraightLength],\
#                             [0,0],'k')
#         ax1.get_xaxis().set_ticks([])
#         ax1.get_yaxis().set_ticks([])
#         ax1.set_frame_on(True)
#         CenterMass = plt.Rectangle((-CenterBoxHalfWidth,-CenterBoxHalfHeight),\
#                                     2*CenterBoxHalfWidth,2*CenterBoxHalfHeight,Color='#4682b4')
#         ax1.add_patch(CenterMass)
#         Mass1 = plt.Rectangle((-SideBoxHalfWidth+RestingLength,-SideBoxHalfHeight),\
#                                     2*SideBoxHalfWidth,2*SideBoxHalfHeight,Color='#4682b4')
#         ax1.add_patch(Mass1)
#         Mass2 = plt.Rectangle((-SideBoxHalfWidth-RestingLength,-SideBoxHalfHeight),\
#                                     2*SideBoxHalfWidth,2*SideBoxHalfHeight,Color='#4682b4')
#         ax1.add_patch(Mass2)
#
#         PositionArrow, = ax1.plot([x1[0],x1[0]],[0,2*CenterBoxHalfHeight],'k')
#         PositionArrowHead, = ax1.plot([x1[0]],[2*CenterBoxHalfHeight],'k^')
#         PositionArrowTail, = ax1.plot([x1[0]],[0],'ko')
#
#         Scale = ax1.plot([-1.1*A,1.1*A],\
#                             [2.75*CenterBoxHalfHeight,2.75*CenterBoxHalfHeight],\
#                                 '0.60')
#         Ticks = np.linspace(-A,A,5)
#         TickHeights = [0.3*CenterBoxHalfHeight,\
#                         0.15*CenterBoxHalfHeight,\
#                         0.3*CenterBoxHalfHeight,\
#                         0.15*CenterBoxHalfHeight,\
#                         0.3*CenterBoxHalfHeight]
#         [ax1.plot([Ticks[i],Ticks[i]],\
#                 [2.75*CenterBoxHalfHeight-TickHeights[i],2.75*CenterBoxHalfHeight],'0.60') \
#                     for i in range(5)]
#
#         Force1Arrow, = ax1.plot([RestingLength+x3[0]+(5/3)*SideBoxHalfWidth,\
#                                     RestingLength + x3[0]+(5/3)*SideBoxHalfWidth\
#                                         +ForceScaling*u1[0]/(max(u1[5000:]+u2[5000:]))],\
#                                             [0,0],'g')
#         Force1ArrowHead, = \
#             ax1.plot([RestingLength + x3[0]+(5/3)*SideBoxHalfWidth\
#                         +ForceScaling*u1[0]/(max(u1[5000:]+u2[5000:]))],[0],'g>')
#         Force2Arrow, =\
#             ax1.plot([x4[0]-RestingLength-(5/3)*SideBoxHalfWidth\
#                         -ForceScaling*u2[0]/(max(u1[5000:]+u2[5000:])),\
#                             x4[0]-RestingLength-(5/3)*SideBoxHalfWidth],[0,0],'r')
#         Force2ArrowHead, = \
#             ax1.plot([x4[0]-RestingLength-(5/3)*SideBoxHalfWidth\
#                         -ForceScaling*u2[0]/(max(u1[5000:]+u2[5000:]))],[0],'r<')
#
#         LowerBound = (np.array(x4[5001:])-RestingLength-(5/3)*SideBoxHalfWidth\
#                         -ForceScaling*np.array(u2[5000:])/(max(u1[5000:]+u2[5000:]))).min()
#         UpperBound = (RestingLength + np.array(x3[5001:])+(5/3)*SideBoxHalfWidth\
#                         +ForceScaling*np.array(u1[5000:])/(max(u1[5000:]+u2[5000:]))).max()
#         Bound = 1.05*np.array([-LowerBound,UpperBound]).max()
#         ax1.set_xlim([-Bound,Bound])
#         ax1.set_ylim([-1.5*CenterBoxHalfHeight,3.25*CenterBoxHalfHeight])
#         ax1.set_aspect('equal')
#
#         #Force 1
#
#         Force1, = ax3.plot([0],[u1[0]],color = 'g')
#         ax3.set_xlim(0,Time[-1])
#         ax3.set_xticks(list(np.linspace(0,Time[-1],5)))
#         ax3.set_xticklabels([str(0),'','','',str(Time[-1])])
#         ax3.set_ylim(0,1.15*max(u1[5000:]+u2[5000:]))
#         if np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                         int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1).shape[0] < 5:
#             ax3.set_yticks(list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                             int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1)))
#             ax3.set_yticklabels([""]*(int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1))
#         else:
#             NumTicks = np.floor(1.15*max(u1[5000:]+u2[5000:]))
#             MaxTick = NumTicks - NumTicks%5
#             TickStep = MaxTick/5
#             Ticks = list(np.linspace(0,TickStep*5,6))
#             ax3.set_yticks(Ticks)
#             ax3.set_yticklabels([""]*len(Ticks))
#         # ax3.set_yticklabels([str(int(el)) for el in \
#         #                         list(np.linspace(0,\
#         #                             np.ceil(max(u1[int(len(u1)/2):])*1.1) - \
#         #                                 np.ceil(max(u1[int(len(u1)/2):])*1.1)%3,4))],\
#         #                                     fontsize=12)
#         ax3.spines['right'].set_visible(False)
#         ax3.spines['top'].set_visible(False)
#         ax3.set_title("Force 1",fontsize=16,fontweight = 4,color = 'g',y = 0.95)
#         # ax3.set_xlabel("Time (s)")
#
#         #Force 2
#
#         Force2, = ax2.plot([0],[u2[0]],color = 'r')
#         ax2.set_xlim(0,Time[-1])
#         ax2.set_xticks(list(np.linspace(0,Time[-1],5)))
#         ax2.set_xticklabels([str(0),'','','',str(Time[-1])])
#         ax2.set_ylim(0,1.15*max(u1[5000:]+u2[5000:]))
#         ax2.set_yticks(list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                         int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1)))
#         ax2.set_yticklabels([str(int(el)) for el in \
#                                 list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                                     int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1))],\
#                                         fontsize=12)
#         if np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                         int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1).shape[0] < 5:
#             ax2.set_yticks(list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                             int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1)))
#             ax2.set_yticklabels([str(int(el)) for el in \
#                                     list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                                         int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1))],\
#                                             fontsize=12)
#         else:
#             NumTicks = np.floor(1.15*max(u1[5000:]+u2[5000:]))
#             MaxTick = NumTicks - NumTicks%5
#             TickStep = MaxTick/5
#             Ticks = list(np.linspace(0,TickStep*5,6))
#             ax2.set_yticks(Ticks)
#             ax2.set_yticklabels([str(tick) for tick in Ticks])
#         ax2.spines['right'].set_visible(False)
#         ax2.spines['top'].set_visible(False)
#         ax2.set_title("Force 2",fontsize=16,fontweight = 4,color = 'r',y = 0.95)
#         # ax2.set_xlabel("Time (s)")
#
#         # Trajectory
#
#         Predicted, = ax4.plot(Time,r(Time),'0.60',linestyle='--')
#         Actual, = ax4.plot([0],[x1[0]],'b')
#         ax4.set_xlim(0,Time[-1])
#         ax4.set_xticks(list(np.linspace(0,Time[-1],5)))
#         ax4.set_xticklabels([str(0),'','','',str(Time[-1])])
#         ax4.set_ylim([-1.25*A,1.25*A])
#         ax4.set_yticks([-A,0,A])
#         ax4.set_xlabel("Time (s)")
#         ax4.set_ylabel("Position of Center Mass (m)")
#         ax4.spines['right'].set_visible(False)
#         ax4.spines['top'].set_visible(False)
#
#         # Error
#         ErrorArray = x1-r(Time)
#         Error, = ax5.plot([0],[ErrorArray[0]],'k')
#         ax5.set_xlim(0,Time[-1])
#         ax5.set_xticks(list(np.linspace(0,Time[-1],5)))
#         ax5.set_xticklabels([str(0),'','','',str(Time[-1])])
#         ax5.set_ylim([ErrorArray.min() - 0.1*(max(ErrorArray)-min(ErrorArray)),\
#                         ErrorArray.max() + 0.1*(max(ErrorArray)-min(ErrorArray))])
#         ax5.set_xlabel("Time (s)")
#         ax5.set_ylabel("Error (m)")
#         ax5.yaxis.set_label_position("right")
#         ax5.yaxis.tick_right()
#         ax5.spines['left'].set_visible(False)
#         ax5.spines['top'].set_visible(False)
#
#         def animate(i):
#             Spring1.set_xdata(np.linspace(x1[i]+CenterBoxHalfWidth+StraightLength,\
#                                     RestingLength+x3[i]-SideBoxHalfWidth-StraightLength,1001))
#             Spring1_left.set_xdata([x1[i]+CenterBoxHalfWidth,x1[i]+CenterBoxHalfWidth+StraightLength])
#             Spring1_right.set_xdata([RestingLength+x3[i]-SideBoxHalfWidth-StraightLength,\
#                                         RestingLength+x3[i]-SideBoxHalfWidth])
#
#             Spring2.set_xdata(np.linspace(-RestingLength+x4[i]+SideBoxHalfWidth+StraightLength,\
#                                     x1[i]-CenterBoxHalfWidth-StraightLength,1001))
#             Spring2_left.set_xdata([x1[i]-CenterBoxHalfWidth-StraightLength,x1[i]-CenterBoxHalfWidth])
#             Spring2_right.set_xdata([-RestingLength+x4[i]+SideBoxHalfWidth,\
#                                         -RestingLength+x4[i]+SideBoxHalfWidth+StraightLength])
#
#             CenterMass.xy = (-CenterBoxHalfWidth + x1[i],-CenterBoxHalfHeight)
#             Mass1.xy = (-SideBoxHalfWidth+RestingLength + x3[i],-SideBoxHalfHeight)
#             Mass2.xy = (-SideBoxHalfWidth-RestingLength + x4[i],-SideBoxHalfHeight)
#             PositionArrow.set_xdata([x1[i],x1[i]])
#             PositionArrowHead.set_xdata([x1[i]])
#             PositionArrowTail.set_xdata([x1[i]])
#             Force1Arrow.set_xdata([RestingLength+x3[i]+(5/3)*SideBoxHalfWidth,\
#                                     RestingLength + x3[i]+(5/3)*SideBoxHalfWidth\
#                                         +ForceScaling*u1[i]/(max(u1[5000:]+u2[5000:]))])
#             Force1ArrowHead.set_xdata([RestingLength + x3[i]+(5/3)*SideBoxHalfWidth\
#                                         +ForceScaling*u1[i]/(max(u1[5000:]+u2[5000:]))])
#             Force2Arrow.set_xdata([x4[i]-RestingLength-(5/3)*SideBoxHalfWidth\
#                                         -ForceScaling*u2[i]/(max(u1[5000:]+u2[5000:])),\
#                                             x4[i]-RestingLength-(5/3)*SideBoxHalfWidth])
#             Force2ArrowHead.set_xdata([x4[i]-RestingLength-(5/3)*SideBoxHalfWidth\
#                                         -ForceScaling*u2[i]/(max(u1[5000:]+u2[5000:]))])
#
#             Force1.set_xdata(Time[:i])
#             Force1.set_ydata(u1[:i])
#
#             Force2.set_xdata(Time[:i])
#             Force2.set_ydata(u2[:i])
#
#             Actual.set_xdata(Time[:i])
#             Actual.set_ydata(x1[:i])
#
#             Error.set_xdata(Time[:i])
#             Error.set_ydata(ErrorArray[:i])
#
#             return Spring1,Spring1_left,Spring1_right,Spring2,Spring2_left,Spring2_right,CenterMass,Mass1,Mass2,Force1,Force2,Actual,Error,PositionArrow,PositionArrowHead,PositionArrowTail,Force1Arrow,Force1ArrowHead,Force2Arrow,Force2ArrowHead,
#
#         # Init only required for blitting to give a clean slate.
#         def init():
#             Spring1, =\
#                 ax1.plot(np.linspace(x1[0]+CenterBoxHalfWidth+StraightLength,\
#                                         RestingLength+x3[0]-SideBoxHalfWidth-StraightLength,1001),\
#                                             Spring_array,'k')
#             Spring1_left, = \
#                 ax1.plot([x1[0]+CenterBoxHalfWidth,x1[0]+CenterBoxHalfWidth+StraightLength],\
#                     [0,0],'k')
#             Spring1_right, = \
#                 ax1.plot([RestingLength+x3[0]-SideBoxHalfWidth-StraightLength,\
#                             RestingLength+x3[0]-SideBoxHalfWidth],[0,0],'k')
#             Spring2, =\
#                 ax1.plot(np.linspace(-RestingLength+x4[0]+SideBoxHalfWidth+StraightLength,\
#                                         x1[0]-CenterBoxHalfWidth-StraightLength,1001),\
#                                             Spring_array,'k')
#             Spring2_left, =\
#                 ax1.plot([x1[0]-CenterBoxHalfWidth-StraightLength,x1[0]-CenterBoxHalfWidth],\
#                             [0,0],'k')
#             Spring2_right, = \
#                 ax1.plot([-RestingLength+x4[0]+SideBoxHalfWidth,\
#                             -RestingLength+x4[0]+SideBoxHalfWidth+StraightLength],[0,0],'k')
#
#             CenterMass = plt.Rectangle((-CenterBoxHalfWidth,-CenterBoxHalfHeight),\
#                                         2*CenterBoxHalfWidth,2*CenterBoxHalfHeight,Color='#4682b4')
#             ax1.add_patch(CenterMass)
#             Mass1 = plt.Rectangle((-SideBoxHalfWidth+RestingLength,-SideBoxHalfHeight),\
#                                         2*SideBoxHalfWidth,2*SideBoxHalfHeight,Color='#4682b4')
#             ax1.add_patch(Mass1)
#             Mass2 = plt.Rectangle((-SideBoxHalfWidth-RestingLength,-SideBoxHalfHeight),\
#                                         2*SideBoxHalfWidth,2*SideBoxHalfHeight,Color='#4682b4')
#             ax1.add_patch(Mass2)
#
#             PositionArrow, = ax1.plot([x1[0],x1[0]],[0,2*CenterBoxHalfHeight],'k')
#             PositionArrowHead, = ax1.plot([x1[0]],[2*CenterBoxHalfHeight],'k^')
#             PositionArrowTail, = ax1.plot([x1[0]],[0],'ko')
#
#             Force1Arrow, = ax1.plot([RestingLength+x3[0]+(5/3)*SideBoxHalfWidth,\
#                                     RestingLength + x3[0]+(5/3)*SideBoxHalfWidth\
#                                         +ForceScaling*u1[0]/(max(u1[5000:]+u2[5000:]))],\
#                                             [0,0],'g')
#             Force1ArrowHead, = \
#                 ax1.plot([RestingLength + x3[0]+(5/3)*SideBoxHalfWidth\
#                         +ForceScaling*u1[0]/(max(u1[5000:]+u2[5000:]))],[0],'g<')
#             Force2Arrow, = ax1.plot([x4[0]-RestingLength-(5/3)*SideBoxHalfWidth\
#                         -ForceScaling*u2[0]/(max(u1[5000:]+u2[5000:])),\
#                             x4[0]-RestingLength-(5/3)*SideBoxHalfWidth],[0,0],'r')
#             Force2ArrowHead, = \
#                 ax1.plot([x4[0]-RestingLength-(5/3)*SideBoxHalfWidth\
#                     -ForceScaling*u2[0]/(max(u1[5000:]+u2[5000:]))],[0],'r>')
#
#             Force1, = ax3.plot([0],[u1[0]],color = 'g')
#             Force2, = ax2.plot([0],[u2[0]],color = 'r')
#             Predicted, = ax4.plot(Time,r(Time),'0.60',linestyle='--')
#             Actual, = ax4.plot([0],[x1[0]],'b')
#             Error, = ax5.plot([0],[ErrorArray[0]],'k')
#
#             Spring1.set_visible(False)
#             Spring1_left.set_visible(False)
#             Spring1_right.set_visible(False)
#             Spring2.set_visible(False)
#             Spring2_left.set_visible(False)
#             Spring2_right.set_visible(False)
#             CenterMass.set_visible(False)
#             Mass1.set_visible(False)
#             Mass2.set_visible(False)
#             PositionArrow.set_visible(False)
#             PositionArrowHead.set_visible(False)
#             PositionArrowTail.set_visible(False)
#             Force1.set_visible(False)
#             Force2.set_visible(False)
#             Predicted.set_visible(False)
#             Actual.set_visible(False)
#             Error.set_visible(False)
#             Force1Arrow.set_visible(False)
#             Force1ArrowHead.set_visible(False)
#             Force2Arrow.set_visible(False)
#             Force2ArrowHead.set_visible(False)
#
#             return Spring1,Spring1_left,Spring1_right,Spring2,Spring2_left,Spring2_right,CenterMass,Mass1,Mass2,Force1,Force2,Actual,Error,PositionArrow,PositionArrowHead,PositionArrowTail,Force1Arrow,Force1ArrowHead,Force2Arrow,Force2ArrowHead,
#
#         ani = animation.FuncAnimation(fig, animate, np.arange(1, len(Time),10), init_func=init,interval=1, blit=False)
#         # if save_as_gif:
#         # 	ani.save('test.gif', writer='imagemagick', fps=30)
#         plt.show()
# def plot_multiple_PDF_frames(response,Time,x1,x3,x4,u1,u2,FileName=None):
#     assert type(response)==bool, "Input must be either True or False."
#     if FileName != None: assert type(FileName)==str, "FileName must be a string"
#     if response == True:
#         import numpy as np
#         import matplotlib.pyplot as plt
#         from matplotlib.patches import Ellipse
#         import matplotlib.patches as patches
#         from scipy import signal
#         from matplotlib.backends.backend_pdf import PdfPages
#         import os.path
#
#         def return_fig(i):
#             fig = plt.figure(figsize=(10,8))
#             ax1 = plt.subplot2grid((3,4),(0,0),colspan=4)
#             ax2 = plt.subplot2grid((3,4),(1,0),colspan=2)
#             ax3 = plt.subplot2grid((3,4),(1,2),colspan=2)
#             ax4 = plt.subplot2grid((3,4),(2,0),colspan=3)
#             ax5 = plt.subplot2grid((3,4),(2,3))
#
#             plt.suptitle("Underdetermined Mass-Spring System",Fontsize=28,y=0.95)
#
#             # Model Drawing
#             IdealBoxScalingFactor = 0.78533496170320571 # Calculated from w = np.pi
#             CurrentTrialScalingFactor = max([max(x3)-min(x1),max(x1)-min(x4)])
#             StraightLength = 0.05*CurrentTrialScalingFactor/IdealBoxScalingFactor
#             RestingLength = max([max(x1)-min(x3),max(x4)-min(x1)])+2*StraightLength\
#                             +0.30*CurrentTrialScalingFactor/IdealBoxScalingFactor
#             CenterBoxHalfWidth = 0.15*CurrentTrialScalingFactor/IdealBoxScalingFactor
#             CenterBoxHalfHeight = 0.2*CurrentTrialScalingFactor/IdealBoxScalingFactor
#             SideBoxHalfWidth = 0.1*CurrentTrialScalingFactor/IdealBoxScalingFactor
#             SideBoxHalfHeight = 0.075*CurrentTrialScalingFactor/IdealBoxScalingFactor
#             ForceScaling = 1*CurrentTrialScalingFactor/IdealBoxScalingFactor
#
#             Spring_array =\
#              SideBoxHalfWidth\
#                 *np.abs(signal.sawtooth(5*2*np.pi*np.linspace(0,1,1001)-np.pi/2))\
#                     -(1/2)*SideBoxHalfWidth
#
#             Spring1, =\
#                 ax1.plot(np.linspace(x1[i]+CenterBoxHalfWidth+StraightLength,\
#                                         RestingLength+x3[i]-SideBoxHalfWidth-StraightLength,1001),\
#                                             Spring_array,'k')
#             Spring1_left, = \
#                 ax1.plot([x1[i]+CenterBoxHalfWidth,x1[i]+CenterBoxHalfWidth+StraightLength],\
#                             [0,0],'k')
#             Spring1_right, = \
#                 ax1.plot([RestingLength+x3[i]-SideBoxHalfWidth-StraightLength,\
#                             RestingLength+x3[i]-SideBoxHalfWidth],\
#                                 [0,0],'k')
#
#             Spring2, =\
#                 ax1.plot(np.linspace(-RestingLength+x4[i]+SideBoxHalfWidth+StraightLength,\
#                                         x1[i]-CenterBoxHalfWidth-StraightLength,1001),\
#                                             Spring_array,'k')
#             Spring2_left, = \
#                 ax1.plot([x1[i]-CenterBoxHalfWidth-StraightLength,x1[i]-CenterBoxHalfWidth],\
#                             [0,0],'k')
#             Spring2_right, = \
#                 ax1.plot([-RestingLength+x4[i]+SideBoxHalfWidth,\
#                             -RestingLength+x4[i]+SideBoxHalfWidth+StraightLength],\
#                                 [0,0],'k')
#             ax1.get_xaxis().set_ticks([])
#             ax1.get_yaxis().set_ticks([])
#             ax1.set_frame_on(True)
#             CenterMass = plt.Rectangle((-CenterBoxHalfWidth + x1[i],-CenterBoxHalfHeight),\
#                                         2*CenterBoxHalfWidth,2*CenterBoxHalfHeight,Color='#4682b4')
#             ax1.add_patch(CenterMass)
#             Mass1 = plt.Rectangle((-SideBoxHalfWidth+RestingLength + x3[i],-SideBoxHalfHeight),\
#                                         2*SideBoxHalfWidth,2*SideBoxHalfHeight,Color='#4682b4')
#             ax1.add_patch(Mass1)
#             Mass2 = plt.Rectangle((-SideBoxHalfWidth-RestingLength + x4[i],-SideBoxHalfHeight),\
#                                         2*SideBoxHalfWidth,2*SideBoxHalfHeight,Color='#4682b4')
#             ax1.add_patch(Mass2)
#
#             PositionArrow, = ax1.plot([x1[i],x1[i]],[0,2*CenterBoxHalfHeight],'k')
#             PositionArrowHead, = ax1.plot([x1[i]],[2*CenterBoxHalfHeight],'k^')
#             PositionArrowTail, = ax1.plot([x1[i]],[0],'ko')
#
#             Scale = ax1.plot([-1.1*A,1.1*A],\
#                                 [2.75*CenterBoxHalfHeight,2.75*CenterBoxHalfHeight],\
#                                     '0.60')
#             Ticks = np.linspace(-A,A,5)
#             TickHeights = [0.3*CenterBoxHalfHeight,\
#                             0.15*CenterBoxHalfHeight,\
#                             0.3*CenterBoxHalfHeight,\
#                             0.15*CenterBoxHalfHeight,\
#                             0.3*CenterBoxHalfHeight]
#             [ax1.plot([Ticks[i],Ticks[i]],\
#                     [2.75*CenterBoxHalfHeight-TickHeights[i],2.75*CenterBoxHalfHeight],'0.60') \
#                         for i in range(5)]
#
#             Force1Arrow, = ax1.plot([RestingLength+x3[i]+(5/3)*SideBoxHalfWidth,\
#                                         RestingLength + x3[i]+(5/3)*SideBoxHalfWidth\
#                                             +ForceScaling*u1[i]/(max(u1[5000:]+u2[5000:]))],\
#                                                     [0,0],'g')
#             Force1ArrowHead, = \
#                 ax1.plot([RestingLength + x3[i]+(5/3)*SideBoxHalfWidth\
#                             +ForceScaling*u1[i]/(max(u1[5000:]+u2[5000:]))],[0],'g>')
#             Force2Arrow, =\
#                 ax1.plot([x4[i]-RestingLength-(5/3)*SideBoxHalfWidth\
#                             -ForceScaling*u2[i]/(max(u1[5000:]+u2[5000:])),\
#                                 x4[i]-RestingLength-(5/3)*SideBoxHalfWidth],[0,0],'r')
#             Force2ArrowHead, = \
#                 ax1.plot([x4[i]-RestingLength-(5/3)*SideBoxHalfWidth\
#                             -ForceScaling*u2[i]/(max(u1[5000:]+u2[5000:]))],[0],'r<')
#
#             LowerBound = (np.array(x4[5001:])-RestingLength-(5/3)*SideBoxHalfWidth\
#                             -ForceScaling*np.array(u2[5000:])/(max(u1[5000:]+u2[5000:]))).min()
#             UpperBound = (RestingLength + np.array(x3[5001:])+(5/3)*SideBoxHalfWidth\
#                             +ForceScaling*np.array(u1[5000:])/(max(u1[5000:]+u2[5000:]))).max()
#             Bound = 1.05*np.array([-LowerBound,UpperBound]).max()
#             ax1.set_xlim([-Bound,Bound])
#             ax1.set_ylim([-1.5*CenterBoxHalfHeight,3.25*CenterBoxHalfHeight])
#             ax1.set_aspect('equal')
#
#             #Force 1
#
#             Force1, = ax3.plot(Time[:i],u1[:i],color = 'g')
#             ax3.set_xlim(0,Time[-1])
#             ax3.set_xticks(list(np.linspace(0,Time[-1],5)))
#             ax3.set_xticklabels([str(0),'','','',str(Time[-1])])
#             ax3.set_ylim(0,1.15*max(u1[5000:]+u2[5000:]))
#             if np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                             int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1).shape[0] < 5:
#                 ax3.set_yticks(list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                                 int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1)))
#                 ax3.set_yticklabels([""]*(int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1))
#             else:
#                 NumTicks = np.floor(1.15*max(u1[5000:]+u2[5000:]))
#                 MaxTick = NumTicks - NumTicks%5
#                 TickStep = MaxTick/5
#                 Ticks = list(np.linspace(0,TickStep*5,6))
#                 ax3.set_yticks(Ticks)
#                 ax3.set_yticklabels([""]*len(Ticks))
#             # ax3.set_yticklabels([str(int(el)) for el in \
#             #                         list(np.linspace(0,\
#             #                             np.ceil(max(u1[int(len(u1)/2):])*1.1) - \
#             #                                 np.ceil(max(u1[int(len(u1)/2):])*1.1)%3,4))],\
#             #                                     fontsize=12)
#             ax3.spines['right'].set_visible(False)
#             ax3.spines['top'].set_visible(False)
#             ax3.set_title("Force 1",fontsize=16,fontweight = 4,color = 'g',y = 0.95)
#             # ax3.set_xlabel("Time (s)")
#
#             #Force 2
#
#             Force2, = ax2.plot(Time[:i],u2[:i],color = 'r')
#             ax2.set_xlim(0,Time[-1])
#             ax2.set_xticks(list(np.linspace(0,Time[-1],5)))
#             ax2.set_xticklabels([str(0),'','','',str(Time[-1])])
#             ax2.set_ylim(0,1.15*max(u1[5000:]+u2[5000:]))
#             ax2.set_yticks(list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                             int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1)))
#             ax2.set_yticklabels([str(int(el)) for el in \
#                                     list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                                         int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1))],\
#                                             fontsize=12)
#             if np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                             int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1).shape[0] < 5:
#                 ax2.set_yticks(list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                                 int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1)))
#                 ax2.set_yticklabels([str(int(el)) for el in \
#                                         list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
#                                             int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1))],\
#                                                 fontsize=12)
#             else:
#                 NumTicks = np.floor(1.15*max(u1[5000:]+u2[5000:]))
#                 MaxTick = NumTicks - NumTicks%5
#                 TickStep = MaxTick/5
#                 Ticks = list(np.linspace(0,TickStep*5,6))
#                 ax2.set_yticks(Ticks)
#                 ax2.set_yticklabels([str(tick) for tick in Ticks])
#             ax2.spines['right'].set_visible(False)
#             ax2.spines['top'].set_visible(False)
#             ax2.set_title("Force 2",fontsize=16,fontweight = 4,color = 'r',y = 0.95)
#             # ax2.set_xlabel("Time (s)")
#
#             # Trajectory
#
#             Predicted, = ax4.plot(Time,r(Time),'0.60',linestyle='--')
#             Actual, = ax4.plot(Time[:i],x1[:i],'b')
#             ax4.set_xlim(0,Time[-1])
#             ax4.set_xticks(list(np.linspace(0,Time[-1],5)))
#             ax4.set_xticklabels([str(0),'','','',str(Time[-1])])
#             ax4.set_ylim([-1.25*A,1.25*A])
#             ax4.set_yticks([-A,0,A])
#             ax4.set_xlabel("Time (s)")
#             ax4.set_ylabel("Position of Center Mass (m)")
#             ax4.spines['right'].set_visible(False)
#             ax4.spines['top'].set_visible(False)
#
#             # Error
#             ErrorArray = x1-r(Time)
#             Error, = ax5.plot(Time[:i],ErrorArray[:i],'k')
#             ax5.set_xlim(0,Time[-1])
#             ax5.set_xticks(list(np.linspace(0,Time[-1],5)))
#             ax5.set_xticklabels([str(0),'','','',str(Time[-1])])
#             ax5.set_ylim([ErrorArray.min() - 0.1*(max(ErrorArray)-min(ErrorArray)),\
#                             ErrorArray.max() + 0.1*(max(ErrorArray)-min(ErrorArray))])
#             ax5.set_xlabel("Time (s)")
#             ax5.set_ylabel("Error (m)")
#             ax5.yaxis.set_label_position("right")
#             ax5.yaxis.tick_right()
#             ax5.spines['left'].set_visible(False)
#             ax5.spines['top'].set_visible(False)
#
#             return(fig)
#         i = 1
#         if FileName == None:
#             FileName = "ReferenceTrackingTest.pdf"
#         else:
#             FileName = FileName + ".pdf"
#         if os.path.exists(FileName) == True:
#             while os.path.exists(FileName) == True:
#                 i += 1
#                 FileName = FileName[:-4]
#                 FileName = FileName	+ "_" + "{:0>2d}".format(i) +".pdf"
#         PDFFile = PdfPages(FileName)
#
#         for i in np.linspace(0,len(Time)-1,201)[:-1]:
#             t_i = int(i)
#             fig = return_fig(t_i)
#             PDFFile.savefig(fig)
#             plt.close("all")
#         PDFFile.close()
#
# N = 10001
# Time = np.linspace(0,10,N)
# x1,x2,x3,x4,x5,x6 = [A],[-1],[0],[0],[0],[0]
# u1,u2 = [],[]
#
# AddNoise = True
# if AddNoise == True:
#     np.random.seed(seed=None)
#     NoiseArray = np.random.normal(loc=0.0,scale=0.2,size=(2,len(Time)))
# else:
#     NoiseArray = np.zeros((2,len(Time)))
#
# def update_policy(t,x1,x2,x3,x4,x5,x6,dt,NoiseArray,e=2):
#     import numpy as np
#     X = [x1[-1],x2[-1],x3[-1],x4[-1],x5[-1],x6[-1]]
#     U = return_U(t,X,e,NoiseArray[:,int(t/dt)])
#     u1.append(U[0])
#     u2.append(U[1])
#     x6.append(x6[-1] + dx6(t,X,U)*dt)
#     x5.append(x5[-1] + dx5(t,X,U)*dt)
#     x4.append(x4[-1] + dx4(t,X)*dt)
#     x3.append(x3[-1] + dx3(t,X)*dt)
#     x2.append(x2[-1] + dx2(t,X)*dt)
#     x1.append(x1[-1] + dx1(t,X)*dt)
#
# for t in Time[1:]:
#     update_policy(t,x1,x2,x3,x4,x5,x6,dt,NoiseArray,e=CocontractionIndex)
#
# plt.figure()
# # plt.title(r'$\dot{x}_{1} = x_{1}^{2} + x_{2}; \hspace{1em} \dot{x}_{2} = u$',\
# #                 fontsize=16,color='gray')
# plt.title("Underdetermined Spring Example",\
#                 fontsize=16,color='gray')
# plt.plot(Time,x1,'b',lw=2)
# plt.plot(Time,r(Time),'r--')
# plt.xlabel("Time (s)")
# plt.ylabel("Desired Measure")
# plt.legend([r"Output $y = x_{1}$",r"Reference $r(t) = " + str(A) + "\sin(" + str(w) + "t)$"])
#
# plt.figure()
# plt.title('Error vs. Time')
# plt.plot(Time, r(Time)-x1,color='r')
# plt.xlabel("Time (s)")
# plt.ylabel("Error")
# #
# # plt.show()
#
# plt.close('all')
# plot_multiple_PDF_frames(False,Time,x1,x3,x4,u1,u2)
# animate_trajectory(True,Time,x1,x3,x4,u1,u2)
