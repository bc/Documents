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
	if Title != '': Title = ' '*(22-len(Title)) + Title + ': '
	statusbar = Title +'[' + '\u25a0'*int((i+1)/(N/50)) + '\u25a1'*(50-int((i+1)/(N/50))) + '] '
	if StartTime != False:
		if i==0:
			time_array = []
			TimeLeft = '--'
		elif i==int(0.02*N):
			time_array.append(time.time()-StartTime)
			TimeLeft = '{0:1.1f}'.format(time_array[-1]*(N/(i+1)))
		elif i%int(0.02*N)==0:
			time_array.append(time.time()-StartTime)
			TimeLeft = '{0:1.1f}'.format(float(interpolate.interp1d(np.arange(len(time_array)),time_array,fill_value='extrapolate')(49)))
		print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete, ' + '{0:1.1f}'.format(time.time() - StartTime) \
			+ 'sec, (est. ' + TimeLeft,' sec left)		\r', end='')
	else:
		print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete           \r',end = '')
def set_link_lengths(New_L1=0.256,New_L2=0.315,EOM = "Uno"):
    """
    Sets the link lengths for a 2 DOF planar reaching task. Changes the values of the link lengths. New_L1 and New_L2 must be numbers. Set EOM to Uno of Zadravec.
    """
    assert type(New_L1) == int or type(New_L1)==float, "New_L1 must be a number."
    assert type(New_L2) == int or type(New_L2)==float, "New_L2 must be a number."
    assert EOM in ["Uno","Zadravec"], "EOM can be either 'Uno' or 'Zadravec'"
    global L1, L2
    if New_L1 != 0.256 or New_L2 != 0.315:
        L1 = New_L1
        L2 = New_L2
    elif EOM == "Uno":
        L1 = 0.256
        L2 = 0.315
    else: # EOM = "Zadravec"
        L1 = 0.298
        L2 = 0.419
def test_global_link_lengths():
    print('L1 = ' + str(L1))
    print('L2 = ' + str(L2))
def create_angle_lists():
	"""
	Creates global lists of angles, angular velocities, and angular accelerations.
	"""
	global A1,A2,Ȧ1,Ȧ2,Ä1,Ä2
	A1,A2=[],[]
	Ȧ1,Ȧ2=[],[]
	Ä1,Ä2=[],[]
def inverse_kinematics(X):
    """
    Takes in either a (2,) or (2,1) list/array with values for x and y (endpoint) and updates the lists A1 and A2 with the current angles (in radians) from the inverse kinematics.
    """
    import numpy as np
    from math import acos, sin, cos, atan, atan2
    from numpy import pi
    assert len(X)==2, "X must be either a (2,) or (2,1) list/array"
    x,y = X
    if np.shape(X) == (2,1): x,y = x[0],y[0]
    # x² + y² > (L₁+L₂)² = L₁² + 2L₁L₂ + L₂² > L₁² + L₂²
    # x² + y² - L₁² - L₂² > 0
    # a₂ = cos⁻¹((x² + y² - L₁² - L₂²)/2L₁L₂) ∊ (0,π/2)
    a2 = acos((x**2 + y**2 - L1**2 - L2**2)/(2*L1*L2))
    a1 = atan2(y,x) - atan2(L2*sin(a2),(L1+L2*cos(a2)))
    global A1,A2
    A1.append(a1)
    A2.append(a2)
def find_X_values(t,Xi,Xf):
    """
    This takes t (a numpy.ndarray of normalized time ∈ [0,1]) and either a (2,) or (2,1) list/array with values for both initial and final x and y (endpoint). To avoid singularities, ||X[i]|| cannot be greater than L1 + L2.
    """
    import numpy as np
    assert len(Xf)==2, "Xf must be either a (2,) or (2,1) list/array"
    assert len(Xi)==2, "Xi must be either a (2,) or (2,1) list/array"
    xi,yi = Xi
    xf,yf = Xf
    if np.shape(Xi) == (2,1): xi,yi = xi[0],yi[0]
    if np.shape(Xf) == (2,1): xf,yf = xf[0],yf[0]
    x = xi + (xf-xi)*(10*t**3 - 15*t**4 + 6*t**5)
    y = yi + (yf-yi)*(10*t**3 - 15*t**4 + 6*t**5)
    X = np.array(np.concatenate(([x],[y]),axis=0))
    assert (sum(X**2)**0.5<L1+L2).all(), "Trajectory creates singularities."
    return(X)
def find_Ẋ_values(t,Xi,Xf,t_end = 1):
    """
    This takes t (a numpy.ndarray of normalized time ∈ [0,1]) and either a (2,) or (2,1) list/array with values for both initial and final x and y (endpoint). To avoid singularities, ||X[i]|| cannot be greater than L1 + L2.
    """
    import numpy as np
    assert len(Xi)==2, "Xi must be either a (2,) or (2,1) list/array"
    assert len(Xf)==2, "Xf must be either a (2,) or (2,1) list/array"
    xi,yi = Xi
    xf,yf = Xf
    if np.shape(Xi) == (2,1): xi,yi = xi[0],yi[0]
    if np.shape(Xf) == (2,1): xf,yf = xf[0],yf[0]
    ẋ = (xf-xi)*(30*t**2 - 60*t**3 + 30*t**4)/t_end
    ẏ = (yf-yi)*(30*t**2 - 60*t**3 + 30*t**4)/t_end
    Ẋ = np.array(np.concatenate(([ẋ],[ẏ]),axis=0))
    return(Ẋ)
def find_Ẍ_values(t,Xi,Xf,t_end=1):
    """
    This takes t (a numpy.ndarray of normalized time ∈ [0,1]) and either a (2,) or (2,1) list/array with values for both initial and final x and y (endpoint).
    """
    import numpy as np
    assert len(Xi)==2, "Xi must be either a (2,) or (2,1) list/array"
    assert len(Xf)==2, "Xf must be either a (2,) or (2,1) list/array"
    xi,yi = Xi
    xf,yf = Xf
    if np.shape(Xi) == (2,1): xi,yi = xi[0],yi[0]
    if np.shape(Xf) == (2,1): xf,yf = xf[0],yf[0]
    ẍ = (xf-xi)*(60*t - 180*t**2 + 120*t**3)/(t_end**2)
    ÿ = (yf-yi)*(60*t - 180*t**2 + 120*t**3)/(t_end**2)
    Ẍ = np.array(np.concatenate(([ẍ],[ÿ]),axis=0))
    return(Ẍ)
def return_X_values(t,Xi,Xf,t_end):
    X = find_X_values(t,Xi,Xf)
    Ẋ = find_Ẋ_values(t,Xi,Xf,t_end=t_end)
    Ẍ = find_Ẍ_values(t,Xi,Xf,t_end=t_end)
    return(X,Ẋ,Ẍ)
def jacobian():
	import numpy as np
	from math import cos, sin
	J = np.matrix([[-L1*sin(A1[-1])-L2*sin(A1[-1]+A2[-1]),-L2*sin(A1[-1]+A2[-1])],\
				[L1*cos(A1[-1])+L2*cos(A1[-1]+A2[-1]),L2*cos(A1[-1]+A2[-1])]])
	return(J)
def update_angular_velocity(Ẋ):
    import numpy as np
    J = jacobian()
    Ȧ = J**(-1)*np.array([Ẋ]).T
    global Ȧ1,Ȧ2
    Ȧ1.append(Ȧ[0,0])
    Ȧ2.append(Ȧ[1,0])
def update_angular_acceleration(Ẋ,Ẍ):
    from math import cos,sin
    import numpy as np
    assert len(Ẋ)==2, "Ẋ must be either a (2,) or (2,1) list/array"
    assert len(Ẍ)==2, "Ẍ must be either a (2,) or (2,1) list/array"
    ẋ,ẏ = Ẋ
    ẍ,ÿ = Ẍ
    if np.shape(Ẋ) == (2,1): ẋ,ẏ = ẋ[0],ẏ[0]
    if np.shape(Ẍ) == (2,1): ẍ,ÿ = ẍ[0],ÿ[0]
    a1,a2 = A1[-1],A2[-1]
    ȧ1,ȧ2 = Ȧ1[-1],Ȧ2[-1]
    ä1 = ( ( (ȧ1+ȧ2)*(ẏ*cos(a1+a2)-ẋ*sin(a1+a2)) \
            + ẍ*cos(a1+a2) \
            + ÿ*sin(a1+a2) ) * L1*sin(a2) \
            - \
            (ẋ*cos(a1+a2) + ẏ*sin(a1+a2)) * L1*cos(a2)*ȧ2 ) \
            / ((L1*sin(a2))**2)
    ä2 = ( ( \
            ((ȧ1+ȧ2)*L2*sin(a1+a2) + ȧ1*L1*sin(a1))*ẋ \
            + (-L2*cos(a1+a2) - L1*cos(a1))*ẍ \
            + (-(ȧ1+ȧ2)*L2*cos(a1+a2) - ȧ1*L1*cos(a1))*ẏ \
            + (-L2*sin(a1+a2)-L1*sin(a1))*ÿ ) \
            * (L1*L2*sin(a2)) \
            - \
            ( (-L2*cos(a1+a2) - L1*cos(a1))*ẋ \
            + (-L2*sin(a1+a2) - L1*sin(a1))*ẏ ) * (L1*L2*cos(a2)*ȧ2) \
            ) \
            / (L1*L2*sin(a2)**2)
    global Ä1,Ä2
    Ä1.append(ä1)
    Ä2.append(ä2)
def update_angle_lists(X,Ẋ,Ẍ):
    import numpy as np
    for i in range(np.shape(X)[1]):
        inverse_kinematics(X[:,i])
        update_angular_velocity(Ẋ[:,i])
        update_angular_acceleration(Ẋ[:,i],Ẍ[:,i])
def reaching_task(dt = 0.001, Xi = [0,0.25], Xf = [0.25,0.50],t_end=1):
    import numpy as np
    set_link_lengths()
    create_angle_lists()
    t = np.arange(0,1+dt,dt)
    X,Ẋ,Ẍ = return_X_values(t,Xi,Xf,t_end)
    update_angle_lists(X,Ẋ,Ẍ)
def plot_resulting_kinematics():
    import matplotlib.pyplot as plt
    import numpy as np
    from math import cos,sin

    dt = 1/(len(A1)-1)
    t = np.arange(0,1+dt,dt)
    x = np.array([L1*cos(A1[i])+L2*cos(A1[i]+A2[i]) for i in range(len(A1))])
    y = np.array([L1*sin(A1[i])+L2*sin(A1[i]+A2[i]) for i in range(len(A1))])

    plt.figure()
    plt.plot(x,y)
    ax1 = plt.gca()
    ax1.set_aspect('equal')
    ax1.set_title("Total Trajectory")
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    ax1.set_ylim(ymin=min(y)-abs(0.05*(max(y)-min(y))),ymax=max(y)+abs(0.05*(max(y)-min(y))))
    ax1.set_xlim(xmin=min(x)-abs(0.05*(max(x)-min(x))),xmax=max(x)+abs(0.05*(max(x)-min(x))))

    plt.figure()
    plt.plot(t,x)
    ax2 = plt.gca()
    ax2.set_title("Trajectory x-Component")
    ax2.set_xlabel('Normalized Time')
    ax2.set_ylabel('x (m)')

    plt.figure()
    plt.plot(t,y)
    ax3 = plt.gca()
    ax3.set_title("Trajectory y-Component")
    ax3.set_xlabel('Normalized Time')
    ax3.set_ylabel('y (m)')

    # plt.show()
def MA_function(Coefficients,Angles = None):
	"""
	Note: Angles should be a number if Coefficients has a length of 5, or a list of length 2 when the Coefficients have lengths 16 or 18. Angles[0] will be the PRIMARY ANGLE for the DOF being considered while Angles[1] will be the secondary angle.
	"""
	import numpy as np
<<<<<<< HEAD
	assert type(Coefficients) == list, "Coefficients must be a list."
	assert len(Coefficients) in [5,16,18], "Coefficients must be a list of length 5, 16, or 18."
	if len(Coefficients)==5:
		MomentArm = (np.matrix(Coefficients)*np.matrix([1,Angles,Angles**2,Angles**3,Angles**4]).T)[0,0]
=======
	import sympy as sp
	global q1,q2,q_alt

	assert type(Coefficients) == list or type(Coefficients) == int or type(Coefficients) == float, "Coefficients must be a list, int, or float."
	if type(Coefficients) not in [int,float]:
		assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."

	if type(Coefficients) in [int,float]:
		q2 = sp.symbols('q2')
		MomentArm = Coefficients
	elif len(Coefficients)==5:
		q2 = sp.symbols('q2')
		MomentArm = (sp.Matrix(Coefficients).T*sp.Matrix([1,q2,q2**2,q2**3,q2**4]))[0,0]
>>>>>>> 7d8b48796ec1a090f503dce659a61ec7087f6ba7
	elif len(Coefficients)==16:
		q2,q_alt = sp.symbols('q2'),sp.symbols('q_alt')
		MomentArm = (sp.Matrix(Coefficients).T*\
						sp.Matrix([1, q2, q_alt, q2*q_alt, q2**2, \
									q_alt**2, (q2**2)*q_alt, q2*(q_alt**2), \
									(q2**2)*(q_alt**2), q2**3, q_alt**3, \
									(q2**3)*q_alt, q2*(q_alt**3), \
									(q2**3)*(q_alt**2), (q2**2)*(q_alt**3), \
									(q2**3)*(q_alt**3)]))[0, 0]
	else: # len(Coefficients)==18
		q2,q_alt = sp.symbols('q2'),sp.symbols('q_alt')
		MomentArm = (sp.Matrix(Coefficients).T*\
						sp.Matrix([1, q2, q_alt, q2*q_alt, q2**2, \
							q_alt**2, (q2**2)*q_alt, q2*(q_alt**2), (q2**2)*(q_alt**2), \
							q2**3, (q2**3)*q_alt, (q2**3)*(q_alt**2), \
							q2**4, (q2**4)*q_alt, (q2**4)*(q_alt**2),  \
							q2**5, (q2**5)*q_alt, (q2**5)*(q_alt**2)]))[0, 0]
	if Angles == None:
		return(MomentArm)
	else:
		if type(Coefficients) in [int,float]:
			return(MomentArm.subs(q2,Angles))
		elif len(Coefficients)==5:
			return(MomentArm.subs(q2,Angles))
		elif len(Coefficients)==16:
			return(MomentArm.subs([(q2,Angles[0]),(q_alt,Angles[1])]))
		else: # len(Coefficients)==18
			return(MomentArm.subs([(q2,Angles[0]),(q_alt,Angles[1])]))
def MA_function_integral(Coefficients,Angles = None):
	"""
	Note: Angles should be a number if Coefficients has a length of 1 or 5, or a list of length 2 when the Coefficients have lengths 16 or 18. a1 will be the PRIMARY ANGLE for the DOF being considered while Angles[1] will be the secondary angle.
	"""
	import numpy as np
<<<<<<< HEAD
=======
	import sympy as sp
>>>>>>> 7d8b48796ec1a090f503dce659a61ec7087f6ba7
	import ipdb
	assert type(Coefficients) == list or type(Coefficients) == int or type(Coefficients) == float, "Coefficients must be a list, int, or float."
	if type(Coefficients) not in [int,float]:
		assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
	if type(Coefficients) in [int,float]:
<<<<<<< HEAD
		MomentArmIntegral = Coefficients*Angles
	elif len(Coefficients)==5:
		integral = np.multiply(np.matrix([1,Angles,Angles**2,Angles**3,Angles**4]),(1,1/2,1/3,1/4,1/5))*Angles
=======
		q2 = sp.symbols('q2')
		MomentArmIntegral = Coefficients*q2
	elif len(Coefficients)==5:
		q2 = sp.symbols('q2')
		integral = np.multiply(sp.Matrix([1,q2,q2**2,q2**3,q2**4]).T,(1,1/2,1/3,1/4,1/5))*q2
>>>>>>> 7d8b48796ec1a090f503dce659a61ec7087f6ba7
		MomentArmIntegral = (np.matrix(Coefficients)*integral.T)[0,0]
	elif len(Coefficients)==16:
		q2,q_alt = sp.symbols('q2'),sp.symbols('q_alt')
		integral = np.multiply(sp.Matrix([1, q2, q_alt, q2*q_alt, q2**2, \
						q_alt**2, (q2**2)*q_alt, q2*(q_alt**2), \
						(q2**2)*(q_alt**2), q2**3, q_alt**3, \
						(q2**3)*q_alt, q2*(q_alt**3), \
						(q2**3)*(q_alt**2), (q2**2)*(q_alt**3), \
						(q2**3)*(q_alt**3)]).T,\
						[1,1/2,1,1/2,1/3,1,1/3,1/2,1/3,1/4,1,1/4,1/2,1/4,1/3,1/4])*q2
		MomentArmIntegral = (np.matrix(Coefficients)*integral.T)[0, 0]
	else: # len(Coefficients)==18
		q2,q_alt = sp.symbols('q2'),sp.symbols('q_alt')
		integral = np.multiply(sp.Matrix([1, q2, q_alt, q2*q_alt, q2**2, \
							q_alt**2, (q2**2)*q_alt, q2*(q_alt**2), (q2**2)*(q_alt**2), \
							q2**3, (q2**3)*q_alt, (q2**3)*(q_alt**2), \
							q2**4, (q2**4)*q_alt, (q2**4)*(q_alt**2),  \
							q2**5, (q2**5)*q_alt, (q2**5)*(q_alt**2)]).T,\
							[1,1/2,1,1/2,1/3,1,1/3,1/2,1/3,1/4,1/4,1/4,1/5,1/5,1/5,1/6,1/6,1/6])*q2
		MomentArmIntegral = (np.matrix(Coefficients)*integral.T)[0, 0]
	if Angles == None:
		return(MomentArmIntegral)
	else:
		if type(Coefficients) in [int,float]:
			return(MomentArmIntegral.subs(q2,Angles))
		elif len(Coefficients)==5:
			return(MomentArmIntegral.subs(q2,Angles))
		elif len(Coefficients)==16:
			return(MomentArmIntegral.subs([(q2,Angles[0]),(q_alt,Angles[1])]))
		else: # len(Coefficients)==18
			return(MomentArmIntegral.subs([(q2,Angles[0]),(q_alt,Angles[1])]))
def global_R_matrix():
	import sympy as sp
	import numpy as np
	from numpy import pi
	global q2,q_alt
	q2,q_alt = sp.symbols('q2'),sp.symbols('q_alt')
	DELTa_Coefficients = [19, 0]
	CB_Coefficients = [20, 0]
	DELTp_Coefficients = [-8, 0]
	BIC_Coefficients = [15, [8.4533,36.6147,2.4777,-19.432,2.0571,0,13.6502,0,0,-5.6172,0,-2.0854,0,0,0,0]]
	TRI_Coefficients = [-15, [-24.5454,-8.8691,9.3509,-1.7518,0]]
	BRA_Coefficients = [0, [16.1991,-16.1463,24.5512,-6.3335,0]]
	BRD_Coefficients = [0, [15.2564,-11.8355,2.8129,-5.7781,44.8143,0,2.9032,0,0,-13.4956,0,-0.3940,0,0,0,0]]
	PRO_Coefficients = [0, [11.0405,-1.0079,0.3933,-10.4824,-12.1639,-0.4369,36.9174,3.5232,-10.4223,21.2604,-37.2444,10.2666,-11.0060,14.5974,-3.9919,1.7526,-2.0089,0.5460]]
	FCR_Coefficients = [0, 14]
	ECRB_Coefficients = [0, [-11.256,17.8548,1.6398,-0.5073,-2.8827,0,-0.0942,0,0,0,0,0,0,0,0,0]]
	ECRL_Coefficients = [0, [-7.7034,16.3913,7.4361,-1.7566,0,-1.3336,0,0.0742,0,0,0,0,0,0,0,0]]
	FCU_Coefficients = [0, 19]
	FDS_Coefficients = [0, 20]
	PL_Coefficients = [0, 25]
	ECU_Coefficients = [0, -23]
	EDM_Coefficients = [0, -10]
	EDC_Coefficients = [0, -20]

	global AllCoefficients
	AllCoefficients = [DELTa_Coefficients, CB_Coefficients, DELTp_Coefficients, BIC_Coefficients, \
						TRI_Coefficients, BRA_Coefficients, BRD_Coefficients, PRO_Coefficients, \
						FCR_Coefficients, ECRB_Coefficients, ECRL_Coefficients, FCU_Coefficients, \
						FDS_Coefficients, PL_Coefficients, ECU_Coefficients, EDM_Coefficients, EDC_Coefficients]

	global RMatrix_Transpose
	RMatrix_Transpose = sp.Matrix([[MA_function(AllCoefficients[i][j]) for j in range(2)] for i in range(17)])
def return_MA_matrix(A1,A2):
	"""
	Notes:

	The angle of pronation is fixed at pi for this reaching paradigm.
	The angle of radial/ulnar deviation is set to zero for the fixed wrist apparatus of the model paradigm.
	These functions have been verified to match the previous posture dependent MAs from Ramsey (2010) - SEE ERRATUM
	"""
	import numpy as np
	from numpy import pi
<<<<<<< HEAD
	# MAs taken from Ramsay (2009), Holzbaur, FVC, Rankin (2012)
	#R_tranpose Column 1
	r1DELTa = 19 # in mm
	r1CB = 20 # in mm
	r1DELTp = -8 # in mm
	r1BIC = 15 # in mm
	r1TRI = -15 # in mm
	r1BRA = 0 # in mm
	r1BRD = 0 # in mm
	r1PRO = 0 # in mm
	r1FCR = 0 # in mm
	r1ECRB = 0 # in mm
	r1ECRL = 0 # in mm
	r1FCU = 0 # in mm
	r1FDS = 0 # in mm
	r1PL = 0 # in mm
	r1ECU = 0 # in mm
	r1EDM = 0 # in mm
	r1EDC = 0 # in mm

	# Note: Ramsay (2009) suggests these muscles have contributions to WFE
	#R_tranpose Column 2
	r2DELTa = 0 # in mm
	r2CB = 0 # in mm
	r2DELTp = 0 # in mm
	r2BIC = MA_function([8.4533,36.6147,2.4777,-19.432,2.0571,0,13.6502,0,0,-5.6172,0,-2.0854,0,0,0,0],[A2,pi]) # in mm
	r2TRI = MA_function([-24.5454,-8.8691,9.3509,-1.7518,0],A2) # in mm
	r2BRA = MA_function([16.1991,-16.1463,24.5512,-6.3335,0],A2) # in mm
	r2BRD = MA_function([15.2564,-11.8355,2.8129,-5.7781,44.8143,0,2.9032,0,0,-13.4956,0,-0.3940,0,0,0,0],[A2,pi]) # in mm
	r2PRO = MA_function([11.0405,-1.0079,0.3933,-10.4824,-12.1639,-0.4369,36.9174,3.5232,-10.4223,21.2604,-37.2444,10.2666,-11.0060,14.5974,-3.9919,1.7526,-2.0089,0.5460],[A2,pi]) # in mm
	r2FCR = 20 # in mm
	r2ECRB = MA_function([-11.256,17.8548,1.6398,-0.5073,-2.8827,0,-0.0942,0,0,0,0,0,0,0,0,0],[A2,pi]) # in mm
	r2ECRL = MA_function([-7.7034,16.3913,7.4361,-1.7566,0,-1.3336,0,0.0742,0,0,0,0,0,0,0,0],[A2,pi]) # in mm
	r2FCU = 20 # in mm
	r2FDS = 20 # in mm
	r2PL = 20 # in mm
	r2ECU = -20 # in mm
	r2EDM = -20 # in mm
	r2EDC = -20 # in mm

	MomentArmMatrix = np.matrix([[r1DELTa,r1CB,r1DELTp,r1BIC,r1TRI,r1BRA,r1BRD,r1PRO,r1FCR,r1ECRB,r1ECRL,r1FCU,r1FDS,r1PL,r1ECU,r1EDM,r1EDC],\
									[r2DELTa,r2CB,r2DELTp,r2BIC,r2TRI,r2BRA,r2BRD,r2PRO,r2FCR,r2ECRB,r2ECRL,r2FCU,r2FDS,r2PL,r2ECU,r2EDM,r2EDC]])
=======
	import sympy as sp
	global q1,q2,q_alt,RMatrix_Transpose
	q1,q2,q_alt = sp.symbols('q1'),sp.symbols('q2'),sp.symbols('q_alt')

	# #R_tranpose Column 1
	# r1DELTa = 19 # in mm
	# r1CB = 20 # in mm
	# r1DELTp = -8 # in mm
	# r1BIC = 15 # in mm
	# r1TRI = -15 # in mm
	# r1BRA = 0 # in mm
	# r1BRD = 0 # in mm
	# r1PRO = 0 # in mm
	# r1FCR = 0 # in mm
	# r1ECRB = 0 # in mm
	# r1ECRL = 0 # in mm
	# r1FCU = 0 # in mm
	# r1FDS = 0 # in mm
	# r1PL = 0 # in mm
	# r1ECU = 0 # in mm
	# r1EDM = 0 # in mm
	# r1EDC = 0 # in mm
	#
	# #R_tranpose Column 2
	# r2DELTa = 0 # in mm
	# r2CB = 0 # in mm
	# r2DELTp = 0 # in mm
	# r2BIC = MA_function([8.4533,36.6147,2.4777,-19.432,2.0571,0,13.6502,0,0,-5.6172,0,-2.0854,0,0,0,0], Angles=[A1,pi]) # in mm
	# r2TRI = MA_function([-24.5454,-8.8691,9.3509,-1.7518,0],Angles=A1) # in mm
	# r2BRA = MA_function([16.1991,-16.1463,24.5512,-6.3335,0],Angles=A1) # in mm
	# r2BRD = MA_function([15.2564,-11.8355,2.8129,-5.7781,44.8143,0,2.9032,0,0,-13.4956,0,-0.3940,0,0,0,0], Angles=[A1,pi]) # in mm
	# r2PRO = MA_function([11.0405,-1.0079,0.3933,-10.4824,-12.1639,-0.4369,36.9174,3.5232,-10.4223,21.2604,-37.2444,10.2666,-11.0060,14.5974,-3.9919,1.7526,-2.0089,0.5460], Angles=[A1,pi]) # in mm
	# r2FCR = 14 # in mm
	# r2ECRB = MA_function([-11.256,17.8548,1.6398,-0.5073,-2.8827,0,-0.0942,0,0,0,0,0,0,0,0,0], Angles=[A1,pi]) # in mm
	# r2ECRL = MA_function([-7.7034,16.3913,7.4361,-1.7566,0,-1.3336,0,0.0742,0,0,0,0,0,0,0,0], Angles=[A1,pi]) # in mm
	# r2FCU = 19 # in mm
	# r2FDS = 20 # in mm
	# r2PL = 25 # in mm
	# r2ECU = -23 # in mm
	# r2EDM = -10 # in mm
	# r2EDC = -20 # in mm
	#
	# MomentArmMatrix = np.matrix([[r1DELTa,r1CB,r1DELTp,r1BIC,r1TRI,r1BRA,r1BRD,r1PRO,r1FCR,r1ECRB,r1ECRL,r1FCU,r1FDS,r1PL,r1ECU,r1EDM,r1EDC],\
	# 								[r2DELTa,r2CB,r2DELTp,r2BIC,r2TRI,r2BRA,r2BRD,r2PRO,r2FCR,r2ECRB,r2ECRL,r2FCU,r2FDS,r2PL,r2ECU,r2EDM,r2EDC]])
	MomentArmMatrix = RMatrix_Transpose.subs([(q2,A2),(q_alt,pi)]).T
>>>>>>> 7d8b48796ec1a090f503dce659a61ec7087f6ba7
	return(MomentArmMatrix)
def calculate_muscle_velocities(A1,A2,Ȧ1,Ȧ2):
	import numpy as np
	import sympy as sp

	q2,q_alt = sp.symbols('q2'),sp.symbols('q_alt')

	R = return_MA_matrix(A1,A2)
	MuscleVelocity = -R.T*(np.matrix([Ȧ1,Ȧ2]).T)
	OptimalMuscleLength = np.array([(98),     (93),      (137),      (116),      (134),    (86),   \
	                  				(173),    (49),       (63),       (59),       (81),    (51),   \
	                   				(84),     (64),       (62),       (68),        (70)] ) # in mm
	NormalizedMuscleVelocity = np.array([[MuscleVelocity[i][0,0]/OptimalMuscleLength[i] for i in range(len(MuscleVelocity))]])
	return(NormalizedMuscleVelocity)
def eccentric_velocities(NormalizedMuscleVelocity):
	import numpy as np
	PositiveMuscleVelocities = [NormalizedMuscleVelocity.T[i,:][NormalizedMuscleVelocity.T[i,:]>0] for i in range(np.shape(NormalizedMuscleVelocity)[1])]
	return(PositiveMuscleVelocities)
def concentric_velocities(NormalizedMuscleVelocity):
	import numpy as np
	NegativeMuscleVelocities = [NormalizedMuscleVelocity.T[i,:][NormalizedMuscleVelocity.T[i,:]<0] for i in range(np.shape(NormalizedMuscleVelocity)[1])]
	return(NegativeMuscleVelocities)
def cost_function(X,type="avg"):
	assert type in ['avg','sos','l2norm','l1norm'], "type must be either 'avg','sos','l2norm', or 'l1norm'"
	if type == 'avg':
		cost = sum(X)/len(X)
	elif type == 'sos':
		cost = sum([el**2 for el in X])
	elif type == 'l1norm':
		cost = sum(X)
	elif type == 'l2norm':
		cost = sum([el**2 for el in X])**0.5
	return(cost)
def eccentric_cost(NormalizedMuscleVelocity,t_end = 1, dt = 0.001,type ='l2norm'):
	import numpy as np
	PositiveMuscleVelocities = eccentric_velocities(NormalizedMuscleVelocity)
	TotalPositiveExcursion = [np.trapz(el,dx=t_end*dt) for el in PositiveMuscleVelocities]
	EccentricCost = cost_function(TotalPositiveExcursion,type=type)
	return(EccentricCost)
def concentric_cost(NormalizedMuscleVelocity,t_end = 1, dt = 0.001,type = 'l2norm'):
	import numpy as np
	NegativeMuscleVelocities = concentric_velocities(NormalizedMuscleVelocity)
	TotalNegativeExcursion = [np.trapz(el,dx=t_end*dt) for el in NegativeMuscleVelocities]
	ConcentricCost = cost_function(TotalNegativeExcursion,type=type)
	return(ConcentricCost)
def return_initial_muscle_length(Coefficients,OptimalMuscleLength):
	"""
	Coefficients should be a list with the coefficients for each joint angle.
	"""
	import numpy as np
	from math import pi
	global A1,A2
	# Shoulder functions are constant. Therefore the input will only be A1[0] or 0.
	Angle1Initial = A1[0]
	Angle1Anatomical = 0
	if type(Coefficients[1]) in [float,int] or len(Coefficients[1]) == 5:
		Angle2Initial = A2[0]
		Angle2Anatomical = 0
	elif len(Coefficients[1]) == 16 or len(Coefficients[1]) == 18:
		# ISSUE! This does not take into account the length change that happens with the pronation!!
		Angle2Initial = [A2[0],pi]
		Angle2Anatomical = [0,pi]

	InitialLength = OptimalMuscleLength \
				- (MA_function_integral(Coefficients[0],Angles = Angle1Initial) - MA_function_integral(Coefficients[0],Angles = Angle1Anatomical)) \
				- (MA_function_integral(Coefficients[1],Angles = Angle2Initial) - MA_function_integral(Coefficients[1],Angles = Angle2Anatomical))
	return(InitialLength)
def calculate_muscle_lengths():
	import numpy as np
	from math import pi
	global MuscleLengths,A1,A2
	MuscleLengths = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

	DELTa_Coefficients = [19, 0]
	CB_Coefficients = [20, 0]
	DELTp_Coefficients = [-8, 0]
	BIC_Coefficients = [15, [8.4533,36.6147,2.4777,-19.432,2.0571,0,13.6502,0,0,-5.6172,0,-2.0854,0,0,0,0]]
	TRI_Coefficients = [-15, [-24.5454,-8.8691,9.3509,-1.7518,0]]
	BRA_Coefficients = [0, [16.1991,-16.1463,24.5512,-6.3335,0]]
	BRD_Coefficients = [0, [15.2564,-11.8355,2.8129,-5.7781,44.8143,0,2.9032,0,0,-13.4956,0,-0.3940,0,0,0,0]]
	PRO_Coefficients = [0, [11.0405,-1.0079,0.3933,-10.4824,-12.1639,-0.4369,36.9174,3.5232,-10.4223,21.2604,-37.2444,10.2666,-11.0060,14.5974,-3.9919,1.7526,-2.0089,0.5460]]
	FCR_Coefficients = [0, 14]
	ECRB_Coefficients = [0, [-11.256,17.8548,1.6398,-0.5073,-2.8827,0,-0.0942,0,0,0,0,0,0,0,0,0]]
	ECRL_Coefficients = [0, [-7.7034,16.3913,7.4361,-1.7566,0,-1.3336,0,0.0742,0,0,0,0,0,0,0,0]]
	FCU_Coefficients = [0, 19]
	FDS_Coefficients = [0, 20]
	PL_Coefficients = [0, 25]
	ECU_Coefficients = [0, -23]
	EDM_Coefficients = [0, -10]
	EDC_Coefficients = [0, -20]

	AllCoefficients = [DELTa_Coefficients, CB_Coefficients, DELTp_Coefficients, BIC_Coefficients, \
						TRI_Coefficients, BRA_Coefficients, BRD_Coefficients, PRO_Coefficients, \
						FCR_Coefficients, ECRB_Coefficients, ECRL_Coefficients, FCU_Coefficients, \
						FDS_Coefficients, PL_Coefficients, ECU_Coefficients, EDM_Coefficients, EDC_Coefficients]
	OptimalMuscleLength = np.array([(98),     (93),      (137),      (116),      (134),    (86),   \
                      				(173),    (49),       (63),       (59),       (81),    (51),   \
                       				(84),     (64),       (62),       (68),        (70)] )
	InitialMuscleLengths = [return_initial_muscle_length(AllCoefficients[i],OptimalMuscleLength[i]) for i in range(17)]

	for i in range(len(A1)):
		for j in range(17):
			# Shoulder functions are constant. Therefore the input will only be A1[0] or 0.
			Angle1Initial = A1[0]
			Angle1Current = A1[i]
			if type(AllCoefficients[j][1]) in [float,int] or len(AllCoefficients[j][1]) == 5:
				Angle2Initial = A2[0]
				Angle2Current = A2[i]
			elif len(AllCoefficients[j][1]) == 16 or len(AllCoefficients[j][1]) == 18:
				# ISSUE! This does not take into account the length change that happens with the pronation!!
				Angle2Initial = [A2[0],pi]
				Angle2Current = [A2[i],pi]
			TemporaryMuscleLength = InitialMuscleLengths[j]  \
									- (MA_function_integral(AllCoefficients[j][0],Angle1Current) - MA_function_integral(AllCoefficients[j][0],Angle1Initial)) \
									- (MA_function_integral(AllCoefficients[j][1],Angle2Current) - MA_function_integral(AllCoefficients[j][1],Angle2Initial))
			MuscleLengths[j].append(TemporaryMuscleLength)
import numpy as np
import matplotlib.pyplot as plt
import math
t_end = 0.5
dt = 0.0001
t = np.arange(0,1+dt,dt)*t_end

# 0.80*(L1+L2) = 0.4586
<<<<<<< HEAD

# Side to Side, Path Length = 0.8⋅(L1+L2) = 0.4586
# Xi = [-0.4586/2,0.1+0.4586/2]
# Xf = [0.4586/2,0.1+0.4586/2]

# Center Reach, Path Length = 0.8⋅(L1+L2) = 0.4586
# Xi = [0,0.1]
# Xf = [0,0.1 + 0.4586]

# Left Diagonal Reach, Path Length = 0.8⋅(L1+L2) = 0.4586
# Xi = [0,0.1]
# Xf = [-0.4586/(2**0.5),0.1 + 0.4586/(2**0.5)]

# Right Diagonal Reach, Path Length = 0.8⋅(L1+L2) = 0.4586
Xi = [0,0.1]
Xf = [0.4586/(2**0.5),0.1 + 0.4586/(2**0.5)]



=======
Xi = [0,0.1]#[-0.4586/2,0.4586/2]
Xf = [0,0.1+0.4586]#[0.4586/2,0.4586/2]# [0.4586/(2**0.5),0.1 + 0.4586/(2**0.5)]
>>>>>>> 7d8b48796ec1a090f503dce659a61ec7087f6ba7
# forward

reaching_task(Xi=Xi, Xf=Xf,dt=dt, t_end=t_end)
calculate_muscle_lengths()
MuscleLengths_Forward = MuscleLengths
NormalizedMuscleVelocity_Forward = calculate_muscle_velocities(A1[0],A2[0],Ȧ1[0],Ȧ2[0])
for i in range(1,len(A1)):
	NormalizedMuscleVelocity_Forward = np.concatenate((NormalizedMuscleVelocity_Forward,calculate_muscle_velocities(A1[i],A2[i],Ȧ1[i],Ȧ2[i])),axis=0)

plt.figure()
color = iter(plt.cm.rainbow(np.linspace(0,1,17)))
[plt.plot(t,NormalizedMuscleVelocity_Forward[:,i],c=next(color)) for i in range(17)]
ax1a = plt.gca()
ax1a.set_xlim(0,t_end*(1.3))
ax1a.set_title('Forward Direction')
if t_end!=1:
	ax1a.set_xlabel('Time (s)')
else:
	ax1a.set_xlabel('Normalized Time')
ax1a.set_ylabel('Normalized Muscle Velocity\nConcentric $\longleftrightarrow$ Eccentric')
ax1a.legend(["DELTa","CB","DELTp","BIC","TRI","BRA","BRD","PRO","FCR","ECRB","ECRL","FCU","FDS","PL","ECU","EDM","EDC"])

plt.figure()
plt.plot(t,A1,'g')
plt.plot(t,A2,'g:')
ax1b = plt.gca()
ax1b.set_xlim(0,t_end*(1.3))
ax1b.set_title('Forward Direction')
if t_end!=1:
	ax1b.set_xlabel('Time (s)')
else:
	ax1b.set_xlabel('Normalized Time')
ax1b.set_ylabel('Joint Angles (in radians)')
ax1b.legend(["Shoulder","Elbow"])

plt.figure()
plt.plot(A1,A2,'g')
ax1c = plt.gca()
ax1c.set_title("Forward Direction\nConfiguration Space Trajectory")
ax1c.set_ylabel("Shouler Angle (in radians)")
ax1c.set_xlabel("Elbow Angle (in radians)")

EccentricCost_Forward = eccentric_cost(NormalizedMuscleVelocity_Forward,t_end=t_end,dt = dt)
ConcentricCost_Forward = concentric_cost(NormalizedMuscleVelocity_Forward,t_end=t_end,dt = dt)

plt.figure()
color = iter(plt.cm.rainbow(np.linspace(0,1,17)))
[plt.plot(t,MuscleLengths_Forward[i],c=next(color)) for i in range(17)]
<<<<<<< HEAD
=======
# [plt.plot(t,MuscleLengths_Forward[i]) for i in range(17)]
>>>>>>> 7d8b48796ec1a090f503dce659a61ec7087f6ba7
ax1d = plt.gca()
ax1d.set_xlim(0,t_end*(1.3))
ax1d.set_title('Forward Direction')
if t_end!=1:
	ax1d.set_xlabel('Time (s)')
else:
	ax1d.set_xlabel('Normalized Time')
ax1d.set_ylabel('Muscle Lengths (in mm)')
ax1d.legend(["DELTa","CB","DELTp","BIC","TRI","BRA","BRD","PRO","FCR","ECRB","ECRL","FCU","FDS","PL","ECU","EDM","EDC"])


# reverse
reaching_task(Xi=Xf, Xf=Xi,dt=dt, t_end=t_end)
calculate_muscle_lengths()
MuscleLengths_Reverse = MuscleLengths
NormalizedMuscleVelocity_Reverse = calculate_muscle_velocities(A1[0],A2[0],Ȧ1[0],Ȧ2[0])
for i in range(1,len(A1)):
	NormalizedMuscleVelocity_Reverse = np.concatenate((NormalizedMuscleVelocity_Reverse,calculate_muscle_velocities(A1[i],A2[i],Ȧ1[i],Ȧ2[i])),axis=0)

plt.figure()
<<<<<<< HEAD
=======
# plt.plot(t,NormalizedMuscleVelocity_Reverse)
>>>>>>> 7d8b48796ec1a090f503dce659a61ec7087f6ba7
color = iter(plt.cm.rainbow(np.linspace(0,1,17)))
[plt.plot(t,NormalizedMuscleVelocity_Reverse[:,i],c=next(color)) for i in range(17)]
ax2a = plt.gca()
ax2a.set_xlim(0,t_end*(1.3))
ax2a.set_title('Reversed Direction')
if t_end!=1:
	ax2a.set_xlabel('Time (s)')
else:
	ax2a.set_xlabel('Normalized Time')
ax2a.set_ylabel('Normalized Muscle Velocity\nConcentric $\longleftrightarrow$ Eccentric')
ax2a.legend(["DELTa","CB","DELTp","BIC","TRI","BRA","BRD","PRO","FCR","ECRB","ECRL","FCU","FDS","PL","ECU","EDM","EDC"])

plt.figure()
plt.plot(t,A1,'b')
<<<<<<< HEAD
plt.plot(t,A2,'b--')
=======
plt.plot(t,A2,'b:')
>>>>>>> 7d8b48796ec1a090f503dce659a61ec7087f6ba7
ax2b = plt.gca()
ax2b.set_xlim(0,t_end*(1.3))
ax2b.set_title('Reversed Direction')
if t_end!=1:
	ax2b.set_xlabel('Time (s)')
else:
	ax2b.set_xlabel('Normalized Time')
ax2b.set_ylabel('Joint Angles (in radians)')
ax2b.legend(["Shoulder","Elbow"])

plt.figure()
plt.plot(A1,A2,'b')
ax2c = plt.gca()
ax2c.set_title("Reversed Direction\nConfiguration Space Trajectory")
ax2c.set_ylabel("Shouler Angle (in radians)")
ax2c.set_xlabel("Elbow Angle (in radians)")

EccentricCost_Reverse = eccentric_cost(NormalizedMuscleVelocity_Reverse,t_end=t_end,dt = dt)
ConcentricCost_Reverse = concentric_cost(NormalizedMuscleVelocity_Reverse,t_end=t_end,dt = dt)

plt.figure()
<<<<<<< HEAD
=======
# [plt.plot(t,MuscleLengths_Reverse[i]) for i in range(17)]
>>>>>>> 7d8b48796ec1a090f503dce659a61ec7087f6ba7
color = iter(plt.cm.rainbow(np.linspace(0,1,17)))
[plt.plot(t,MuscleLengths_Reverse[i],c=next(color)) for i in range(17)]
ax2d = plt.gca()
ax2d.set_xlim(0,t_end*(1.3))
ax2d.set_title('Reverse Direction')
if t_end!=1:
	ax2d.set_xlabel('Time (s)')
else:
	ax2d.set_xlabel('Normalized Time')
ax2d.set_ylabel('Muscle Lengths (in mm)')
ax2d.legend(["DELTa","CB","DELTp","BIC","TRI","BRA","BRD","PRO","FCR","ECRB","ECRL","FCU","FDS","PL","ECU","EDM","EDC"])


plt.figure()
plt.plot(np.arange(0,2*math.ceil(max([EccentricCost_Forward,ConcentricCost_Forward]))+1,1), np.arange(0,2*math.ceil(max([EccentricCost_Forward,ConcentricCost_Forward]))+1,1),'0.75',linestyle='--')
plt.scatter([EccentricCost_Forward],[ConcentricCost_Forward],facecolor ='g', edgecolor = 'k', marker='o',s=100)
plt.scatter([EccentricCost_Reverse],[ConcentricCost_Reverse],facecolor ='b', edgecolor = 'k', marker='o',s=100)
ax3 = plt.gca()
ax3.set_xlim(0,2*math.ceil(max([EccentricCost_Forward,EccentricCost_Reverse])))
ax3.set_ylim(0,2*math.ceil(max([EccentricCost_Forward,EccentricCost_Reverse])))
ax3.set_aspect('equal')
ax3.set_title('Eccentric vs. Concentric Cost\nFor Forward and Reverse Movments')
ax3.set_xlabel('Eccentric Cost (in $\hat{l}_{o}$)')
ax3.set_ylabel('Concentric Cost (in $\hat{l}_{o}$)')
ax3.legend(['Equal Cost Line','Forward Direction','Reversed Direction'])
<<<<<<< HEAD
#

# Using Default Reaching Directions Listed Above...
# CostValues ={'CenterForward': np.array([[ 1.3220062 ,  2.10011198]]),\
# 				'CenterReverse': np.array([[ 2.10011198,  1.3220062 ]]),\
# 					'LeftForward': np.array([[ 1.24643599,  1.88616011]]),\
# 						'LeftReverse': np.array([[ 1.88616011,  1.24643599]]),\
# 							'RightForward': np.array([[ 1.1047386 ,  1.92105444]]),\
# 								'RightReverse': np.array([[ 1.92105444,  1.1047386 ]]),\
# 									'SideForward': np.array([[ 0.37931692,  0.53156831]]),\
# 										'SideReverse': np.array([[ 0.53156831,  0.37931692]])}
# plt.figure()
# MaximumCost = max([max(CostValues[key][0]) for key in CostValues.keys()])
# plt.plot(np.arange(0,2*math.ceil(MaximumCost)+1,1), np.arange(0,2*math.ceil(MaximumCost)+1,1),'0.75',linestyle='--')
# facecolor = iter(['k','#9f9f9f','b','#9f9f9f','r','#9f9f9f','g','#9f9f9f'])
# edgecolor = iter(['k','k','k','b','k','r','k','g'])
# [plt.scatter(CostValues[key][0,0],CostValues[key][0,1],c=next(facecolor),edgecolor=next(edgecolor),marker='o',s=75) for key in CostValues.keys()]
# ax4 = plt.gca()
# ax4.set_xlim(0,math.ceil(MaximumCost))
# ax4.set_ylim(0,math.ceil(MaximumCost))
# ax4.set_aspect('equal')
# ax4.set_title('Eccentric vs. Concentric Cost\nFor Forward and Reverse Movments')
# ax4.set_xlabel('Eccentric Cost (in $\hat{l}_{o}$)')
# ax4.set_ylabel('Concentric Cost (in $\hat{l}_{o}$)')
# # ax4.legend(['Equal Cost Line','Center (Forward)','Center (Reverse)','Left (Forward)','Left (Reverse)','Right (Forward)','Right (Reverse)','Side-to-Side (Forward)','Side-to-Side (Reverse)'],loc='upper right')
#
=======

>>>>>>> 7d8b48796ec1a090f503dce659a61ec7087f6ba7
plt.show()
