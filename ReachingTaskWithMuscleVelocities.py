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
	if Title != '' : Title = ' '*(22-len(Title)) + Title + ' : '
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
def set_link_lengths(New_L1=0.256,New_L2=0.405,EOM = "Uno"):
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
def set_initial_and_final_positions(ReachType):
	"""
	This takes in a string -- either 'Center','Right','Left', or 'Sideways' -- and returns the necessary initial and final positions for the movement, based on Flash/Hogan (1985).

	Parameters:
	# 0.80*(L1+L2) = 0.4586
	# 0.35/(L1+L2) = 0.5295
	# Sternum at -0.177 = -L1*(0.129/0.186) <-- Anthropomorphic Ratio
	"""

	assert ReachType.capitalize() in ['Center','Right','Left','Sideways'], \
		"ReachType must be either 'Center','Right','Left', or 'Sideways'."


	# Side to Side, Path Length = 0.5295⋅(L1+L2) = 0.35
	if ReachType.capitalize() == 'Sideways':
		Xi = [-0.177-0.35/2,0.1+0.35/2]
		Xf = [-0.177+0.35/2,0.1+0.35/2]

	# Center Reach, Path Length = 0.5295⋅(L1+L2) = 0.35
	elif ReachType.capitalize() == 'Center':
		Xi = [-0.177,0.1]
		Xf = [-0.177,0.1 + 0.35]

	# # Left Diagonal Reach, Path Length = 0.5295⋅(L1+L2) = 0.35
	elif ReachType.capitalize() == 'Left':
		Xi = [-0.177,0.1]
		Xf = [-0.177-0.35/(2**0.5),0.1 + 0.35/(2**0.5)]

	# # Right Diagonal Reach, Path Length = 0.5295⋅(L1+L2) = 0.35
	elif ReachType.capitalize() == 'Right':
		Xi = [-0.177,0.1]
		Xf = [-0.177+0.35/(2**0.5),0.1 + 0.35/(2**0.5)]

	return(Xi,Xf)
def find_X_values(t,Xi,Xf):
	"""
	This takes t (a numpy.ndarray of normalized time ∈ [0,1]) and either a (2,) or (2,1) list/array with values for both initial and final x and y (endpoint). To avoid singularities, ||X[i]|| cannot be greater than L1 + L2.
	"""
	import numpy as np
	import ipdb
	assert len(Xf)==2, "Xf must be either a (2,) or (2,1) list/array"
	assert len(Xi)==2, "Xi must be either a (2,) or (2,1) list/array"
	xi,yi = Xi
	xf,yf = Xf
	if np.shape(Xi) == (2,1): xi,yi = xi[0],yi[0]
	if np.shape(Xf) == (2,1): xf,yf = xf[0],yf[0]
	x = xi + (xf-xi)*(10*t**3 - 15*t**4 + 6*t**5)
	y = yi + (yf-yi)*(10*t**3 - 15*t**4 + 6*t**5)
	X = np.array(np.concatenate(([x],[y]),axis=0))
	# ipdb.set_trace()
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
def calculate_torques(EOM="Uno"):
	from math import cos,sin
	import numpy as np
	assert EOM in ["Uno","Zadravec"], "EOM can be either 'Uno' or 'Zadravec'"
	if EOM == "Uno": # Uno et al. Biological Cybernetics (1989)
	    m1,m2 = 1.02,1.16 # kg
	    c1,c2 = 0.104,0.165 # m
	    I1,I2 = 0.0167,0.0474 # kg⋅m²
	    b11,b12,b21,b22 = 0.8,0,0,0.8
	    α = I1 + I2 + m2*(L1**2)
	    β = m2*L1*c2
	    δ = I2
	else: # Zadravec, Biocybernetics and Biomedical Engineering (2013)
	    m1,m2 = 2.089,1.912 # kg
	    c1,c2 = 0.152,0.181 # m
	    I1,I2 = 0.0159,0.0257 # kg⋅m²
	    b11,b12,b21,b22 = 0.74,0.10,0.10,0.82
	    α = I1 + I2 + m1*(c1**2) + m2*(L1**2 + c2**2)
	    β = m2*L1*c2
	    δ = I2 + m2*(c2**2)
	C = lambda a1,a2,ȧ1,ȧ2: \
	        np.matrix([ [-β*ȧ2*sin(a2),     -β*(ȧ1 + ȧ2)*sin(a2)],
	                    [β*ȧ1*sin(a2),      0]]) # kg⋅m² (N⋅m⋅s²)
	M = lambda a1,a2: \
	        np.matrix([ [α + 2*β*cos(a2),   δ + β*cos(a2)],
	                    [δ + β*cos(a2),     δ]]) # kg⋅m² (N⋅m⋅s²)
	B = np.matrix([ [b11, b12],\
	                [b21, b22]]) # kg⋅m²/s (N⋅m⋅s)
	global T1,T2
	T1,T2 = [],[]
	for i in range(len(A1)):
	    Ȧ = np.matrix([[Ȧ1[i]],[Ȧ2[i]]])
	    Ä = np.matrix([[Ä1[i]],[Ä2[i]]])
	    T = M(A1[i],A2[i])*Ä + C(A1[i],A2[i],Ȧ1[i],A2[i])*Ȧ + B*Ȧ
	    T1.append(T[0,0])
	    T2.append(T[1,0])
def calculate_potential_variability(i,δT1,δT2,dt=0.001,EOM = 'Uno',Group = 'flexor'):
	from math import cos,sin
	import numpy as np
	assert EOM in ["Uno","Zadravec"], "EOM can be either 'Uno' or 'Zadravec'"

	global A1,A2,Ȧ1,Ȧ2,T1,T2

	if EOM == "Uno": # Uno et al. Biological Cybernetics (1989)
	    m1,m2 = 1.02,1.16 # kg
	    c1,c2 = 0.104,0.165 # m
	    I1,I2 = 0.0167,0.0474 # kg⋅m²
	    b11,b12,b21,b22 = 0.8,0,0,0.8
	    α = I1 + I2 + m2*(L1**2)
	    β = m2*L1*c2
	    δ = I2
	else: # Zadravec, Biocybernetics and Biomedical Engineering (2013)
	    m1,m2 = 2.089,1.912 # kg
	    c1,c2 = 0.152,0.181 # m
	    I1,I2 = 0.0159,0.0257 # kg⋅m²
	    b11,b12,b21,b22 = 0.74,0.10,0.10,0.82
	    α = I1 + I2 + m1*(c1**2) + m2*(L1**2 + c2**2)
	    β = m2*L1*c2
	    δ = I2 + m2*(c2**2)

	if Group == 'flexor':
		δT1== δT1
		δT2== δT2
	elif Group == 'extensor':
		δT1== -δT1
		δT2== -δT2

	# if Group[0] == 'flexor':
	# 	δτ1== δτ1
	# elif Group[0] == 'extensor':
	# 	δτ1== -δτ1
	#
	# if Group[1] == 'flexor':
	# 	δτ2== δτ2
	# elif Group[1] == 'extensor':
	# 	δτ2== -δτ2

	ä1 = lambda a1,a2,ȧ1,ȧ2,τ1,τ2,δτ1,δτ2: \
				(-(β**2)*(ȧ1**2)*sin(a2)*cos(a2)	- ȧ1*(b11*δ - b21*(δ + β*cos(a2))) \
				+ δ*β*(ȧ1 + ȧ2)**2*sin(a2) - ȧ2*(b12*δ - b22*(δ + β*cos(a2))) \
				+ (τ1 + δτ1)*δ - (τ2 + δτ2)*(δ + β*cos(a2))) / \
				(α*δ - δ**2 -β**2*cos(a2)**2)
	ä2 = lambda a1,a2,ȧ1,ȧ2,τ1,τ2,δτ1,δτ2: \
				-((β**2)*(2*ȧ1**2 + 2*ȧ1*ȧ2 + ȧ2**2)*sin(a2)*cos(a2) \
				- ȧ1*(b11*(δ + β*cos(a2))  - b21*(α + 2*β*cos(a2))) \
				+ δ*β*(ȧ1 + ȧ2)**2*sin(a2) + β*ȧ1**2*(α-δ)*sin(a2) \
				- ȧ2*(b12*(δ + β*cos(a2)) - b22*(α + 2*β*cos(a2))) \
				+ (τ1 + δτ1)*(δ + β*cos(a2)) - (τ2 + δτ2)*(α + 2*β*cos(a2))) / \
				(α*δ - δ**2 -β**2*cos(a2)**2)
	def next_angle_from_forward_euler(i,δT1,δT2,dt):
		import numpy as np
		global A1,A2,Ȧ1,Ȧ2,T1,T2
		next_ȧ1	= Ȧ1[i] + ä1(A1[i],A2[i],Ȧ1[i],Ȧ2[i],T1[i],T2[i],δT1[i],δT2[i])*dt
		next_a1 = A1[i] + next_ȧ1*dt
		next_ȧ2	= Ȧ2[i] + ä2(A1[i],A2[i],Ȧ1[i],Ȧ2[i],T1[i],T2[i],δT1[i],δT2[i])*dt
		next_a2 = A2[i] + next_ȧ2*dt
		return(next_a1,next_a2)

	# testA1,testA2 = [],[]
	# for i in range(len(A1)):
	# 	next_a1,next_a2 = next_angle_from_forward_euler(A1[i],A2[i],Ȧ1[i],Ȧ2[i],T1[i],T2[i],0,0,0.001)
	# 	testA1.append(next_a1)
	# 	testA2.append(next_a2)
	next_a1,next_a2 = next_angle_from_forward_euler(i,δT1,δT2,dt)
	X = np.array([	[L1*cos(A1[i]) + L2*cos(A1[i] + A2[i])],\
					[L1*sin(A1[i]) + L2*sin(A1[i] + A2[i])] ])
	new_X = np.array([	[L1*cos(next_a1) + L2*cos(next_a1 + next_a2)],\
						[L1*sin(next_a1) + L2*sin(next_a1 + next_a2)] ])
	# global A1_dev,A2_dev
	# A1_dev.append(A1[i]-next_a1)
	# A2_dev.append(A2[i]-next_a2)
	variation = float((sum((X-new_X)**2)**0.5))
	return(variation)
def reaching_task(dt = 0.001, Xi = [0,0.25], Xf = [0.25,0.50],t_end=1):
	import numpy as np
	set_link_lengths()
	create_angle_lists()
	t = np.arange(0,1+dt,dt)
	X,Ẋ,Ẍ = return_X_values(t,Xi,Xf,t_end)
	update_angle_lists(X,Ẋ,Ẍ)
	calculate_torques()
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

	src = Parameters['src']
	Coefficients = Parameters['MA']
	eq = Parameters['eq']
	dof = Parameters['dof']
	threshold = Parameters['threshold']

	global q1,q2,q_PS
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
			q = q2
		else:
			q = q1
		MomentArm = (np.matrix(Coefficients)*np.matrix([1,q,q**2,q**3,q**4,q**5]).T)[0,0]
	elif src.capitalize() == 'Est' :
		MomentArm = Coefficients
	else: #src.capitalize() == 'Ramsay'
		q = q2
		assert type(Coefficients) == list, "Coefficients must be a list."
		assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
		assert eq in [1,2,3], "eq must be either 1, 2, or 3 when using Ramsay (2009)."
		if eq == 1:
			assert len(Coefficients) == 5, "For Eq. 1, Coefficients must be 5 elements long."
			MomentArm = (sp.Matrix(Coefficients).T*sp.Matrix([1,q,q**2,q**3,q**4]))[0,0]
		elif eq == 2:
			assert len(Coefficients)==16, "For Eq. 2, Coefficients must be 16 elements long."
			MomentArm = (sp.Matrix(Coefficients).T*\
							sp.Matrix([1, q, q_PS, q*q_PS, q**2, \
										q_PS**2, (q**2)*q_PS, q*(q_PS**2), \
										(q**2)*(q_PS**2), q**3, q_PS**3, \
										(q**3)*q_PS, q*(q_PS**3), \
										(q**3)*(q_PS**2), (q**2)*(q_PS**3), \
										(q**3)*(q_PS**3)]))[0, 0]
		else: # eq == 3
			assert len(Coefficients)==18, "For Eq. 3, Coefficients must be 18 elements long."
			MomentArm = (sp.Matrix(Coefficients).T*\
							sp.Matrix([1, q, q_PS, q*q_PS, q**2, \
								q_PS**2, (q**2)*q_PS, q*(q_PS**2), (q**2)*(q_PS**2), \
								q**3, (q**3)*q_PS, (q**3)*(q_PS**2), \
								q**4, (q**4)*q_PS, (q**4)*(q_PS**2),  \
								q**5, (q**5)*q_PS, (q**5)*(q_PS**2)]))[0, 0]
	if threshold == None:
		return(MomentArm)
	else:
		assert type(threshold) in [int,float], "threshold must be a number."
		MomentArm = sp.Piecewise((MomentArm,q<threshold),(MomentArm.subs(q,threshold),q>=threshold))
		return(MomentArm)
def MA_function_integral(Coefficients,Angles = None):
	"""
	Note: Angles should be a number if Coefficients has a length of 1 or 5, or a list of length 2 when the Coefficients have lengths 16 or 18. a1 will be the PRIMARY ANGLE for the DOF being considered while Angles[1] will be the secondary angle.
	"""
	import numpy as np
	import sympy as sp
	import ipdb
	assert type(Coefficients) == list or type(Coefficients) == int or type(Coefficients) == float, "Coefficients must be a list, int, or float."
	if type(Coefficients) not in [int,float]:
		assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
	if type(Coefficients) in [int,float]:
		q2 = sp.symbols('q2')
		MomentArmIntegral = Coefficients*q2
	elif len(Coefficients)==5:
		q2 = sp.symbols('q2')
		integral = np.multiply(sp.Matrix([1,q2,q2**2,q2**3,q2**4]).T,(1,1/2,1/3,1/4,1/5))*q2
		MomentArmIntegral = (np.matrix(Coefficients)*integral.T)[0,0]
	elif len(Coefficients)==16:
		q2,q_PS = sp.symbols('q2'),sp.symbols('q_PS')
		integral = np.multiply(sp.Matrix([1, q2, q_PS, q2*q_PS, q2**2, \
						q_PS**2, (q2**2)*q_PS, q2*(q_PS**2), \
						(q2**2)*(q_PS**2), q2**3, q_PS**3, \
						(q2**3)*q_PS, q2*(q_PS**3), \
						(q2**3)*(q_PS**2), (q2**2)*(q_PS**3), \
						(q2**3)*(q_PS**3)]).T,\
						[1,1/2,1,1/2,1/3,1,1/3,1/2,1/3,1/4,1,1/4,1/2,1/4,1/3,1/4])*q2
		MomentArmIntegral = (np.matrix(Coefficients)*integral.T)[0, 0]
	else: # len(Coefficients)==18
		q2,q_PS = sp.symbols('q2'),sp.symbols('q_PS')
		integral = np.multiply(sp.Matrix([1, q2, q_PS, q2*q_PS, q2**2, \
							q_PS**2, (q2**2)*q_PS, q2*(q_PS**2), (q2**2)*(q_PS**2), \
							q2**3, (q2**3)*q_PS, (q2**3)*(q_PS**2), \
							q2**4, (q2**4)*q_PS, (q2**4)*(q_PS**2),  \
							q2**5, (q2**5)*q_PS, (q2**5)*(q_PS**2)]).T,\
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
			return(MomentArmIntegral.subs([(q2,Angles[0]),(q_PS,Angles[1])]))
		else: # len(Coefficients)==18
			return(MomentArmIntegral.subs([(q2,Angles[0]),(q_PS,Angles[1])]))
def Pigeon_coeff_conversion(Coefficients):
	"""
	Takes in Coefficient values from Pigeon (1996) -- which take in angles in degrees -- and coverts them into the properly scaled coefficients for radians, additionally scaled by the magnitude listed in the paper.

	Note that the coefficients listed in Pigeon (1996) are given in decending order (i.e., c₅,c₄,c₃,c₂,c₁,c₀). However to maintain continuity with the equations given in Ramsay (2009), we list coefficients in order of increasing power (i.e., c₀,c₁,c₂,c₃,c₄,c₅).
	"""
	import numpy as np
	assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
	assert type(Coefficients)==list, 'Coefficients must be a 6 element list.'
	rad_conversion = np.multiply(Coefficients,\
			np.array([1,(180/np.pi),(180/np.pi)**2,(180/np.pi)**3,(180/np.pi)**4,(180/np.pi)**5]))
	new_Coefficients = np.multiply(rad_conversion,[1,1e-1,1e-3,1e-5,1e-7,1e-9])
	return(new_Coefficients)
def global_R_matrix():
	"""
	Notes:
	Coefficients from observation, Ramsay, FVC, Holtzbaur, Pigeon, Kuechle, or Banks.

	BRA EFE MA for Ramsay has R² = 0.990 whereas Pigeon has R² = 0.9988. Curve appears to be a better fit, as it experiences its smallest MA when Elbow angle = 0. Coefficients and equation number/type are listed below to test either implementation.

	src = 'Ramsay', eq = 1, Coefficients = [16.1991,-16.1463,24.5512,-6.3335,0], threshold = None
	src = 'Pigeon', dof = 'elbow', eq = None, Coefficients = Pigeon_coeff_conversion([5.5492,2.3080,2.3425,-2.0530,0,0]), threshold = None

	BRD for Ramsay has R² = 0.988 whereas Pigeon has R² = 0.9989. Pigeon, however, only takes elbow angle into account, whereas Ramsay takes in variable PS angles. Coefficients and equation number/type are listed below to test either implementation.

	src = 'Ramsay', eq = 2, Coefficients = [15.2564,-11.8355,2.8129,-5.7781,44.8143,0,2.9032,0,0,-13.4956,0,-0.3940,0,0,0,0]
	src = 'Pigeon', dof = 'elbow', eq = None, Coefficients = Pigeon_coeff_conversion([19.490,1.6681,10.084,-6.5171,0,0])

	FCR EFE MA is not listed in Ramsay but Pigeon has a quadratic function with R² = 0.9975. Pigeon only takes elbow angle into account. Coefficients and equation number/type are listed below to test either implementation. EFE MA was estimated to be constant and 10 mm for this muscle. If you use Pigeon, make sure to only accept positive moment arm values, as this model fails outside the ROM. One option is to set the MA to the constant value (i.e., MA[theta>140°] = MA[140°])

	src = 'est', eq = 'constant', Coefficients = [10], clip = None
	src = 'Pigeon', dof = 'elbow', eq = None, Coefficients = Pigeon_coeff_conversion([0.9351,0.5375,-0.3627,0,0,0]),threshold = 2.86

	ECRL EFE MA for Ramsay has R² = 0.978 whereas Pigeon has R² = 0.9986. Pigeon, however, only takes elbow angle into account, whereas Ramsay takes in variable PS angles. Additionally, Pigeon only considers elbow angles between 0 and 140 degrees and exhibits a decrease in MA as elbow angle approaches the upper bound of the ROM. This should (intiutively speaking) make the extensors MA largest, but Pigeon exhibits a drop off that may make it less favorable for movements at the boundary of the ROM. Coefficients and equation number/type are listed below to test either implementation.

	src = 'Ramsay', eq = 2, Coefficients = [-7.7034,16.3913,7.4361,-1.7566,0,-1.3336,0,0.0742,0,0,0,0,0,0,0,0]
	src = 'Pigeon', dof = 'elbow', eq = None, Coefficients = Pigeon_coeff_conversion([4.7304,1.2590,4.4347,-3.0229,0,0])

	ECU EFE MA is not listed in Ramsay but Pigeon has a quadratic function with R² = 0.9966. Pigeon only takes elbow angle into account. Coefficients and equation number/type are listed below to test either implementation. EFE MA was estimated to be constant and -10 mm for this muscle. If you use Pigeon, make sure to only accept negative moment arm values, as this model fails outside the ROM. One option is to set the MA to the constant value (i.e., MA[theta>140°] = MA[140°])

	src = 'est', eq = 'constant', Coefficients = [-10]
	src = 'Pigeon', dof = 'elbow', eq = None, Coefficients = Pigeon_coeff_conversion([-2.1826,-1.7386,1.1491,0,0,0])

	BIC EFE MA for Ramsay has R² = 0.985 whereas Pigeon has R² = 0.9918. Pigeon, however, only takes elbow angle into account, whereas Ramsay takes in variable PS angles. It appears that because Pigeon uses an average of fully pronated and fully supinated MAs, the BIC moment arm is similar but subject to variation as the level of PS is changed. Coefficients and equation number/type are listed below to test either implementation. (NOTE: BIC becomes slightly negative when q2 > 3.021. If trajectory has elbow angles exceding this value, enter a threshold of 3.021 into the model.)

	Additionally, the SFE MA for the BIC is held constant in Pigeon at 29.21 mm while it was estimated as 15 mm.

	src = 'Ramsay', eq = 2, Coefficients = [8.4533,36.6147,2.4777,-19.432,2.0571,0,13.6502,0,0,-5.6172,0,-2.0854,0,0,0,0], threshold = 3.021
	src = 'Pigeon', dof = 'elbow', eq = None, Coefficients = Pigeon_coeff_conversion([14.660,4.5322,1.8047,-2.9883,0,0]), threshold = 2.9326

	TRI EFE MA for Ramsay has R² = 0.997 whereas Pigeon has R² = 0.9904. Pigeon appears to really fail when the elbow angle is greater than 140°. For this reason, Ramsay should be used. However the approach of fixing the MA for values greater than 140° can be adopted for completeness. Coefficients and equation number/type are listed below to test either implementation.

	Additionally, the SFE MA for the TRI is held constant in Pigeon at -25.40 mm while it was estimated as -15 mm.

	src = 'Ramsay', eq = 1, Coefficients = [-24.5454,-8.8691,9.3509,-1.7518,0]
	src = 'Pigeon', dof = 'elbow', eq = None, Coefficients = Pigeon_coeff_conversion([-23.287,-3.0284,12.886,-19.092,13.277,-3.5171])

	DELTa SFE MA is listed as 33.02 mm in Pigeon and estimated as 19 mm. Using Pigeon Coefficients convention, Kuechle (1997) has the DELTp MA for [-140,90] as Pigeon_coeff_conversion([ 1.34293189,  0.20316226, -0.02339031,  0.27807828,  0.,  0.]). This will yield a piecewise function that creates jumps with the new velocity formulation. Instead, we are going to try Pigeon_coeff_conversion([ 1.27928795,  0.20480346,  0.08917734,  0.32207214, -0.23928223,  0.]) so that the function is within range and continuous during the ROM. Threshold (pi/2) has been removed for this new MA function.

	DELTp SFE MA is listed as -78.74 mm in Pigeon and estimated as -8 mm. Using Pigeon Coefficients convention, Kuechle (1997) has the DELTp MA for [-140,90] as Pigeon_coeff_conversion([ 2.28547177,  0.39721238, -0.33900829, -0.36146546,  0.,  0.]). This will yield a piecewise function that creates jumps with the new velocity formulation. Instead, we are going to try Pigeon_coeff_conversion([-2.38165173, -0.4486164 ,  0.58655808,  0.65003255, -0.82736695,0.20812998]) so that the function is within range and continuous during the ROM. Threshold (pi/2) has been removed for this new MA function.

	PC (Clavicle attachment of Pectoralis) SFE MA is listed as 50.80 mm in Pigeon.

	CB SFE MA was estimated in Holzbaur (2005) as 20 mm while Bassett (1990) estimates from 7 cadavers the MA to be 36 mm.
	"""
	import sympy as sp
	import numpy as np
	from numpy import pi

	global q1,q2,q_PS
	q1,q2,q_PS = sp.symbols('q1'),sp.symbols('q2'),sp.symbols('q_PS')

	# Coefficients from observation, Ramsay, Pigeon, FVC, Holtzbaur, or Banks.
	# Moment arms are in mm. Mass is in grams. threshold is in radians.
	PC_Coefficients = {\
		'Shoulder' : {\
			'MA' : [50.80,0,0,0,0,0],\
			'src' : 'Pigeon', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : None, \
		'Actual No' : None, \
		'Corrected No' : None, \
		'Relative Abundance' : None,\
		'Optimal Muscle Length' : None,\
	    'Group' : 'flexor'}
	DELTa_Coefficients = {\
		'Shoulder' : {\
			'MA' : Pigeon_coeff_conversion([ 1.27928795,  0.20480346,  0.08917734,  0.32207214, -0.23928223,  0.        ]),\
			'src' : 'Pigeon', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 355.7/3, \
		'Actual No' : 182/3, \
		'Corrected No' : 426.3/3, \
		'Relative Abundance' : 0.43,\
		'Optimal Muscle Length' : 98,\
	    'Group' : 'flexor'}
	CB_Coefficients = {\
		'Shoulder' : {\
			'MA' : 20, \
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 39.8, \
		'Actual No' : 123, \
		'Corrected No' : 147.3, \
		'Relative Abundance' : 0.83,\
		'Optimal Muscle Length' : 93,\
	    'Group' : 'flexor'}
	DELTp_Coefficients = {\
		'Shoulder' : {\
			'MA' : Pigeon_coeff_conversion([-2.38165173, -0.4486164 ,  0.58655808,  0.65003255, -0.82736695,0.20812998]), \
			'src' : 'Pigeon', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 355.7/3, \
		'Actual No' : 182/3, \
		'Corrected No' : 426.3/3, \
		'Relative Abundance' : 0.43,\
		'Optimal Muscle Length' : 137,\
	    'Group' : 'extensor'}
	BIC_Coefficients = {\
		'Shoulder' : {\
			'MA' : [29.21,0,0,0,0,0],\
			'src' : 'Pigeon', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : [8.4533,36.6147,2.4777,-19.432,2.0571,0,13.6502,0,0,-5.6172,0,-2.0854,0,0,0,0],\
			'src' : 'Ramsay', 'eq' : 2, 'threshold' : 3.021, \
			'dof' : 'Elbow'}, \
		'Mass' : 163.8,\
		'Actual No' : 320,\
		'Corrected No' : 292.6,\
		'Relative Abundance' : 1.1,\
		'Optimal Muscle Length' : 116,\
	    'Group' : 'flexor'}
	TRI_Coefficients = {\
		'Shoulder' : {\
			'MA' : [-25.40,0,0,0,0,0], \
			'src' : 'Pigeon', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : [-24.5454,-8.8691,9.3509,-1.7518,0],\
			'src' : 'Ramsay', 'eq' : 1, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : (94.2+138.4+92.5), \
		'Actual No' : (200+222+98),\
		'Corrected No' : (223.7+269.6+221.8),\
		'Relative Abundance' : (0.89+0.82+0.44)/3,\
		'Optimal Muscle Length' : 134,\
	    'Group' : 'extensor'}
	BRA_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : [16.1991,-16.1463,24.5512,-6.3335,0],\
			'src' : 'Ramsay', 'eq' : 1, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 141,\
		'Actual No' : 256,\
		'Corrected No' : 272.1,\
		'Relative Abundance' : 0.94,\
		'Optimal Muscle Length' : 86,\
	    'Group' : 'flexor'}
	BRD_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0, \
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : 	[15.2564,-11.8355,2.8129,-5.7781,44.8143,0,2.9032,0,0,-13.4956,0,-0.3940,0,0,0,0],\
			'src' : 'Ramsay', 'eq' : 2, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 64.7,\
		'Actual No' : 70,\
		'Corrected No' : 190.2,\
		'Relative Abundance' : 0.37,\
		'Optimal Muscle Length' : 173,\
	    'Group' : 'flexor'}
	PRO_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : 	[11.0405,-1.0079,0.3933,-10.4824,-12.1639,-0.4369,36.9174,3.5232,-10.4223,21.2604,-37.2444,10.2666,-11.0060,14.5974,-3.9919,1.7526,-2.0089,0.5460],\
			'src' : 'Ramsay', 'eq' : 3,'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 38.8, \
		'Actual No' : 187.6, \
		'Corrected No' : 185.5, \
		'Relative Abundance' : 1.3,\
		'Optimal Muscle Length' : 49,\
	    'Group' : 'flexor'}
	FCR_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0, \
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : Pigeon_coeff_conversion([0.9351,0.5375,-0.3627,0,0,0]),\
			'src' : 'Pigeon', 'eq' : None, 'threshold' : 2.86, \
			'dof' : 'Elbow'}, \
		'Mass' : 28.7, \
		'Actual No' : 129, \
		'Corrected No' : 125.7, \
		'Relative Abundance' : 1.0,\
		'Optimal Muscle Length' : 63,\
	    'Group' : 'flexor'}
	ECRB_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : [-11.256,17.8548,1.6398,-0.5073,-2.8827,0,-0.0942,0,0,0,0,0,0,0,0,0],\
			'src' : 'Ramsay', 'eq' : 2, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 32.1, \
		'Actual No' : 102, \
		'Corrected No' : 132.7, \
		'Relative Abundance' : 0.77,\
		'Optimal Muscle Length' : 59,\
	    'Group' : 'extensor'}
	ECRL_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : Pigeon_coeff_conversion([4.7304,1.2590,4.4347,-3.0229,0,0]),\
			'src' : 'Pigeon', 'eq' : None, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 44.3, \
		'Actual No' : 74, \
		'Corrected No' : 155.2, \
		'Relative Abundance' : 0.48,\
		'Optimal Muscle Length' : 81,\
	    'Group' : 'extensor'}
	FCU_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : 5,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 36.5,\
		'Actual No' : 175,\
		'Corrected No' : 141.2,\
		'Relative Abundance' : 1.2,\
		'Optimal Muscle Length' : 51,\
	    'Group' : 'flexor'}
	FDS_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : 5,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 95.2,\
		'Actual No' : 356,\
		'Corrected No' : 224.9,\
		'Relative Abundance' : 1.6,\
		'Optimal Muscle Length' : 84,\
	    'Group' : 'flexor'}
	PL_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : 10,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : None, \
		'Actual No' : None, \
		'Corrected No' : None, \
		'Relative Abundance' : None,\
		'Optimal Muscle Length' : 64,\
	    'Group' : 'flexor'}
	ECU_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : Pigeon_coeff_conversion([-2.1826,-1.7386,1.1491,0,0,0]),\
			'src' : 'Pigeon', 'eq' : None, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 25.2,\
		'Actual No' : 157,\
		'Corrected No' : 118,\
		'Relative Abundance' : 1.3,\
		'Optimal Muscle Length' : 62,\
	    'Group' : 'extensor'}
	EDM_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0, \
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : -5, \
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 6.2, \
		'Actual No' : 53, \
		'Corrected No' : 59.8, \
		'Relative Abundance' : 0.89,\
		'Optimal Muscle Length' : 68,\
	    'Group' : 'extensor'}
	EDC_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0, \
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : -5,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 42.8, \
		'Actual No' : 219, \
		'Corrected No' : 152.6, \
		'Relative Abundance' : 1.4,\
		'Optimal Muscle Length' : 70,\
	    'Group' : 'extensor'}
	AN_Coefficients = {\
		'Shoulder' : {\
			'MA' : 0,\
			'src' : 'Est', 'eq' : None, 'threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA' : Pigeon_coeff_conversion([-5.3450,-2.2841e-1,8.4297e-3,-14.329e-5,10.448e-7,-2.736e-9]),\
			'src' : 'Pigeon', 'eq' : None, 'threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : None, \
		'Actual No' : None, \
		'Corrected No' : None, \
		'Relative Abundance' : None,\
		'Optimal Muscle Length' : None,\
	    'Group' : 'extensor'}

	global AllCoefficients
	# AllCoefficients = [DELTa_Coefficients, CB_Coefficients, DELTp_Coefficients, BIC_Coefficients, \
	# 					TRI_Coefficients, BRA_Coefficients, BRD_Coefficients, PRO_Coefficients, \
	# 					FCR_Coefficients, ECRB_Coefficients, ECRL_Coefficients, FCU_Coefficients, \
	# 					FDS_Coefficients, PL_Coefficients, ECU_Coefficients, EDM_Coefficients, EDC_Coefficients]
	AllCoefficients = {'DELTa' : DELTa_Coefficients, 'CB' : CB_Coefficients, \
						'DELTp' : DELTp_Coefficients, 'BIC' : BIC_Coefficients, \
						'TRI' : TRI_Coefficients, 'BRA' : BRA_Coefficients, \
						'BRD' : BRD_Coefficients, 'PRO' : PRO_Coefficients, \
						'FCR' : FCR_Coefficients, 'ECRB' : ECRB_Coefficients, \
						'ECRL' : ECRL_Coefficients, 'FCU' : FCU_Coefficients, \
						'FDS' : FDS_Coefficients, 'ECU' : ECU_Coefficients, \
						'EDM' : EDM_Coefficients, 'EDC' : EDC_Coefficients}
	MuscleList = AllCoefficients.keys()
	# AllCoefficients = [DELTa_Coefficients, CB_Coefficients, DELTp_Coefficients, BIC_Coefficients, \
	# 					TRI_Coefficients, BRA_Coefficients, BRD_Coefficients, PRO_Coefficients, \
	# 					FCR_Coefficients, ECRB_Coefficients, ECRL_Coefficients, FCU_Coefficients, \
	# 					FDS_Coefficients, ECU_Coefficients, EDM_Coefficients, EDC_Coefficients]

	global RMatrix_Transpose
	RMatrix_Transpose = sp.Matrix([[MA_function(AllCoefficients[muscle][dof]) for dof in ['Shoulder','Elbow']] for muscle in MuscleList])
def return_MA_matrix(A1,A2):
	"""
	Notes:

	The angle of pronation is fixed at pi/2 for this reaching paradigm.
	The angle of radial/ulnar deviation is set to zero for the fixed wrist apparatus of the model paradigm.
	These functions have been verified to match the previous posture dependent MAs from Ramsey (2010) - SEE ERRATUM
	"""
	import numpy as np
	from numpy import pi
	import sympy as sp
	global q1,q2,q_PS,RMatrix_Transpose
	MomentArmMatrix = RMatrix_Transpose.subs([(q1,A1),(q2,A2),(q_PS,pi/2)]).T
	return(MomentArmMatrix)
def diagonal(X):
	import sympy as sp
	import numpy as np
	return(sp.Matrix(np.diagflat(X)))
def calculate_muscle_velocities(A1,A2,Ȧ1,Ȧ2):
	"""
	v = (-dR/dϑ)⋅(dϑ/dt)⋅ϑ + (-R)⋅(dϑ/dt)
	(dsⱼ/dϑ₁)⋅(dϑ₁/dt) = √((dϑ₁/dt)²⋅(dr₁ⱼ/dϑ₁)² + (dϑ₁/dt)²⋅(r₁ⱼ)²)
	(dsⱼ/dt) = sum([√(el) for el in diag(dϑ/dt)(J²_{Rⱼ})(dϑ/dt) + diag(Rⱼ)(diag(dϑ/dt)²)(dϑ₁/dt)])

	Need to hard code whether the rotation causes lengthening or shortening. This is determined by the sign of -r₁ⱼ⋅dϑ₁/dt for each d.o.f.
	"""
	import numpy as np
	import sympy as sp
	global q1,q2,q_PS,AllCoefficients,RMatrix_Transpose
	OptimalMuscleLength = np.array([AllCoefficients[key]['Optimal Muscle Length'] for key in AllCoefficients.keys()])
	R = return_MA_matrix(A1,A2)
	# J²_{Rⱼ} --> The squared jacobian of each column of the R matrix w.r.t. {ϑ₁,...,ϑₙ}
	SquaredJacobians = [RMatrix_Transpose[i,:].jacobian([q1,q2]).subs([(q1,A1),(q2,A2),(q_PS,np.pi/2)])**2 for i in range(np.shape(RMatrix_Transpose)[0])]
	sign_convention = lambda j: -np.sign(np.multiply([Ȧ1,Ȧ2],R[:,j].T))
	MuscleVelocity = [float(sign_convention(j)\
		*np.matrix([el**0.5 for el in \
			diagonal([Ȧ1,Ȧ2])*SquaredJacobians[j]*np.matrix([[Ȧ1],[Ȧ2]]) \
				+ diagonal(R[:,j])*(diagonal([Ȧ1,Ȧ2])**2)*R[:,j]]).T) \
					for j in range(np.shape(R)[1])]
	# dR = np.concatenate([sp.diff(RMatrix_Transpose[:,0],q1).subs(q1,A1),sp.diff(RMatrix_Transpose[:,1],q2).subs([(q2,A2),(q_PS,np.pi/2)])],axis=1)
	# MuscleVelocity = -R.T*(np.matrix([Ȧ1,Ȧ2]).T)#-dR*(np.matrix(np.multiply([A1,A2],[Ȧ1,Ȧ2])).T)
	NormalizedMuscleVelocity = np.array([[MuscleVelocity[i]/OptimalMuscleLength[i] for i in range(len(MuscleVelocity))]])
	return(NormalizedMuscleVelocity)
def calculate_weighted_unscaled_potential_torque_variations_REACHING_TASK(A1,A2,NormalizedMuscleVelocity):
	"""
	v = (-dR/dϑ)⋅(dϑ/dt)⋅ϑ + (-R)⋅(dϑ/dt)
	(dsⱼ/dϑ₁)⋅(dϑ₁/dt) = √((dϑ₁/dt)²⋅(dr₁ⱼ/dϑ₁)² + (dϑ₁/dt)²⋅(r₁ⱼ)²)
	(dsⱼ/dt) = sum([√(el) for el in diag(dϑ/dt)(J²_{Rⱼ})(dϑ/dt) + diag(Rⱼ)(diag(dϑ/dt)²)(dϑ₁/dt)])

	Need to hard code whether the rotation causes lengthening or shortening. This is determined by the sign of -r₁ⱼ⋅dϑ₁/dt for each d.o.f.
	"""
	import numpy as np
	import sympy as sp
	global AllCoefficients
	OptimalMuscleLength = np.array([AllCoefficients[key]['Optimal Muscle Length'] for key in AllCoefficients.keys()])
	R = return_MA_matrix(A1,A2)
	# # J²_{Rⱼ} --> The squared jacobian of each column of the R matrix w.r.t. {ϑ₁,...,ϑₙ}
	# SquaredJacobians = [RMatrix_Transpose[i,:].jacobian([q1,q2]).subs([(q1,A1),(q2,A2),(q_PS,np.pi/2)])**2 for i in range(np.shape(RMatrix_Transpose)[0])]
	# sign_convention = lambda j: -np.sign(np.multiply([Ȧ1,Ȧ2],R[:,j].T))
	# MuscleVelocity = [float(sign_convention(j)\
	# 	*np.matrix([el**0.5 for el in \
	# 		diagonal([Ȧ1,Ȧ2])*SquaredJacobians[j]*np.matrix([[Ȧ1],[Ȧ2]]) \
	# 			+ diagonal(R[:,j])*(diagonal([Ȧ1,Ȧ2])**2)*R[:,j]]).T) \
	# 				for j in range(np.shape(R)[1])]
	MuscleList = list(AllCoefficients.keys())

	WeightedNormalizedMuscleVelocity_shoulder = []
	WeightedNormalizedMuscleVelocity_elbow = []

	for i in range(len(AllCoefficients)):
		WeightedNormalizedMuscleVelocity_shoulder.append(NormalizedMuscleVelocity[i]*(abs(R.T[i,0] )/1000)*AllCoefficients[MuscleList[i]]['Corrected No'])

		WeightedNormalizedMuscleVelocity_elbow.append(NormalizedMuscleVelocity[i]*(abs(R.T[i,1] )/1000)*AllCoefficients[MuscleList[i]]['Corrected No'])

	return(np.array([WeightedNormalizedMuscleVelocity_shoulder]),np.array([WeightedNormalizedMuscleVelocity_elbow]))
def calculate_weighted_muscle_velocities_REACHING_TASK(A1,A2,Ȧ1,Ȧ2):
	"""
	v = (-dR/dϑ)⋅(dϑ/dt)⋅ϑ + (-R)⋅(dϑ/dt)
	(dsⱼ/dϑ₁)⋅(dϑ₁/dt) = √((dϑ₁/dt)²⋅(dr₁ⱼ/dϑ₁)² + (dϑ₁/dt)²⋅(r₁ⱼ)²)
	(dsⱼ/dt) = sum([√(el) for el in diag(dϑ/dt)(J²_{Rⱼ})(dϑ/dt) + diag(Rⱼ)(diag(dϑ/dt)²)(dϑ₁/dt)])

	Need to hard code whether the rotation causes lengthening or shortening. This is determined by the sign of -r₁ⱼ⋅dϑ₁/dt for each d.o.f.
	"""
	import numpy as np
	import sympy as sp
	global q1,q2,q_PS,AllCoefficients,RMatrix_Transpose
	OptimalMuscleLength = np.array([AllCoefficients[key]['Optimal Muscle Length'] for key in AllCoefficients.keys()])
	R = return_MA_matrix(A1,A2)
	# J²_{Rⱼ} --> The squared jacobian of each column of the R matrix w.r.t. {ϑ₁,...,ϑₙ}
	SquaredJacobians = [RMatrix_Transpose[i,:].jacobian([q1,q2]).subs([(q1,A1),(q2,A2),(q_PS,np.pi/2)])**2 for i in range(np.shape(RMatrix_Transpose)[0])]
	sign_convention = lambda j: -np.sign(np.multiply([Ȧ1,Ȧ2],R[:,j].T))
	MuscleVelocity = [float(sign_convention(j)\
		*np.matrix([el**0.5 for el in \
			diagonal([Ȧ1,Ȧ2])*SquaredJacobians[j]*np.matrix([[Ȧ1],[Ȧ2]]) \
				+ diagonal(R[:,j])*(diagonal([Ȧ1,Ȧ2])**2)*R[:,j]]).T) \
					for j in range(np.shape(R)[1])]
	# import numpy as np
	# import sympy as sp
	# import ipdb
	# global q1,q2,q_PS,AllCoefficients
	# R = return_MA_matrix(A1,A2)
	# #dR = np.concatenate([sp.diff(RMatrix_Transpose[:,0],q1).subs(q1,A1),sp.diff(RMatrix_Transpose[:,1],q2).subs([(q2,A2),(q_PS,np.pi/2)])],axis=1)
	# MuscleVelocity = -R.T*(np.matrix([Ȧ1,Ȧ2]).T)#-dR*(np.matrix(np.multiply([A1,A2],[Ȧ1,Ȧ2])).T)#
	MuscleList = list(AllCoefficients.keys())
	# WeightedNormalizedMuscleVelocity = np.array([[float(MuscleVelocity[i,0])/AllCoefficients[list(MuscleList)[i]]['Optimal Muscle Length'] for i in range(len(MuscleVelocity))]])
	WeightedNormalizedMuscleVelocity = []
	for i in range(len(AllCoefficients)):
		if MuscleList[i] not in ['BIC','TRI']:
			WeightedNormalizedMuscleVelocity.append(MuscleVelocity[i]*abs(sum(R.T[i,:] )/1000)*AllCoefficients[MuscleList[i]]['Corrected No']/AllCoefficients[MuscleList[i]]['Optimal Muscle Length'])
		else:
			WeightedNormalizedMuscleVelocity.append(MuscleVelocity[i]*abs((sum(R.T[i,:] )/2)/1000)*AllCoefficients[MuscleList[i]]['Corrected No']/AllCoefficients[MuscleList[i]]['Optimal Muscle Length'])
	return(np.array([WeightedNormalizedMuscleVelocity]))
def eccentric_velocities(NormalizedMuscleVelocity):
	import numpy as np
	PositiveMuscleVelocities = []
	for i in range(np.shape(NormalizedMuscleVelocity)[1]):
		PositiveEntries = np.multiply(np.array([[1]*np.shape(NormalizedMuscleVelocity)[0]]),\
										NormalizedMuscleVelocity.T[i,:]>=0)
		PositiveMuscleVelocities.append(np.multiply(NormalizedMuscleVelocity.T[i,:],\
														PositiveEntries))
	return(PositiveMuscleVelocities)
def concentric_velocities(NormalizedMuscleVelocity):
	import numpy as np
	NegativeMuscleVelocities = [NormalizedMuscleVelocity.T[i,:][NormalizedMuscleVelocity.T[i,:]<=0] for i in range(np.shape(NormalizedMuscleVelocity)[1])]
	return(NegativeMuscleVelocities)
def cost_function(X,type="avg"):
	assert type in ['avg','sos','l2norm','l1norm'], "type must be either 'avg','sos','l2norm', or 'l1norm'"
	if type == 'avg' :
		cost = abs(sum(X)/len(X))
	elif type == 'sos' :
		cost = sum([el**2 for el in X])
	elif type == 'l1norm' :
		cost = abs(sum(X))
	elif type == 'l2norm' :
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
import time
from matplotlib.backends.backend_pdf import PdfPages

	# Save Output Figures?

SaveOutputFigures = False

	# DescriptiveTitle should be something to identify the trial either by reaching location, (i.e., Left, Right, Center, or Sideways) or by what has changed in the most current iteration (e.g., CB_Ramsay, DELTa_Est, etc.). Spaces will be replaced by '_' symbols for the filename but kept for figure titles.

ReachType = ['Sideways','Center','Left','Right'][1]
DescriptiveTitle = ReachType + ' Reach'

	# Open PDF files to write plots to

if SaveOutputFigures == True:
	pdf_forward = PdfPages(DescriptiveTitle.replace(' ','_') + '_Forward.pdf')
	pdf_reverse = PdfPages(DescriptiveTitle.replace(' ','_') + '_Reverse.pdf')
	pdf_bar = PdfPages(DescriptiveTitle.replace(' ','_') + '_bar.pdf')

	# Set reaching duration and sampling period

t_end = 1
dt = 0.001
t = np.arange(0,1+dt,dt)*t_end

	# Define the Model and establish flexor/extensor list

global_R_matrix()
n_muscles=len(AllCoefficients)
n_extensors = sum([AllCoefficients[key]['Group']=='extensor' for key in AllCoefficients.keys()])
n_flexors = sum([AllCoefficients[key]['Group']=='flexor' for key in AllCoefficients.keys()])
flexor_cmap=plt.get_cmap('autumn')
extensor_cmap = plt.get_cmap('YlGnBu')
flexor_colors = iter(flexor_cmap(np.linspace(0,1,n_flexors)))
extensor_colors = iter(list(reversed(extensor_cmap(np.linspace(0,1,n_flexors)))))
flexor_colors_list = []
extensor_colors_list = []
FlexorOrderedMuscleList = []
ExtensorOrderedMuscleList = []
for muscle in AllCoefficients:
   if AllCoefficients[muscle]['Group'] == 'flexor':
        flexor_colors_list.append(next(flexor_colors))
        FlexorOrderedMuscleList.append(muscle)
   else:
        extensor_colors_list.append(next(extensor_colors))
        ExtensorOrderedMuscleList.append(muscle)
OrderedColorsList = flexor_colors_list + list(reversed(extensor_colors_list))
OrderedMuscleList = FlexorOrderedMuscleList + list(reversed(ExtensorOrderedMuscleList))
OrderNumber = [list(AllCoefficients.keys()).index(el) for el in OrderedMuscleList]

	# Define Initial and Final positions

Xi,Xf = set_initial_and_final_positions(ReachType)

	# Forward Direction

reaching_task(Xi=Xi, Xf=Xf,dt=dt, t_end=t_end)

	# Initialize Arrays before loop

# calculate_muscle_lengths()
# MuscleLengths_Forward = MuscleLengths
# WeightedNormalizedMuscleVelocity_Forward = calculate_weighted_muscle_velocities_REACHING_TASK(A1[0],A2[0],Ȧ1[0],Ȧ2[0])
NormalizedMuscleVelocity_Forward = calculate_muscle_velocities(A1[0],A2[0],Ȧ1[0],Ȧ2[0])
WeightedPotentialTorqueVariation_shoulder_Forward,\
WeightedPotentialTorqueVariation_elbow_Forward = \
	calculate_weighted_unscaled_potential_torque_variations_REACHING_TASK(A1[0],A2[0],NormalizedMuscleVelocity_Forward[0])

	# Iterate over the rest of the movement.

# global A1_dev, A2_dev
# A1_dev,A2_dev = [],[]

StartTime = time.time()
statusbar(0,len(A1),Title = 'Forward Vm',StartTime=StartTime)
for i in range(1,len(A1)):
	statusbar(i,len(A1),Title = 'Forward Vm',StartTime=StartTime)
	# WeightedNormalizedMuscleVelocity_Forward = np.concatenate((WeightedNormalizedMuscleVelocity_Forward,calculate_weighted_muscle_velocities_REACHING_TASK(A1[i],A2[i],Ȧ1[i],Ȧ2[i])),axis=0)
	NormalizedMuscleVelocity_Forward = np.concatenate( \
				(NormalizedMuscleVelocity_Forward,\
					calculate_muscle_velocities(A1[i],A2[i],Ȧ1[i],Ȧ2[i])),\
						axis=0)
	WeightedPotentialTorqueVariation_temp = \
		calculate_weighted_unscaled_potential_torque_variations_REACHING_TASK(A1[i],A2[i],\
			NormalizedMuscleVelocity_Forward[i])
	WeightedPotentialTorqueVariation_shoulder_Forward = \
			np.concatenate((WeightedPotentialTorqueVariation_shoulder_Forward,\
						WeightedPotentialTorqueVariation_temp[0]),axis=0)
	WeightedPotentialTorqueVariation_elbow_Forward = \
			np.concatenate((WeightedPotentialTorqueVariation_elbow_Forward,\
						WeightedPotentialTorqueVariation_temp[1]),axis=0)

print('\n')

	# Adding the weighted values regains the Weighted Normalized Muscle Velocities. But it is important to average the BIARTICULATING MUSCLES as shown below when considering the average afferent response.

WeightedNormalizedMuscleVelocity_Forward = WeightedPotentialTorqueVariation_elbow_Forward + WeightedPotentialTorqueVariation_shoulder_Forward

	# Correct for BIC and TRI by averaging the contributions from the MAs.

WeightedNormalizedMuscleVelocity_Forward[:,3:5] = WeightedNormalizedMuscleVelocity_Forward[:,3:5]/2

	# Calculate only the lengthening Torque Variation Contributions for both the shoulder and the elbow.

EccentricTorqueVariations_Shoulder =\
 	eccentric_velocities(WeightedPotentialTorqueVariation_shoulder_Forward)
EccentricTorqueVariations_Elbow =\
 	eccentric_velocities(WeightedPotentialTorqueVariation_elbow_Forward)

	# Concatenate, Sum, and Scale the total lengthening Torque Variations for both Shoulder and Elbow.

ScalingFactor = 2000
TotalEccentricTorqueVariations_Shoulder =\
	ScalingFactor*np.concatenate(EccentricTorqueVariations_Shoulder,axis=0)
TotalEccentricTorqueVariations_Elbow =\
	ScalingFactor*np.concatenate(EccentricTorqueVariations_Elbow,axis=0)

	# Calculate the time series of potential endpoint variation.
MuscleList = list(AllCoefficients.keys())
TotalEndpointPotentialVariation = []

for j in range(len(MuscleList)):
	TotalEndpointPotentialVariation.append( \
		[calculate_potential_variability(\
			i,\
			TotalEccentricTorqueVariations_Shoulder[j,:],\
			TotalEccentricTorqueVariations_Elbow[j,:],\
			EOM='Uno',\
			dt=dt,\
			Group = AllCoefficients[MuscleList[j]]['Group']) \
				for i in range(len(A1))])

	# Plot the resulting potential variation

[plt.plot(TotalEndpointPotentialVariation[i]) for i in range(16)]

	# Plot Afferent-Weighted Muscle Velocity (Forward)

fig1e = plt.figure()
[plt.plot(t,WeightedNormalizedMuscleVelocity_Forward[:,i]) for i in OrderNumber]
ax1e = plt.gca()
ax1e.set_xlim(0,t_end*(1.3))
ax1e.set_ylim(-12,12)
if t_end!=1:
	ax1e.set_xlabel('Time (s)')
else:
	ax1e.set_xlabel('Normalized Time')
ax1e.set_ylabel('Afferent-Weighted $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
ax1e.set_title(DescriptiveTitle+'\nAfferent-Weighted Normalized Muscle Velocity')
[j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax1e.lines)]
ax1e.legend(OrderedMuscleList)

plt.show()
#
#
# 	# Plot Normalized Muscle Velocity (Forward)
#
# fig1a = plt.figure()
# [plt.plot(t,NormalizedMuscleVelocity_Forward[:,i]) for i in OrderNumber]
# ax1a = plt.gca()
# ax1a.set_xlim(0,t_end*(1.3))
# ax1a.set_ylim(-1.5,1.5)
# ax1a.set_title(DescriptiveTitle+'\nNormalized Muscle Velocity')
# if t_end!=1:
# 	ax1a.set_xlabel('Time (s)')
# else:
# 	ax1a.set_xlabel('Normalized Time')
# ax1a.set_ylabel('Normalized Muscle Velocity\nConcentric $\longleftrightarrow$ Eccentric')
# [j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax1a.lines)]
# ax1a.legend(OrderedMuscleList)
#
#
# 	# Plot Joint Angles (Forward)
#
# fig1b = plt.figure()
# plt.plot(t,A1,'g')
# plt.plot(t,A2,'g:')
# ax1b = plt.gca()
# ax1b.set_xlim(0,t_end*(1.3))
# ax1b.set_title(DescriptiveTitle+'\nAngle vs. Time')
# if t_end!=1:
# 	ax1b.set_xlabel('Time (s)')
# else:
# 	ax1b.set_xlabel('Normalized Time')
# ax1b.set_ylabel('Joint Angles (in radians)')
# ax1b.legend(["Shoulder","Elbow"])
#
#
# 	# Plot Configuration Space (Forward)
#
# fig1c = plt.figure()
# plt.plot(A1,A2,'g')
# ax1c = plt.gca()
# ax1c.set_title(DescriptiveTitle+'\nConfiguration Space')
# ax1c.set_ylabel("Shouler Angle (in radians)")
# ax1c.set_xlabel("Elbow Angle (in radians)")
#
#
# 	# Calculate Costs from Afferent-Weighted Muscle Velocities
#
# EccentricCost_Forward = eccentric_cost(WeightedNormalizedMuscleVelocity_Forward,t_end=t_end,dt = dt,type='l1norm')
# ConcentricCost_Forward = concentric_cost(WeightedNormalizedMuscleVelocity_Forward,t_end=t_end,dt = dt,type='l1norm')
#
# # # Plot Muscle Lengths (Forward)
# #
# # fig1d = plt.figure()
# # [plt.plot(t,MuscleLengths_Forward[i]) for i in OrderNumber]
# # ax1d = plt.gca()
# # ax1d.set_xlim(0,t_end*(1.3))
# # ax1d.set_title(DescriptiveTitle+'\nMuscle Lengths')
# # if t_end!=1:
# # 	ax1d.set_xlabel('Time (s)')
# # else:
# # 	ax1d.set_xlabel('Normalized Time')
# # ax1d.set_ylabel('Muscle Lengths (in mm)')
# # [j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax1d.lines)]
# # ax1d.legend(OrderedMuscleList)
#
#
# 	# Reverse Direction
#
# reaching_task(Xi=Xf, Xf=Xi,dt=dt, t_end=t_end)
# # calculate_muscle_lengths()
# # MuscleLengths_Reverse = MuscleLengths
# WeightedNormalizedMuscleVelocity_Reverse = calculate_weighted_muscle_velocities_REACHING_TASK(A1[0],A2[0],Ȧ1[0],Ȧ2[0])
# NormalizedMuscleVelocity_Reverse = calculate_muscle_velocities(A1[0],A2[0],Ȧ1[0],Ȧ2[0])
# StartTime = time.time()
# statusbar(0,len(A1),Title = 'Reverse Vm',StartTime=StartTime)
# for i in range(1,len(A1)):
#     statusbar(i,len(A1),Title = 'Reverse Vm',StartTime=StartTime)
#     WeightedNormalizedMuscleVelocity_Reverse = np.concatenate((WeightedNormalizedMuscleVelocity_Reverse,calculate_weighted_muscle_velocities_REACHING_TASK(A1[i],A2[i],Ȧ1[i],Ȧ2[i])),axis=0)
#     NormalizedMuscleVelocity_Reverse = np.concatenate((NormalizedMuscleVelocity_Reverse,calculate_muscle_velocities(A1[i],A2[i],Ȧ1[i],Ȧ2[i])),axis=0)
# print('\n')
#
# 	# Plot Afferent-Weighted Muscle Velocity (Reverse)
#
# fig2e = plt.figure()
# [plt.plot(t,WeightedNormalizedMuscleVelocity_Reverse[:,i]) for i in OrderNumber]
# ax2e = plt.gca()
# ax2e.set_xlim(0,t_end*(1.3))
# ax2e.set_ylim(-12,12)
# ax2e.set_title(DescriptiveTitle+'\nAfferent-Weighted Normalized Muscle Velocity')
# if t_end!=1:
# 	ax2e.set_xlabel('Time (s)')
# else:
# 	ax2e.set_xlabel('Normalized Time')
# ax2e.set_ylabel('Afferent-Weighted $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
# [j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax2e.lines)]
# ax2e.legend(OrderedMuscleList)
#
#
# 	# Plot Normalized Muscle Velocity (Reverse)
#
# fig2a = plt.figure()
# [plt.plot(t,NormalizedMuscleVelocity_Reverse[:,i]) for i in OrderNumber]
# ax2a = plt.gca()
# ax2a.set_xlim(0,t_end*(1.3))
# ax2a.set_ylim(-1.5,1.5)
# ax2a.set_title(DescriptiveTitle+'\nNormalized Muscle Velocity')
# if t_end!=1:
# 	ax2a.set_xlabel('Time (s)')
# else:
# 	ax2a.set_xlabel('Normalized Time')
# ax2a.set_ylabel('Normalized Muscle Velocity\nConcentric $\longleftrightarrow$ Eccentric')
# [j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax2a.lines)]
# ax2a.legend(OrderedMuscleList)
#
#
# 	# Plot Joint Angles (Reverse)
#
# fig2b = plt.figure()
# plt.plot(t,A1,'b')
# plt.plot(t,A2,'b:')
# ax2b = plt.gca()
# ax2b.set_xlim(0,t_end*(1.3))
# ax2b.set_title(DescriptiveTitle+'\nAngle vs. Time')
# if t_end!=1:
# 	ax2b.set_xlabel('Time (s)')
# else:
# 	ax2b.set_xlabel('Normalized Time')
# ax2b.set_ylabel('Joint Angles (in radians)')
# ax2b.legend(["Shoulder","Elbow"])
#
#
# 	# Plot Configuration Space (Reverse)
#
# fig2c = plt.figure()
# plt.plot(A1,A2,'b')
# ax2c = plt.gca()
# ax2c.set_title(DescriptiveTitle+'\nConfiguration Space')
# ax2c.set_ylabel("Shouler Angle (in radians)")
# ax2c.set_xlabel("Elbow Angle (in radians)")
#
#
# 	# Calculate Costs from Afferent-Weighted Muscle Velocities
#
# EccentricCost_Reverse = eccentric_cost(WeightedNormalizedMuscleVelocity_Reverse,t_end=t_end,dt = dt,type='l1norm')
# ConcentricCost_Reverse = concentric_cost(WeightedNormalizedMuscleVelocity_Reverse,t_end=t_end,dt = dt,type='l1norm')
#
# # # Plot Muscle Lengths (Reverse)
# #
# # fig2d = plt.figure()
# # [plt.plot(t,MuscleLengths_Reverse[i]) for i in OrderNumber]
# # ax2d = plt.gca()
# # ax2d.set_xlim(0,t_end*(1.3))
# # ax2d.set_title(DescriptiveTitle+'\nMuscle Lengths')
# # if t_end!=1:
# # 	ax2d.set_xlabel('Time (s)')
# # else:
# # 	ax2d.set_xlabel('Normalized Time')
# # ax2d.set_ylabel('Muscle Lengths (in mm)')
# # [j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax2d.lines)]
# # ax2d.legend(OrderedMuscleList)
#
#
# # plt.figure()
# # plt.plot(np.arange(0,2*math.ceil(max([EccentricCost_Forward,ConcentricCost_Forward]))+1,1), np.arange(0,2*math.ceil(max([EccentricCost_Forward,ConcentricCost_Forward]))+1,1),'0.75',linestyle='--')
# # plt.scatter([EccentricCost_Forward],[ConcentricCost_Forward],facecolor ='g', edgecolor = 'k', marker='o',s=100)
# # plt.scatter([EccentricCost_Reverse],[ConcentricCost_Reverse],facecolor ='b', edgecolor = 'k', marker='o',s=100)
# # ax3 = plt.gca()
# # ax3.set_xlim(0,2*math.ceil(max([EccentricCost_Forward,EccentricCost_Reverse])))
# # ax3.set_ylim(0,2*math.ceil(max([EccentricCost_Forward,EccentricCost_Reverse])))
# # ax3.set_aspect('equal')
# # ax3.set_title('Eccentric vs. Concentric Cost\nFor Forward and Reverse Movments')
# # ax3.set_xlabel('Eccentric Cost (in $\hat{l}_{o}$)')
# # ax3.set_ylabel('Concentric Cost (in $\hat{l}_{o}$)')
# # ax3.legend(['Equal Cost Line','Forward Direction','Reversed Direction'])
#
#
# # Using Default Reaching Directions Listed Above...
# # CostValues ={'CenterForward' : np.array([[ 1.3220062 ,  2.10011198]]),\
# # 				'CenterReverse' : np.array([[ 2.10011198,  1.3220062 ]]),\
# # 					'LeftForward' : np.array([[ 1.24643599,  1.88616011]]),\
# # 						'LeftReverse' : np.array([[ 1.88616011,  1.24643599]]),\
# # 							'RightForward' : np.array([[ 1.1047386 ,  1.92105444]]),\
# # 								'RightReverse' : np.array([[ 1.92105444,  1.1047386 ]]),\
# # 									'SideForward' : np.array([[ 0.37931692,  0.53156831]]),\
# # 										'SideReverse' : np.array([[ 0.53156831,  0.37931692]])}
# # plt.figure()
# # MaximumCost = max([max(CostValues[key][0]) for key in CostValues.keys()])
# # plt.plot(np.arange(0,2*math.ceil(MaximumCost)+1,1), np.arange(0,2*math.ceil(MaximumCost)+1,1),'0.75',linestyle='--')
# # facecolor = iter(['k','#9f9f9f','b','#9f9f9f','r','#9f9f9f','g','#9f9f9f'])
# # edgecolor = iter(['k','k','k','b','k','r','k','g'])
# # [plt.scatter(CostValues[key][0,0],CostValues[key][0,1],c=next(facecolor),edgecolor=next(edgecolor),marker='o',s=75) for key in CostValues.keys()]
# # ax4 = plt.gca()
# # ax4.set_xlim(0,math.ceil(MaximumCost))
# # ax4.set_ylim(0,math.ceil(MaximumCost))
# # ax4.set_aspect('equal')
# # ax4.set_title('Eccentric vs. Concentric Cost\nFor Forward and Reverse Movments')
# # ax4.set_xlabel('Eccentric Cost (in $\hat{l}_{o}$)')
# # ax4.set_ylabel('Concentric Cost (in $\hat{l}_{o}$)')
# # # ax4.legend(['Equal Cost Line','Center (Forward)','Center (Reverse)','Left (Forward)','Left (Reverse)','Right (Forward)','Right (Reverse)','Side-to-Side (Forward)','Side-to-Side (Reverse)'],loc='upper right')
# #
#
# 	# Plot bar graph comparing the two directions
#
# fig2 = plt.figure()
# plt.bar(np.arange(2),[EccentricCost_Reverse,EccentricCost_Forward])
# ax4 = plt.gca()
# ax4.set_xticklabels(('Reverse','Forward'))
# ax4.set_ylim(0,18)
# ax4.set_yticks([0,9,18])
# ax4.set_yticklabels(['0','','18'])
# ax4.set_title(DescriptiveTitle + '\nEccentric Cost for Forward and Reverse Movements')
# ax4.set_ylabel('Sum of Afferent-Weighted Muscle Lengthening')
# ax4.set_xticks([0,1])
#
#
# 	# Need to write PDF files and close the PDF files/figures
# if SaveOutputFigures == True:
# 	pdf_forward.savefig(fig1e)
# 	pdf_forward.savefig(fig1a)
# 	pdf_forward.savefig(fig1b)
# 	pdf_forward.savefig(fig1c)
# 	# pdf_forward.savefig(fig1d)
# 	pdf_reverse.savefig(fig2e)
# 	pdf_reverse.savefig(fig2a)
# 	pdf_reverse.savefig(fig2b)
# 	pdf_reverse.savefig(fig2c)
# 	# pdf_reverse.savefig(fig2d)
# 	pdf_bar.savefig(fig2)
#
# 	pdf_forward.close()
# 	pdf_reverse.close()
# 	pdf_bar.close()
#
# 		# Close any open plots
#
# 	plt.close('all')
#
# plt.show()
