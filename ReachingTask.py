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
    Changes the values of the link lengths. New_L1 and New_L2 must be numbers.
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
def find_Ẋ_values(t,Xi,Xf):
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
    ẋ = (xf-xi)*(30*t**2 - 60*t**3 + 30*t**4)
    ẏ = (yf-yi)*(30*t**2 - 60*t**3 + 30*t**4)
    Ẋ = np.array(np.concatenate(([ẋ],[ẏ]),axis=0))
    return(Ẋ)
def find_Ẍ_values(t,Xi,Xf):
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
    ẍ = (xf-xi)*(60*t - 180*t**2 + 120*t**3)
    ÿ = (yf-yi)*(60*t - 180*t**2 + 120*t**3)
    Ẍ = np.array(np.concatenate(([ẍ],[ÿ]),axis=0))
    return(Ẍ)
def return_X_values(t,Xi,Xf):
    X = find_X_values(t,Xi,Xf)
    Ẋ = find_Ẋ_values(t,Xi,Xf)
    Ẍ = find_Ẍ_values(t,Xi,Xf)
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
def reaching_task(dt = 0.001, Xi = [0,0.25], Xf = [0.25,0.50]):
    import numpy as np
    set_link_lengths()
    create_angle_lists()
    t = np.arange(0,1+dt,dt)
    X,Ẋ,Ẍ = return_X_values(t,Xi,Xf)
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
def calculate_torque(EOM="Uno"):
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
    global τ1,τ2
    τ1,τ2 = [],[]
    for i in range(len(A1)):
        Ȧ = np.matrix([[Ȧ1[i]],[Ȧ2[i]]])
        Ä = np.matrix([[Ä1[i]],[Ä2[i]]])
        τ = M(A1[i],A2[i])*Ä + C(A1[i],A2[i],Ȧ1[i],A2[i])*Ȧ + B*Ȧ
        τ1.append(τ[0,0])
        τ2.append(τ[1,0])
def plot_resulting_torques():
    import matplotlib.pyplot as plt
    import numpy as np
    from math import cos,sin

    dt = 1/(len(A1)-1)
    t = np.arange(0,1+dt,dt)
    plt.figure()
    plt.plot(t,τ1,'r',t,τ2,'b')
    ax=plt.gca()
    ax.set_xlabel('Normalized Time')
    ax.set_ylabel('Joint Torques (N⋅m)')
    ax.legend(['Shoulder','Elbow'])
    # plt.show()
def constant_moment_arm_matrix(R=None):
    import numpy as np
    if R!= None:
        assert np.shape(R)==(2,4),"R must be a (2,4) Moment Arm Matrix"
        # assert np.shape(R)==(2,6),"R must be a (2,6) Moment Arm Matrix"
        r11,r21 = R[0,0], R[1,0] # m
        r12,r22 = R[0,1], R[1,1] # m
        r13,r23 = R[0,2], R[1,2] # m
        r14,r24 = R[0,3], R[1,3] # m
        # r15,r25 = R[0,4], R[1,4] # m
        # r16,r26 = R[0,5], R[1,5] # m
    else:
        r11,r21 = -0.055, 0 # m
        r12,r22 = 0.055, 0 # m
        r13,r23 = 0, -0.045 # m
        r14,r24 = 0, 0.045 # m
        # r15,r25 = -0.055, -0.045 # m
        # r16,r26 = 0.055, 0.045 # m
    R = np.matrix([ [r11, r12, r13, r14],\
                    [r21, r22, r23, r24]])
    # R = np.matrix([ [r11, r12, r13, r14, r15, r16],\
    #                 [r21, r22, r23, r24, r25, r26]])
    return(R)
def A_eq_matrix():
    import numpy as np
    R = constant_moment_arm_matrix()
    A_eq = R
    return(A_eq)
def A_ub_matrix():
    import numpy as np
    global f1max,f2max,f3max,f4max
    f1max,f2max,f3max,f4max = 1500,1500,1500,1500
    global cc_level
    cc_level = 0.5 # cc_level ∊ [0,1]
    cc_ub = 1+cc_level
    cc_lb = 1-cc_level
    CC_lb = np.array([[((cc_lb/2)-1)/f1max, cc_lb/(2*f2max), 0, 0],\
                    [0, 0, ((cc_lb/2)-1)/f3max, cc_lb/(2*f4max)]])
    CC_ub = np.array([[-((cc_ub/2)-1)/f1max, -cc_ub/(2*f2max), 0, 0],\
                    [0, 0, -((cc_ub/2)-1)/f3max, -cc_ub/(2*f4max)]])
    I = np.identity(4)
    A_ub = np.concatenate((CC_lb,CC_ub,I,-I),axis=0)
    return(A_ub)
def b_eq_vector(τ1,τ2):
    import numpy as np
    b_eq = np.array([τ1,τ2])
    return(b_eq)
def b_ub_vector():
    import numpy as np
    b_ub = np.array([0,0,\
                    0,0,\
                    f1max,f2max,f3max,f4max,\
                    0,0,0,0])
    return(b_ub)
def optimize_muscle_forces():
    import numpy as np
    from scipy.optimize import linprog
    f1,f2,f3,f4=[],[],[],[]
    c = np.array([1,1,1,1])
    R = constant_moment_arm_matrix()
    A_ub = A_ub_matrix()
    b_ub = b_ub_vector()
    for i in range(len(τ1)):
        A_eq = A_eq_matrix()
        b_eq = b_eq_vector(τ1[i],τ2[i])
        F = linprog(c,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq)
        f1.append(F['x'][0])
        f2.append(F['x'][1])
        f3.append(F['x'][2])
        f4.append(F['x'][3])
        statusbar(i,len(τ1),Title="Reaching Task")
    return(f1,f2,f3,f4)
def plot_muscle_forces(f1,f2,f3,f4):
    import numpy as np
    import matplotlib.pyplot as plt
    dt = 1/(len(f1)-1)
    t = np.arange(0,1+dt,dt)
    plt.figure()
    plt.plot(t,f1,'r',t,f2,'r--',t,f3,'b',t,f4,'b--')
    ax = plt.gca()
    ax.set_title(str(cc_level*100)+"% Co-Contraction \n(Solid Lines: Extensors, Dashed Lines: Flexors)")
    ax.set_xlabel('Normalized Time')
    ax.set_ylabel('Muscle Force (N)')
    ax.legend(["Muscle 1","Muscle 2", "Muscle 3", "Muscle 4"],loc="upper left",fancybox=True)
def show_plots():
    import matplotlib.pyplot as plt
    plt.show()
reaching_task()
calculate_torque(EOM="Zadravec")
f1,f2,f3,f4 = optimize_muscle_forces()
plot_resulting_kinematics()
plot_resulting_torques()
plot_muscle_forces(f1,f2,f3,f4)
show_plots()
