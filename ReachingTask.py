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
def set_link_lengths(New_L1=0.256,New_L2=0.315):
    """
    Changes the values of the link lengths. New_L1 and New_L2 must be numbers.
    """
    assert type(New_L1) == int or type(New_L1)==float, "New_L1 must be a number."
    assert type(New_L2) == int or type(New_L2)==float, "New_L2 must be a number."
    global L1, L2
    L1 = New_L1
    L2 = New_L2
def test_global_link_lengths():
    print('L1 = ' + str(L1))
    print('L2 = ' + str(L2))
def create_angle_lists():
    global Q1,Q2,ddt_Q1,ddt_Q2,d2dt2_Q1,d2dt2_Q2
    Q1=[]
    Q2=[]
    ddt_Q1 = []
    ddt_Q2 = []
    d2dt2_Q1 = []
    d2dt2_Q2 = []
def inverse_kinematics(X):
    """
    Takes in either a (2,) or (2,1) list/array with values for x and y (endpoint) and updates the lists Q1 and Q2 with the current angles (in radians) from the inverse kinematics.
    """
    import numpy as np
    from math import acos, sin, cos, atan
    from numpy import pi
    assert len(X)==2, "X must be either a (2,) or (2,1) list/array"
    x,y = X
    if np.shape(X) == (2,1): x,y = x[0],y[0]
    q2 = acos((x**2 + y**2 - L1**2 - L2**2)/(2*L1*L2))
    if x!= 0:
        q1 = atan(y/x) - atan((L2*sin(q2))/(L1+L2*cos(q2)))
    elif x == 0 and y > 0:
        q1 = pi/2 - atan((L2*sin(q2))/(L1+L2*cos(q2)))
    elif x == 0 and y < 0:
        q1 = -pi/2 - atan((L2*sin(q2))/(L1+L2*cos(q2)))
    if x < 0 and y > 0: q1 += pi
    # This accounts for atan(y/x) ∊ [-π,π] when then workspace requires y>0. For x,y < 0, this does not occur, and the range of motion prohibits x > 0 when y < 0.
    global Q1,Q2
    Q1.append(q1)
    Q2.append(q2)
def find_X_values(t,X_initial,X_final):
    """
    This takes t (a numpy.ndarray of normalized time ∈ [0,1]) and either a (2,) or (2,1) list/array with values for both initial and final x and y (endpoint). To avoid singularities, ||X[i]|| cannot be greater than L1 + L2.
    """
    import numpy as np
    assert len(X_final)==2, "X_final must be either a (2,) or (2,1) list/array"
    assert len(X_initial)==2, "X_initial must be either a (2,) or (2,1) list/array"
    xi,yi = X_initial
    xf,yf = X_final
    if np.shape(X_initial) == (2,1): xi,yi = xi[0],yi[0]
    if np.shape(X_final) == (2,1): xf,yf = xf[0],yf[0]
    x = xi + (xf-xi)*(10*t**3 - 15*t**4 + 6*t**5)
    y = yi + (yf-yi)*(10*t**3 - 15*t**4 + 6*t**5)
    X = np.array(np.concatenate(([x],[y]),axis=0))
    assert (sum(X**2)**0.5<L1+L2).all(), "Trajectory creates singularities."
    return(X)
def find_X_dot_values(t,X_initial,X_final):
    """
    This takes t (a numpy.ndarray of normalized time ∈ [0,1]) and either a (2,) or (2,1) list/array with values for both initial and final x and y (endpoint). To avoid singularities, ||X[i]|| cannot be greater than L1 + L2.
    """
    import numpy as np
    assert len(X_final)==2, "X_final must be either a (2,) or (2,1) list/array"
    assert len(X_initial)==2, "X_initial must be either a (2,) or (2,1) list/array"
    xi,yi = X_initial
    xf,yf = X_final
    if np.shape(X_initial) == (2,1): xi,yi = xi[0],yi[0]
    if np.shape(X_final) == (2,1): xf,yf = xf[0],yf[0]
    ddt_x = (xf-xi)*(30*t**2 - 60*t**3 + 30*t**4)
    ddt_y = (yf-yi)*(30*t**2 - 60*t**3 + 30*t**4)
    ddt_X = np.array(np.concatenate(([ddt_x],[ddt_y]),axis=0))
    return(ddt_X)
def find_X_d2dot_values(t,X_initial,X_final):
    """
    This takes t (a numpy.ndarray of normalized time ∈ [0,1]) and either a (2,) or (2,1) list/array with values for both initial and final x and y (endpoint).
    """
    import numpy as np
    assert len(X_final)==2, "X_final must be either a (2,) or (2,1) list/array"
    assert len(X_initial)==2, "X_initial must be either a (2,) or (2,1) list/array"
    xi,yi = X_initial
    xf,yf = X_final
    if np.shape(X_initial) == (2,1): xi,yi = xi[0],yi[0]
    if np.shape(X_final) == (2,1): xf,yf = xf[0],yf[0]
    d2dt2_x = (xf-xi)*(60*t - 180*t**2 + 120*t**3)
    d2dt2_y = (yf-yi)*(60*t - 180*t**2 + 120*t**3)
    d2dt2_X = np.array(np.concatenate(([d2dt2_x],[d2dt2_y]),axis=0))
    return(d2dt2_X)
def jacobian():
	import numpy as np
	from math import cos, sin
	J = np.matrix([[-L1*sin(Q1[-1])-L2*sin(Q1[-1]+Q2[-1]),-L2*sin(Q1[-1]+Q2[-1])],\
				[L1*cos(Q1[-1])+L2*cos(Q1[-1]+Q2[-1]),L2*cos(Q1[-1]+Q2[-1])]])
	return(J)
def update_angular_velocity(ddt_X):
    import numpy as np
    J = jacobian()
    ddt_Q = J**(-1)*np.array([ddt_X]).T
    global ddt_Q1,ddt_Q2
    ddt_Q1.append(ddt_Q[0,0])
    ddt_Q2.append(ddt_Q[1,0])
def update_angular_acceleration(ddt_X,d2dt2_X):
    from math import cos,sin
    import numpy as np
    assert len(ddt_X)==2, "ddt_X must be either a (2,) or (2,1) list/array"
    assert len(d2dt2_X)==2, "d2dt2_X must be either a (2,) or (2,1) list/array"
    ddt_x,ddt_y = ddt_X
    d2dt2_x,d2dt2_y = d2dt2_X
    if np.shape(ddt_X) == (2,1): ddt_x,ddt_y = ddt_x[0],ddt_y[0]
    if np.shape(d2dt2_X) == (2,1): d2dt2_x,d2dt2_y = d2dt2_x[0],d2dt2_y[0]
    q1,q2 = Q1[-1],Q2[-1]
    w1,w2 = ddt_Q1[-1],ddt_Q2[-1]
    a1 = ( ( (w1+w2)*(ddt_y*cos(q1+q2)-ddt_x*sin(q1+q2)) \
            + d2dt2_x*cos(q1+q2) \
            + d2dt2_y*sin(q1+q2) ) * L1*sin(q2) \
            - \
            (ddt_x*cos(q1+q2) + ddt_y*sin(q1+q2)) * L1*cos(q2)*w2 ) \
            / ((L1*sin(q2))**2)
    a2 = ( ( \
            ((w1+w2)*L2*sin(q1+q2) + w1*L1*sin(q1))*ddt_x \
            + (-L2*cos(q1+q2) - L1*cos(q1))*d2dt2_x \
            + (-(w1+w2)*L2*cos(q1+q2) - w1*L1*cos(q1))*ddt_y \
            + (-L2*sin(q1+q2)-L1*sin(q1))*d2dt2_y ) \
            * (L1*L2*sin(q2)) \
            - \
            ( (-L2*cos(q1+q2) - L1*cos(q1))*ddt_x \
            + (-L2*sin(q1+q2) - L1*sin(q1))*ddt_y ) * (L1*L2*cos(q2)*w2) \
            ) \
            / (L1*L2*sin(q2)**2)
    global d2dt2_Q1,d2dt2_Q2
    d2dt2_Q1.append(a1)
    d2dt2_Q2.append(a2)
def update_angle_lists(X,ddt_X,d2dt2_X):
    import numpy as np
    for i in range(np.shape(X)[1]):
        inverse_kinematics(X[:,i])
        update_angular_velocity(ddt_X[:,i])
        update_angular_acceleration(ddt_X[:,i],d2dt2_X[:,i])
def reaching_task(dt = 0.001, X_initial = [0,0.25], X_final = [0.25,0.50]):
    import numpy as np
    set_link_lengths()
    create_angle_lists()
    t = np.arange(0,1+dt,dt)
    X = find_X_values(t,X_initial,X_final)
    ddt_X = find_X_dot_values(t,X_initial,X_final)
    d2dt2_X = find_X_d2dot_values(t,X_initial,X_final)
    update_angle_lists(X,ddt_X,d2dt2_X)
def plot_resulting_kinematics():
    import matplotlib.pyplot as plt
    import numpy as np
    from math import cos,sin

    dt = 1/(len(Q1)-1)
    t = np.arange(0,1+dt,dt)
    x = np.array([L1*cos(Q1[i])+L2*cos(Q1[i]+Q2[i]) for i in range(len(Q1))])
    y = np.array([L1*sin(Q1[i])+L2*sin(Q1[i]+Q2[i]) for i in range(len(Q1))])

    plt.figure()
    plt.plot(x,y)
    ax1 = plt.gca()
    ax1.set_aspect('equal')

    plt.figure()
    plt.plot(t,x)
    ax2 = plt.gca()

    plt.figure()
    plt.plot(t,y)
    ax3 = plt.gca()

    plt.show()
def calculate_torque():
    from math import cos,sin
    import numpy as np
    m1,m2 = 1.02,1.16 # kg
    c1,c2 = 0.104,0.165 # m
    I1,I2 = 0.0167,0.0474 # kg⋅m²
    b1,b2 = 0.8,0.8 # kg⋅m²/s
    global τ1,τ2
    τ1,τ2 = [],[]
    for i in range(len(Q1)):
        I = np.matrix([[I1 + I2 + 2*m2*L1*c2*cos(Q2[i]) + m2*L1**2, I2 + m2*L1*c2*cos(Q2[i])],[I2 + m2*L1*c2*cos(Q2[i]), I2]])
        a = np.matrix([[-m2*L1*c2*(2*ddt_Q1[i]+ddt_Q2[i])*ddt_Q2[i]*sin(Q2[i]) + b1*ddt_Q1[i]],\
        [m2*L1*c2*ddt_Q1[i]**2*sin(Q2[i])+b2*ddt_Q2[i]]])
        Torques = I*np.matrix([[d2dt2_Q1[i]],[d2dt2_Q2[i]]]) + a
        τ1.append(Torques[0,0])
        τ2.append(Torques[1,0])
def plot_resulting_torques():
    import matplotlib.pyplot as plt
    import numpy as np
    from math import cos,sin

    dt = 1/(len(Q1)-1)
    t = np.arange(0,1+dt,dt)
    plt.figure()
    plt.plot(t,τ1,'r',t,τ2,'b')
    plt.show()
reaching_task()
calculate_torque()
# plot_resulting_kinematics()
plot_resulting_torques()
