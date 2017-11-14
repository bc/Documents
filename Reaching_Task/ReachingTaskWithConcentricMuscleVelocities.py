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
    Sets the global link lengths for a 2 DOF planar reaching task. Changes the values of the link lengths. New_L1 and New_L2 must be numbers. Set EOM to Uno of Zadravec.
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
	Takes in a (2,N) list/array with values for x and y (endpoint) and maps to the  A1 and A2 with the current angles (in radians) from the inverse kinematics.
	"""
	import numpy as np
	from math import acos, sin, cos, atan, atan2
	from numpy import pi
	assert np.shape(X)[0]==2, "X must be either a (2,len(X)) list/array"
	x,y = np.split(X,2,axis=0) # Returns 2 (1,N) arrays
	"""
	Math Logic:
	x² + y² > (L₁+L₂)² = L₁² + 2L₁L₂ + L₂² > L₁² + L₂²
	x² + y² - L₁² - L₂² > 0
	a₂ = cos⁻¹((x² + y² - L₁² - L₂²)/2L₁L₂) ∊ (0,π/2)

	Map functions take in (N,1) list/arrays -- i.e., only length-1 arrays can be converted to Python scalars needed for the map function. Therefore, the transpose is needed to change (1,N) to (N,1)
	"""
	a2 = lambda x,y: acos((x**2 + y**2 - L1**2 - L2**2)/(2*L1*L2))
	a1 = lambda x,y,a2: atan2(y,x) - atan2(L2*sin(a2),(L1+L2*cos(a2)))
	global A1,A2
	A2 = np.array(list(map(a2,x.T,y.T)), dtype='float64', ndmin=2) # Returns a (1,N) array
	A1 = np.array(list(map(a1,x.T,y.T,A2.T)), dtype='float64', ndmin=2) # Returns a (1,N) array
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
		Xi = [-0.177,0.15]
		Xf = [-0.177,0.15 + 0.35]

	# # Left Diagonal Reach, Path Length = 0.5295⋅(L1+L2) = 0.35
	elif ReachType.capitalize() == 'Left':
		Xi = [-0.177,0.15]
		Xf = [-0.177-0.35/(2**0.5),0.15 + 0.35/(2**0.5)]

	# # Right Diagonal Reach, Path Length = 0.5295⋅(L1+L2) = 0.35
	elif ReachType.capitalize() == 'Right':
		Xi = [-0.177,0.15]
		Xf = [-0.177+0.35/(2**0.5),0.15 + 0.35/(2**0.5)]

	return(Xi,Xf)
def find_X_values(t,Xi,Xf):
	"""
	This takes t (a numpy.ndarray (dtype = 'float64', ndmin = 2) of normalized time ∈ [0,1]) and either a (2,) or (2,1) list/array with values for both initial and final x and y (endpoint). To avoid singularities, ||X[i]|| cannot be greater than L1 + L2.

	Returns a (2,t.shape) array of dtype = 'float64', ndmin = 2.
	"""
	import numpy as np
	assert type(t) == np.ndarray, "t must be a (1,N) 'float64' numpy.ndarray"
	assert np.shape(t) == (1,t.size), "t must be a (1,N) 'float64' numpy.ndarray"
	assert t.dtype == 'float64', "t must be a (1,N) 'float64' numpy.ndarray"
	assert len(Xf)==2, "Xf must be either a (2,) or (2,1) list/array"
	assert len(Xi)==2, "Xi must be either a (2,) or (2,1) list/array"
	xi,yi = Xi
	xf,yf = Xf
	if np.shape(Xi) == (2,1): xi,yi = xi[0],yi[0]
	if np.shape(Xf) == (2,1): xf,yf = xf[0],yf[0]
	x = xi + (xf-xi)*(10*t**3 - 15*t**4 + 6*t**5) # Returns an (1,N) array
	y = yi + (yf-yi)*(10*t**3 - 15*t**4 + 6*t**5) # Returns an (1,N) array
	X = np.concatenate((x,y),axis=0) # Returns an (2,N) array

		# Test of singularities
	assert (sum(X**2)**0.5<L1+L2).all(), "Trajectory creates singularities."

	return(X)
def find_Ẋ_values(t,Xi,Xf,t_end = 1):
	"""
	This takes t (a numpy.ndarray (dtype = 'float64', ndmin = 2) of normalized time ∈ [0,1]) and either a (2,) or (2,1) list/array with values for both initial and final x and y (endpoint).

	Returns a (2,t.shape) array of dtype = 'float64', ndmin = 2.
	"""
	import numpy as np
	assert type(t) == np.ndarray, "t must be a (1,N) 'float64' numpy.ndarray"
	assert np.shape(t) == (1,t.size), "t must be a (1,N) 'float64' numpy.ndarray"
	assert t.dtype == 'float64', "t must be a (1,N) 'float64' numpy.ndarray"
	assert len(Xf)==2, "Xf must be either a (2,) or (2,1) list/array"
	assert len(Xi)==2, "Xi must be either a (2,) or (2,1) list/array"
	xi,yi = Xi
	xf,yf = Xf
	if np.shape(Xi) == (2,1): xi,yi = xi[0],yi[0]
	if np.shape(Xf) == (2,1): xf,yf = xf[0],yf[0]
	ẋ = (xf-xi)*(30*t**2 - 60*t**3 + 30*t**4)/t_end # Returns an (1,N) array
	ẏ = (yf-yi)*(30*t**2 - 60*t**3 + 30*t**4)/t_end # Returns an (1,N) array
	Ẋ = np.concatenate((ẋ,ẏ),axis=0) # Returns an (2,N) array
	return(Ẋ)
def find_Ẍ_values(t,Xi,Xf,t_end=1):
	"""
	This takes t (a numpy.ndarray (dtype = 'float64', ndmin = 2) of normalized time ∈ [0,1]) and either a (2,) or (2,1) list/array with values for both initial and final x and y (endpoint).

	Returns a (2,t.shape) array of dtype = 'float64', ndmin = 2.
	"""
	import numpy as np
	assert type(t) == np.ndarray, "t must be a (1,N) 'float64' numpy.ndarray"
	assert np.shape(t) == (1,t.size), "t must be a (1,N) 'float64' numpy.ndarray"
	assert t.dtype == 'float64', "t must be a (1,N) 'float64' numpy.ndarray"
	assert len(Xf)==2, "Xf must be either a (2,) or (2,1) list/array"
	assert len(Xi)==2, "Xi must be either a (2,) or (2,1) list/array"
	xi,yi = Xi
	xf,yf = Xf
	if np.shape(Xi) == (2,1): xi,yi = xi[0],yi[0]
	if np.shape(Xf) == (2,1): xf,yf = xf[0],yf[0]
	ẍ = (xf-xi)*(60*t - 180*t**2 + 120*t**3)/(t_end**2) # Returns an (1,N) array
	ÿ = (yf-yi)*(60*t - 180*t**2 + 120*t**3)/(t_end**2) # Returns an (1,N) array
	Ẍ = np.concatenate((ẍ,ÿ),axis=0) # Returns an (2,N) array
	return(Ẍ)
def return_X_values(t,Xi,Xf,t_end):
	X = find_X_values(t,Xi,Xf)
	Ẋ = find_Ẋ_values(t,Xi,Xf,t_end=t_end)
	Ẍ = find_Ẍ_values(t,Xi,Xf,t_end=t_end)
	return(X,Ẋ,Ẍ)
def jacobian(a1,a2):
	import numpy as np
	from math import cos, sin
	J = np.matrix([[-L1*sin(a1)-L2*sin(a1+a2),-L2*sin(a1+a2)],\
								[L1*cos(a1)+L2*cos(a1+a2),L2*cos(a1+a2)]], dtype = 'float64')
	return(J)
def inverse_jacobian(a1,a2):
	import numpy as np
	from math import cos, sin
	det_J = L1*L2*cos(a1 + a2)*sin(a1) - L1*L2*sin(a1 + a2)*cos(a1)
	J_inv = (1/det_J)*np.matrix([[- L2*cos(a1 + a2), 						- L2*sin(a1 + a2)],\
								[L1*cos(a1)+L2*cos(a1 + a2),	L2*sin(a1 + a2) + L1*sin(a1)]], \
								dtype = 'float64')
	return(J_inv)
def update_angular_velocity(Ẋ):
	import numpy as np
	global A1,A2 # both are (1,N) arrays
	"""
	Note: numpy.float64 does not support inverse matrix operators. Luckily, this is an always invertible 2x2 matrix (J is nonsingular within the ROM because a₂ ∊ (0,π/2)), so the inverse Jacobian function has been mapped instead.
	"""
	# J = list(map(jacobian,A1.T,A2.T)) # Returns a list of shape (N,2,2)
	J_inv = list(map(inverse_jacobian,A1.T,A2.T)) # Returns a list of shape (N,2,2)
	"""
	In order to map properly, X (shape (2,N)) must be be first split into N (2,1) arrays (length-1). Therefore, the third argument to map() does not need to be transposed. Similar logic follows for why J_inv is not transposed, as it is has N (2,2) matrix arrays. Changing the output to an array aids in creating arrays Ȧ1 and Ȧ2.
	"""
	Ȧ = list(map(lambda j_inv,ẋ: np.array(j_inv*ẋ), \
						J_inv, np.split(Ẋ,np.shape(Ẋ)[1],axis=1)))
	global Ȧ1,Ȧ2
	Ȧ1,Ȧ2 = np.split(np.concatenate(Ȧ, axis=1), 2, axis=0) # both are (1,N) arrays
def update_angular_acceleration(Ẋ,Ẍ):
	from math import cos,sin
	import numpy as np

	assert np.shape(Ẋ)[0]==2, "Ẋ must be a (2,len(Ẋ)) list/array"
	assert np.shape(Ẍ)[0]==2, "Ẍ must be a (2,len(Ẍ)) list/array"
	ẍ,ÿ = np.split(Ẍ,2,axis=0)
	global A1,A2,Ȧ1,Ȧ2
	"""
		dȦ1/dt = (δȦ1/δA1)*(dA1/dt)
					+ (δȦ1/δA2)*(dA2/dt)
						+ (δȦ1/δẊ[0])*(dẊ[0]/dt)
							+ (δȦ1/δẊ[1])*(dẊ[1]/dt)

		dȦ2/dt = (δȦ2/δA1)*(dA1/dt)
					+ (δȦ2/δA2)*(dA2/dt)
						+ (δȦ2/δẊ[0])*(dẊ[0]/dt)
							+ (δȦ2/δẊ[1])*(dẊ[1]/dt)
	"""
	ä1 = lambda a1,a2,ȧ1,ȧ2,ẋ,ẏ,ẍ,ÿ: \
		((-sin(a1+a2)*ẋ+cos(a1+a2)*ẏ)/(L1*sin(a2)))*ȧ1 \
		+ (((-sin(a1+a2)*ẋ+cos(a1+a2)*ẏ)*L1*sin(a2) - \
		(cos(a1+a2)*ẋ+sin(a1+a2)*ẏ)*L1*cos(a2))/((L1**2)*(sin(a2)**2)))*ȧ2 \
		+ (cos(a1+a2)/(L1*sin(a2)))*ẍ \
		+ (sin(a1+a2)/(L1*sin(a2)))*ÿ
	ä2 = lambda a1,a2,ȧ1,ȧ2,ẋ,ẏ,ẍ,ÿ: \
		(((L1*sin(a1)+L2*sin(a1+a2))*ẋ + (-L1*cos(a1)-L2*cos(a1+a2))*ẏ)/(L1*L2*sin(a2)))*ȧ1 \
		+ (((L2*sin(a1+a2)*ẋ + (-L2*cos(a1+a2))*ẏ)*(L1*L2*sin(a2)) \
		- ((-L1*cos(a1)-L2*cos(a1+a2))*ẋ + (-L1*sin(a1)-L2*sin(a1+a2))*ẏ)*(L1*L2*cos(a2)))\
		/((L1*L2*sin(a2))**2))*ȧ2 \
		+ ((-L1*cos(a1)-L2*cos(a1+a2))/(L1*L2*sin(a2)))*ẍ + \
		((-L1*sin(a1)-L2*sin(a1+a2))/(L1*L2*sin(a2)))*ÿ

	global Ä1,Ä2
	Ä1 = np.array(list(map(ä1,A1.T,A2.T,Ȧ1.T,Ȧ2.T,Ẋ[0].T,Ẋ[1].T,Ẍ[0].T,Ẍ[1].T))\
					,dtype = 'float64',ndmin=2).T # returns a (1,N) array
	Ä2 = np.array(list(map(ä2,A1.T,A2.T,Ȧ1.T,Ȧ2.T,Ẋ[0].T,Ẋ[1].T,Ẍ[0].T,Ẍ[1].T))\
					,dtype = 'float64',ndmin=2).T # returns a (1,N) array
def update_angle_lists(X,Ẋ,Ẍ):
	"""
	Takes in three (2,N) endpoint arrays and returns global lists for angles 1 and 2 of shape (1,N).
	"""
	import numpy as np
	inverse_kinematics(X)
	update_angular_velocity(Ẋ)
	update_angular_acceleration(Ẋ,Ẍ)
def test_inverse_kinematics(X,Ẋ,Ẍ):
	import numpy as np
	from math import cos,sin
	import matplotlib.pyplot as plt
	recovered_x = np.array(list(map(lambda a1,a2: L1*cos(a1) + L2*cos(a1+a2),A1.T,A2.T)),\
	 						dtype='float64',ndmin=2)
	recovered_y = np.array(list(map(lambda a1,a2: L1*sin(a1) + L2*sin(a1+a2),A1.T,A2.T)),\
	 						dtype='float64',ndmin=2)

	recovered_dX = np.concatenate([np.array(jacobian(A1[0,i],A2[0,i])*np.matrix([[Ȧ1[0,i]],[Ȧ2[0,i]]])) for i in range(np.shape(A1)[1])],axis=1)
	recovered_dx,recovered_dy = np.split(recovered_dX,2,axis=0)

	recovered_ddx = np.array(list(map(lambda a1,a2,ȧ1,ȧ2,ä1,ä2:\
	 	-L1*cos(a1)*(ȧ1**2)-L1*sin(a1)*ä1-L2*cos(a1+a2)*((ȧ1+ȧ2)**2)-L2*sin(a1+a2)*(ä1+ä2),\
		A1.T,A2.T,Ȧ1.T,Ȧ2.T,Ä1.T,Ä2.T)),dtype = 'float64',ndmin=2)
	recovered_ddy = np.array(list(map(lambda a1,a2,ȧ1,ȧ2,ä1,ä2:\
	 	-L1*sin(a1)*(ȧ1**2)+L1*cos(a1)*ä1-L2*sin(a1+a2)*((ȧ1+ȧ2)**2)+L2*cos(a1+a2)*(ä1+ä2),\
		A1.T,A2.T,Ȧ1.T,Ȧ2.T,Ä1.T,Ä2.T)),dtype = 'float64',ndmin=2)

	plt.figure()
	ax1 = plt.gca()
	ax1.set_title('Recovered Position')
	plt.plot(X[0],'r', recovered_x[0],'b')
	plt.plot(X[1],'r', recovered_y[0],'b')

	plt.figure()
	ax2 = plt.gca()
	ax2.set_title('Recovered Velocity')
	plt.plot(Ẋ[0],'r', recovered_dx[0],'b')
	plt.plot(Ẋ[1],'r', recovered_dy[0],'b')

	plt.figure()
	ax3 = plt.gca()
	ax3.set_title('Recovered Acceleration')
	plt.plot(Ẍ[0],'r', recovered_ddx,'b')
	plt.plot(Ẍ[1],'r', recovered_ddy,'b')
	plt.show()
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
	C_matrix = lambda a1,a2,ȧ1,ȧ2: \
	        np.matrix([ [-β*ȧ2*sin(a2),     -β*(ȧ1 + ȧ2)*sin(a2)],
	                    [β*ȧ1*sin(a2),      0]]) # kg⋅m² (N⋅m⋅s²)
	M_matrix = lambda a1,a2: \
	        np.matrix([ [α + 2*β*cos(a2),   δ + β*cos(a2)],
	                    [δ + β*cos(a2),     δ]],\
						dtype = 'float64') # kg⋅m² (N⋅m⋅s²)
	B_matrix = np.matrix([ [b11, b12],\
	                	[b21, b22]]) # kg⋅m²/s (N⋅m⋅s)
	global T1,T2
	T1,T2 = [],[]
	Ȧ = np.swapaxes(np.array(np.concatenate((Ȧ1,Ȧ2),axis=0),ndmin=3),0,2)
	Ä = np.swapaxes(np.array(np.concatenate((Ä1,Ä2),axis=0),ndmin=3),0,2)
	M = np.array(list(map(M_matrix,A1.T,A2.T)))
	C = np.array(list(map(C_matrix,A1.T,A2.T,Ȧ1.T,Ȧ2.T)))
	MÄ = np.array(list(map(lambda m,ä: np.matrix(m)*ä,M,Ä)))
	CȦ = np.array(list(map(lambda c,ȧ: np.matrix(np.array(c,ndmin=2,dtype='float64'))*ȧ,C,Ȧ)))
	BȦ = np.array(list(map(lambda ȧ: B_matrix*ȧ,Ȧ)))
	T = MÄ + CȦ + BȦ
	# returns a (N,1,2) 3D array. Therefore we must transpose it first and select first element
	T1,T2 = np.split(T.T[0],2,axis=0) # returns 2 (1,N) arrays
def calculate_potential_variability(X,δT1,δT2,dt=0.001,EOM = 'Uno',scheme = "Total"):
	from math import cos,sin
	import numpy as np
	import ipdb

	assert EOM.capitalize() in ["Uno","Zadravec"], "EOM can be either 'Uno' or 'Zadravec'"
	assert scheme.capitalize() in ["Total","Individual"], "scheme must be either 'Total' or 'Individual'"

	global A1,A2,Ȧ1,Ȧ2,T1,T2,MomentArmMatrix

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

	M_inv_matrix = lambda a1,a2: (1/(δ*(α-δ)-((β*cos(a2))**2)))* \
	        np.matrix([ [δ,   -(δ + β*cos(a2))],
	                    [-(δ + β*cos(a2)),     α + 2*β*cos(a2)]],\
						dtype = 'float64') # kg⋅m² (N⋅m⋅s²)
	C_matrix = lambda a1,a2,ȧ1,ȧ2: \
			np.matrix([ [-β*ȧ2*sin(a2),     -β*(ȧ1 + ȧ2)*sin(a2)],
						[β*ȧ1*sin(a2),      0]]) # kg⋅m² (N⋅m⋅s²)
	B_matrix = np.matrix([ [b11, b12],\
						[b21, b22]]) # kg⋅m²/s (N⋅m⋅s)
	A = np.swapaxes(np.array(np.concatenate((A1,A2),axis=0),ndmin=3),0,2)
	Ȧ = np.swapaxes(np.array(np.concatenate((Ȧ1,Ȧ2),axis=0),ndmin=3),0,2)
	Ä = np.swapaxes(np.array(np.concatenate((Ä1,Ä2),axis=0),ndmin=3),0,2)
	M_inv = np.array(list(map(M_inv_matrix,A1.T,A2.T)))
	C = np.array(list(map(C_matrix,A1.T,A2.T,Ȧ1.T,Ȧ2.T)))
	CȦ = np.array(list(map(lambda c,ȧ: np.matrix(np.array(c,ndmin=2,dtype='float64'))*ȧ,C,Ȧ)))
	BȦ = np.array(list(map(lambda ȧ: B_matrix*ȧ,Ȧ)))
	T = np.swapaxes(np.array(np.concatenate((T1,T2),axis=0),ndmin=3),0,2)

	# Test to make sure that the forward mapping to torque and the inverse mapping to angular acceleration are consistent.

	assert (abs(np.array([np.matrix(M_inv[i])*(T[i]-CȦ[i]-BȦ[i]) for i in range(1001)]).T -Ä.T)<1e-12).all(), "There is an issue with the mapping from Torques to Angular Accelerations."

	# The sign of the MA for each muscle at each joint will sufficiently determine the sign of of the torque as these muscles are already eccentrically contracting (i.e., -r₁ⱼ⋅dϑ₁/dt > 0, therefore if dϑ₁/dt > 0, then r₁ⱼ < 0 and the torque will be negative and vice versa).

	δT = (np.concatenate((δT1,δT2),axis=1)*np.sign(MomentArmMatrix.T)).T
	# returns an (N,2,m) array.

	from itertools import chain, combinations
	def all_subsets(ss):
		return chain(*map(lambda x: combinations(ss, x), range(0, len(ss)+1)))
	LengtheningMuscleFilter = list(\
								filter(lambda i: not ((δT1[i]<=0).all() and (δT2[i]<=0).all()),\
									range(16))\
										)
	AllSubsets = []
	for subset in all_subsets(LengtheningMuscleFilter):
		AllSubsets.append(list(subset))

	if scheme.capitalize() == 'Individual':
		next_Ä = np.concatenate(\
					[np.array([np.matrix(M_inv[i])*(T[i]+δT[i,:,j][:,np.newaxis]-CȦ[i]-BȦ[i]) \
						for i in range(δT.shape[0])]).T for j in range(δT.shape[2])]) # returns (m,2,1001) array
	elif scheme.capitalize() == 'Total':
		next_Ä = np.array([np.matrix(M_inv[i])*(T[i]+δT.sum(axis=2)[:,np.newaxis].swapaxes(1,2)[i]\
							-CȦ[i]-BȦ[i]) \
								for i in range(δT.shape[0])]).T # returns (1,2,N) array


	"""
	A,Ȧ,Ä have shape (N,2,1), but next_Ä has shape (x,2,N), where x is determined by the torque combinations. We have chosen to allow for variable torque contributes in case we wish to observe the effects of individual muscle, muscle groups, etc. Therefore we use next_Ä[j].T[:,np.newaxis].swapaxes(1,2) (where j is in range(x)) to take (N,2) arrays, add a new axis in the middle, and then swap 2nd and 3rd axis to regain (N,2,1) arrays for each j in range(x). Taking the transpose at the end returns an easy to split (1,2,N) array.
	"""
	resulting_A = np.concatenate([(A+(Ȧ+((Ä + next_Ä[j].T[:,np.newaxis].swapaxes(1,2))/2)*dt/2)*dt).T \
									for j in range(next_Ä.shape[0])])
		# returns (x,1,N) arrays for the resulting A1 and A2.
	resulting_A1,resulting_A2 = np.split(resulting_A,resulting_A.shape[1],axis=1)
		# returns (x,N) arrays for the resulting deviation at step n+1.
	resulting_x = np.concatenate([np.array(list(map(lambda a1,a2: L1*cos(a1) + L2*cos(a1+a2),\
										resulting_A1[j].T,resulting_A2[j].T)),\
	 										dtype='float64',ndmin=2) \
												for j in range(resulting_A.shape[0])])
	resulting_y = np.concatenate([np.array(list(map(lambda a1,a2: L1*sin(a1) + L2*sin(a1+a2),\
										resulting_A1[j].T,resulting_A2[j].T)),\
	 										dtype='float64',ndmin=2) \
												for j in range(resulting_A.shape[0])])
	"""
	If scheme is "Total", then the output will be a (1,N) array of total resulting potential variability vs. time. If scheme is "Individual", then the output will be a (m,N) array of the individual muscle potential variability contribution versus time.
	"""
	potential_variability = ((X[0][1:]-resulting_x[:,:-1])**2 +\
	 								(X[1][1:]-resulting_y[:,:-1])**2)**0.5 # returns (x,N-1) array because we do not have any variability at the first timestep.
	potential_variability = \
			np.concatenate((np.array([0]*potential_variability.shape[0],ndmin=2).T,\
								potential_variability),\
									axis=1) # returns an (x,N) array that adds zero values infront of index 0.
	# ipdb.set_trace()
	return(potential_variability)
def reaching_task(dt = 0.001, Xi = [0,0.25], Xf = [0.25,0.50],t_end=1,EOM='Uno'):
	import numpy as np
	set_link_lengths()
	create_angle_lists()
	t = np.array(np.arange(0,1+dt,dt),dtype='float64',ndmin=2)
	X,Ẋ,Ẍ = return_X_values(t,Xi,Xf,t_end)
	update_angle_lists(X,Ẋ,Ẍ)
	calculate_torques(EOM=EOM)
	return(X)
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
		MomentArm = (np.matrix(Coefficients,dtype='float64')\
						*np.matrix([1,q,q**2,q**3,q**4,q**5]).T)[0,0]
	elif src.capitalize() == 'Est' :
		MomentArm = np.array(Coefficients,dtype='float64')
	else: #src.capitalize() == 'Ramsay'
		q = q2
		assert type(Coefficients) == list, "Coefficients must be a list."
		assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
		assert eq in [1,2,3], "eq must be either 1, 2, or 3 when using Ramsay (2009)."
		if eq == 1:
			assert len(Coefficients) == 5, "For Eq. 1, Coefficients must be 5 elements long."
			MomentArm = (sp.Matrix(Coefficients,dtype='float64').T\
							*sp.Matrix([1,q,q**2,q**3,q**4]))[0,0]
		elif eq == 2:
			assert len(Coefficients)==16, "For Eq. 2, Coefficients must be 16 elements long."
			MomentArm = (sp.Matrix(Coefficients,dtype='float64').T*\
							sp.Matrix([1, q, q_PS, q*q_PS, q**2, \
										q_PS**2, (q**2)*q_PS, q*(q_PS**2), \
										(q**2)*(q_PS**2), q**3, q_PS**3, \
										(q**3)*q_PS, q*(q_PS**3), \
										(q**3)*(q_PS**2), (q**2)*(q_PS**3), \
										(q**3)*(q_PS**3)]))[0, 0]
		else: # eq == 3
			assert len(Coefficients)==18, "For Eq. 3, Coefficients must be 18 elements long."
			MomentArm = (sp.Matrix(Coefficients,dtype='float64').T*\
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
			np.array([1,(180/np.pi),(180/np.pi)**2,(180/np.pi)**3,(180/np.pi)**4,(180/np.pi)**5],dtype = 'float64'))
	new_Coefficients =\
	 	np.multiply(rad_conversion,np.array([1,1e-1,1e-3,1e-5,1e-7,1e-9],dtype='float64'))
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

	PC (Clavicle attachment of Pectoralis) SFE MA is listed as 50.80 mm in Pigeon. APPROX OPTIMAL MUSCLE LENGTH! NEED TO FIND ACTUAL NUMBER. We used the Banks numbers for mass, afferent number, corrected number and relative abundance as the stretch will likely affect the whole muscle.

	CB SFE MA was estimated in Holzbaur (2005) as 20 mm while Bassett (1990) estimates from 7 cadavers the MA to be 36 mm.
	"""
	import sympy as sp
	from sympy.utilities import lambdify
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
		'Mass' :  295.6, \
		'Actual No' : 450, \
		'Corrected No' : 389.7, \
		'Relative Abundance' : 1.2,\
		'Optimal Muscle Length' : 150, \
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
	AllCoefficients = {'PC': PC_Coefficients,\
						'DELTa' : DELTa_Coefficients, 'CB' : CB_Coefficients, \
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

	global RMatrix_Transpose, dRMatrix_Transpose
	RMatrix_Transpose = sp.Matrix([[MA_function(AllCoefficients[muscle][dof]) for dof in ['Shoulder','Elbow']] for muscle in MuscleList])
	dRMatrix_Transpose = sp.Matrix(np.concatenate((sp.diff(RMatrix_Transpose[:,0],q1),\
											sp.diff(RMatrix_Transpose[:,1],q2)),axis=1))
	RMatrix_Transpose = lambdify([q1,q2,q_PS],RMatrix_Transpose)
	dRMatrix_Transpose = lambdify([q1,q2,q_PS],dRMatrix_Transpose)
	# returns an (n,m) matrix when n is the number of muscles and m is the number of DOFS. We chose to return R.T because this is commonly utilized in muscle velocity calculations.
def plot_MA_values_for_muscle_(j):
	import matplotlib.pyplot as plt
	import numpy as np
	global MomentArmMatrix
	assert j in range(MomentArmMatrix.shape[2]), "j must be the index of one of the muscles used (j ∊ {0,1,...,m-1})"
	plt.figure()
	ax = plt.gca()
	plt.plot(MomentArmMatrix[:,:,j])
	ax.set_xticks([0,MomentArmMatrix.shape[0]])
	ax.set_xticklabels(['Start','Finish'])
	ax.set_ylabel('Moment Arm (in mm)')
	plt.show()
def return_MA_matrix():
	"""
	Notes:

	The angle of pronation is fixed at pi/2 for this reaching paradigm.
	The angle of radial/ulnar deviation is set to zero for the fixed wrist apparatus of the model paradigm.
	These functions have been verified to match the previous posture dependent MAs from Ramsey (2010) - SEE ERRATUM
	"""
	import numpy as np
	from numpy import pi
	import sympy as sp

	global A1,A2,RMatrix_Transpose,dRMatrix_Transpose
	global MomentArmMatrix,dMomentArmMatrix
	MomentArmMatrix = np.array(list(map(lambda A1,A2: \
						np.float64(RMatrix_Transpose(A1,A2,pi/2).T),\
						A1.T,A2.T)))
	dMomentArmMatrix = np.array(list(map(lambda A1,A2:\
						np.float64(dRMatrix_Transpose(A1,A2,pi/2).T),\
						A1.T,A2.T)))
	# returns two matrices of size (N,2,m)
def plot_moment_arm_functions(OrderNumber, ROM=[[-10,130],[0,160]]):
	import matplotlib.pyplot as plt
	import numpy as np
	from numpy import pi

	global A1,RMatrix_Transpose

	A1_rom = np.linspace(ROM[0][0]*pi/180,ROM[0][1]*pi/180,A1.shape[1])
	A2_rom = np.linspace(ROM[1][0]*pi/180,ROM[1][1]*pi/180,A1.shape[1])

	MomentArmMatrix = np.array(list(map(lambda A1,A2: \
						np.float64(RMatrix_Transpose(A1,A2,pi/2).T),\
						A1_rom,A2_rom)))

	fig1 = plt.figure()
	[plt.plot(A1_rom,MomentArmMatrix[:,0,i].T) for i in OrderNumber]
	ax1 = plt.gca()
	ax1.set_ylim(-30,70)
	ax1.set_xlim(ROM[0][0]*pi/180,ROM[0][1]*pi/180*1.3)
	ax1.set_xlabel('Shoulder Angle')
	ax1.set_ylabel('Moment Arm (mm)')
	ax1.set_title('Moment Arm Functions: Shoulder')
	[j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax1.lines)]
	ax1.legend(OrderedMuscleList)

	fig2 = plt.figure()
	[plt.plot(A2_rom,MomentArmMatrix[:,1,i].T) for i in OrderNumber]
	ax2 = plt.gca()
	ax2.set_ylim(-30,70)
	ax2.set_xlim(ROM[1][0]*pi/180,ROM[1][1]*pi/180*1.3)
	ax2.set_xlabel('Elbow Angle')
	ax2.set_ylabel('Moment Arm (mm)')
	ax2.set_title('Moment Arm Functions: Elbow')
	[j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax2.lines)]
	ax2.legend(OrderedMuscleList)

	plt.show()
def calculate_muscle_velocities():
	"""
	v = (-dR/dϑ)⋅(dϑ/dt)⋅ϑ + (-R)⋅(dϑ/dt)
	(dsⱼ/dϑ₁)⋅(dϑ₁/dt) = √((dϑ₁/dt)²⋅(dr₁ⱼ/dϑ₁)² + (dϑ₁/dt)²⋅(r₁ⱼ)²)
	(dsⱼ/dt) = sum([√(el) for el in diag(dϑ/dt)(J²_{Rⱼ})(dϑ/dt) + diag(Rⱼ)(diag(dϑ/dt)²)(dϑ₁/dt)])

	Need to hard code whether the rotation causes lengthening or shortening. This is determined by the sign of -r₁ⱼ⋅dϑ₁/dt for each d.o.f.
	"""
	import numpy as np
	global Ȧ1,Ȧ2,AllCoefficients,MomentArmMatrix,dMomentArmMatrix
	Ȧ = np.swapaxes(np.array(np.concatenate((Ȧ1,Ȧ2),axis=0),ndmin=3),0,2)
	MuscleVelocity = (abs(Ȧ.T)*(((dMomentArmMatrix.T**2)+(MomentArmMatrix.T**2))**0.5)\
						*np.sign(-Ȧ.T*MomentArmMatrix.T))\
							.sum(axis=1)[:,np.newaxis]
	OptimalMuscleLength = np.array([AllCoefficients[key]['Optimal Muscle Length'] for key in AllCoefficients.keys()],dtype='float64')[:,np.newaxis,np.newaxis]
	NormalizedMuscleVelocity = MuscleVelocity/OptimalMuscleLength
	# returns a (m,1,N) array where m is the number of muscles and N is the number of timesteps.
	return(NormalizedMuscleVelocity)
def calculate_weighted_unscaled_potential_torque_variations_REACHING_TASK():
	"""
	v = (-dR/dϑ)⋅(dϑ/dt)⋅ϑ + (-R)⋅(dϑ/dt)
	(dsⱼ/dϑ₁)⋅(dϑ₁/dt) = √((dϑ₁/dt)²⋅(dr₁ⱼ/dϑ₁)² + (dϑ₁/dt)²⋅(r₁ⱼ)²)
	(dsⱼ/dt) = sum([√(el) for el in diag(dϑ/dt)(J²_{Rⱼ})(dϑ/dt) + diag(Rⱼ)(diag(dϑ/dt)²)(dϑ₁/dt)])

	Need to hard code whether the rotation causes lengthening or shortening. This is determined by the sign of -r₁ⱼ⋅dϑ₁/dt for each d.o.f.
	"""
	import numpy as np
	import ipdb
	global Ȧ1,Ȧ2,AllCoefficients,MomentArmMatrix,dMomentArmMatrix

	MuscleList = list(AllCoefficients.keys())
	Ȧ = np.swapaxes(np.array(np.concatenate((Ȧ1,Ȧ2),axis=0),ndmin=3),0,2)

	"""
	Note: Ȧ is a 3D array of size (N,2,1). Thus, Ȧ.T has shape (1,2,N). Numpy.ndarray multiplication of arrays with IDENTICAL shapes will result in an element by element multiplication similar to np.multiply. Thus we multiply -Ȧ.T[0] by the (2,N) array that is the result of the norm of dMomentArmMatrix.T[j] and MomentArmMatrix.T[j] (each of shape (2,N)). This is equivalent to -ȧ₁*√(ṙ₁ⱼ² + r₁ⱼ²) for all t.
	"""

	MuscleVelocity = \
					(abs(Ȧ.T)*(((dMomentArmMatrix.T**2)+(MomentArmMatrix.T**2))**0.5)\
						*np.sign(-Ȧ.T*MomentArmMatrix.T))\
							.sum(axis=1)[:,np.newaxis]

	"""
	Normalizing these muscle velocities can similarly be done by dividing the (m,1,N) MuscleVelocity array by the OptimalMuscleLength array of shape (m,1,1).
	"""
	CurveCoefficient = 2
	OptimalMuscleLength = np.array([AllCoefficients[key]['Optimal Muscle Length'] for key in AllCoefficients.keys()],dtype='float64')[:,np.newaxis,np.newaxis]
	ConcentricForceVelocityCurve = 1 - (abs(MuscleVelocity)/(5*OptimalMuscleLength))**CurveCoefficient
	# NormalizedMuscleVelocity = MuscleVelocity/OptimalMuscleLength

	"""
	Next we need to multiply each muscle velocity, AT EACH TIMESTEP, by the appropriate scaling factor, which may be a function of joint angle (i.e., movement step number).
	"""

	CorrectedAfferentNumber = np.array([AllCoefficients[key]['Corrected No'] for key in AllCoefficients.keys()],dtype='float64')[:,np.newaxis] # returns an (m,1) array

	"""
	Note: (NormalizedMuscleVelocity*MomentArmMatrix.T/1000).T returns an (N,2,m) matrix. Therefore, we must construct a (2,m) matrix of Corrected Afferented Numbers that have identical rows. This is done by concatenating CorrectedAfferentNumber with itself and then transposing the result. This total product will yield a (m,2,N) array of both shoulder and elbow weighted muscle velocities. This will then be split by np.split() to give (m,1,N) arrays that are then squeezed into (m,N) arrays for the individual joint weighted muscle velocities.

	MomentArmMatrix must be rectified so that the contraction direction is maintained.
	"""

	# MAWeightedMuscleVelocity = \
	# 	((NormalizedMuscleVelocity*abs(MomentArmMatrix).T/1000).T\
	# 		*np.concatenate((CorrectedAfferentNumber,CorrectedAfferentNumber),axis=1).T).T
	SignConvention = np.sign(-Ȧ.T*MomentArmMatrix.T)
	SignConvention[:,:,0],SignConvention[:,:,-1] = SignConvention[:,:,1],SignConvention[:,:,-2]
	MAWeightedMuscleVelocity = \
		(ConcentricForceVelocityCurve*abs(MomentArmMatrix).T/1000)*SignConvention
	MAWeightedMuscleVelocity_Shoulder,MAWeightedMuscleVelocity_Elbow = \
	 	np.split(MAWeightedMuscleVelocity,2,axis=1)
	MAWeightedMuscleVelocity_Shoulder,MAWeightedMuscleVelocity_Elbow = \
		MAWeightedMuscleVelocity_Shoulder.squeeze(),MAWeightedMuscleVelocity_Elbow.squeeze()
	# returns two (m,N) arrays.
	return(MAWeightedMuscleVelocity_Shoulder,MAWeightedMuscleVelocity_Elbow)
def calculate_weighted_muscle_velocities_REACHING_TASK():
	"""
	v = (-dR/dϑ)⋅(dϑ/dt)⋅ϑ + (-R)⋅(dϑ/dt)
	(dsⱼ/dϑ₁)⋅(dϑ₁/dt) = √((dϑ₁/dt)²⋅(dr₁ⱼ/dϑ₁)² + (dϑ₁/dt)²⋅(r₁ⱼ)²)
	(dsⱼ/dt) = sum([√(el) for el in diag(dϑ/dt)(J²_{Rⱼ})(dϑ/dt) + diag(Rⱼ)(diag(dϑ/dt)²)(dϑ₁/dt)])

	Need to hard code whether the rotation causes lengthening or shortening. This is determined by the sign of -r₁ⱼ⋅dϑ₁/dt for each d.o.f. CORRECTION: by removing the angular velocity from the square root, the sign of the rotation is preserved -- regaining the sign convention.
	"""

	import numpy as np
	import ipdb
	global Ȧ1,Ȧ2,AllCoefficients,MomentArmMatrix,dMomentArmMatrix

	MuscleList = list(AllCoefficients.keys())
	Ȧ = np.swapaxes(np.array(np.concatenate((Ȧ1,Ȧ2),axis=0),ndmin=3),0,2)

	"""
	Note: Ȧ is a 3D array of size (N,2,1). Thus, Ȧ.T has shape (1,2,N). Numpy.ndarray multiplication of arrays with IDENTICAL shapes will result in an element by element multiplication similar to np.multiply. Thus we multiply -Ȧ.T[0] by the (2,N) array that is the result of the norm of dMomentArmMatrix.T[j] and MomentArmMatrix.T[j] (each of shape (2,N)). This is equivalent to -ȧ₁*√(ṙ₁ⱼ² + r₁ⱼ²) for all t.

	Normalizing these muscle velocities can similarly be done by dividing the (m,1,N) MuscleVelocity array by the OptimalMuscleLength array of shape (m,1,1).
	"""

	MuscleVelocity = \
					(abs(Ȧ.T)*(((dMomentArmMatrix.T**2)+(MomentArmMatrix.T**2))**0.5)\
						*np.sign(-Ȧ.T*MomentArmMatrix.T))\
							.sum(axis=1)[:,np.newaxis]

	OptimalMuscleLength = np.array([AllCoefficients[key]['Optimal Muscle Length'] for key in AllCoefficients.keys()],dtype='float64')[:,np.newaxis,np.newaxis]
	NormalizedMuscleVelocity = MuscleVelocity/OptimalMuscleLength

	"""
	Next we need to multiply each muscle velocity, AT EACH TIMESTEP, by the appropriate scaling factor, which may be a function of joint angle (i.e., movement step number).
	"""

	CorrectedAfferentNumber = np.array([AllCoefficients[key]['Corrected No'] for key in AllCoefficients.keys()],dtype='float64')[:,np.newaxis]
	WeightedMuscleVelocity = (NormalizedMuscleVelocity*abs(MomentArmMatrix).T/1000).sum(axis=1)  \
								/((abs(MomentArmMatrix)>0).T.sum(axis=1))  \
									*CorrectedAfferentNumber

	"""
	If you wanted to change the weighting scheme, use the above notation whereby the numpy.ndarrays have the same primary shape (i.e., length). Removing the .sum(axis=1) for the first line of the WeightedMuscleVelocity equation will return the individual weighted muscle velocity components per joint.
	"""
	return(WeightedMuscleVelocity)
def eccentric_velocities(NormalizedMuscleVelocity):
	"""
	Returns a (m,1,N) array of positive muscle velocities.
	"""
	import numpy as np
	def positive_entries(NormalizedMuscleVelocity):
		return(np.multiply(np.array([[1]*np.shape(NormalizedMuscleVelocity)[1]]),\
										NormalizedMuscleVelocity>0))
	PositiveMuscleVelocities = np.array(list(map(lambda Vm: np.multiply(positive_entries(Vm),Vm),\
									np.split(NormalizedMuscleVelocity,\
									len(NormalizedMuscleVelocity),axis=0))))
	return(PositiveMuscleVelocities)
def concentric_velocities(NormalizedMuscleVelocity):
	"""
	Returns a (m,1,N) array of negative muscle velocities.
	"""
	def negative_entries(NormalizedMuscleVelocity):
		return(np.multiply(np.array([[1]*np.shape(NormalizedMuscleVelocity)[1]]),\
										NormalizedMuscleVelocity<0))
	NegativeMuscleVelocities = np.array(list(map(lambda Vm: np.multiply(negative_entries(Vm),Vm),\
									np.split(NormalizedMuscleVelocity,\
									len(NormalizedMuscleVelocity),axis=0))))
	return(NegativeMuscleVelocities)
def cost_function(X,costtype="avg"):
	"""
	X must be an numpy.ndarray or size (m,). Returns a scalar.
	"""
	assert costtype in ['avg','sos','l2norm','l1norm'], "costtype must be either 'avg','sos','l2norm', or 'l1norm'"
	if costtype == 'avg' :
		cost = abs(X.mean())
	elif costtype == 'sos' :
		cost = (X**2).sum()
	elif costtype == 'l1norm' :
		cost = abs(X.sum())
	elif costtype == 'l2norm' :
		cost = (X**2).sum()**0.5
	return(cost)
def eccentric_cost(NormalizedMuscleVelocity,t_end = 1, dt = 0.001,costtype ='l2norm'):
	import numpy as np
	PositiveMuscleVelocities = eccentric_velocities(NormalizedMuscleVelocity)
	TotalPositiveExcursion = np.trapz(PositiveMuscleVelocities,dx=t_end*dt)
	EccentricCost = cost_function(TotalPositiveExcursion,costtype=costtype)
	return(EccentricCost)
def concentric_cost(NormalizedMuscleVelocity,t_end = 1, dt = 0.001,costtype = 'l2norm'):
	import numpy as np
	NegativeMuscleVelocities = concentric_velocities(NormalizedMuscleVelocity)
	TotalNegativeExcursion = np.trapz(NegativeMuscleVelocities,dx=t_end*dt)
	ConcentricCost = cost_function(TotalNegativeExcursion,costtype=costtype)
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
import ipdb
import sys

	# Save Output Figures?

SaveOutputFigures = False

	# DescriptiveTitle should be something to identify the trial either by reaching location, (i.e., Left, Right, Center, or Sideways) or by what has changed in the most current iteration (e.g., CB_Ramsay, DELTa_Est, etc.). Spaces will be replaced by '_' symbols for the filename but kept for figure titles.

ReachType = ['Sideways','Center','Left','Right'][1]
DescriptiveTitle = ReachType + ' Reach'

	# Define the Equations of Motion to be used. Can be either 'Zadravec' or 'Uno'

EOM = 'Uno'

	# Define Moment Arm Matrix Functions for the model and establish flexor/extensor list

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

	# Open PDF files to write plots to

if SaveOutputFigures == True:
	pdf_forward = PdfPages(DescriptiveTitle.replace(' ','_') + '_Forward.pdf')
	pdf_reverse = PdfPages(DescriptiveTitle.replace(' ','_') + '_Reverse.pdf')
	pdf_compare = PdfPages(DescriptiveTitle.replace(' ','_') + '_Compare.pdf')

	# Set reaching duration and sampling period

t_end = 1
dt = 0.001
t = np.array(np.arange(0,1+dt,dt,dtype = 'float64')*t_end, ndmin = 2) # returns a (1,(1+dt)/dt) array of ndmin = 2.

	# Define Initial and Final positions

Xi,Xf = set_initial_and_final_positions(ReachType)

###################################################################################################

	# Forward Direction

X_Forward = reaching_task(Xi=Xi, Xf=Xf,dt=dt, t_end=t_end,EOM = EOM)
return_MA_matrix()

	# Calculate resulting muscle velocities, lengths, etc.

NormalizedMuscleVelocity_Forward = calculate_muscle_velocities()
WeightedPotentialTorqueVariation_shoulder_Forward,\
WeightedPotentialTorqueVariation_elbow_Forward = \
calculate_weighted_unscaled_potential_torque_variations_REACHING_TASK()
WeightedNormalizedMuscleVelocity_Forward = calculate_weighted_muscle_velocities_REACHING_TASK()
# calculate_muscle_lengths()
# MuscleLengths_Forward = MuscleLengths

	# Calculate only the lengthening Torque Variation Contributions for both the shoulder and the elbow.

ScalingFactor = 10
# ScalingFactor is in units of N⋅s in order for its product to result in N⋅m⋅s⋅(Af#⋅lô/s) (Note: Af# and lô are unitless measures of corrected afferented number and normalized muscle length, respectively.

EccentricTorqueVariations_Shoulder_Forward = \
	ScalingFactor*concentric_velocities(WeightedPotentialTorqueVariation_shoulder_Forward)
EccentricTorqueVariations_Elbow_Forward =\
 	ScalingFactor*concentric_velocities(WeightedPotentialTorqueVariation_elbow_Forward)

PotentialVariability_Forward = calculate_potential_variability(X_Forward,\
									EccentricTorqueVariations_Shoulder_Forward,\
										EccentricTorqueVariations_Elbow_Forward,\
											dt=dt,EOM = EOM,scheme = "Total")

	# Plot Afferent-Weighted Muscle Velocity (Forward)

fig1e = plt.figure()
[plt.plot(t.T,WeightedNormalizedMuscleVelocity_Forward[i].T) for i in OrderNumber]
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

	# Plot Normalized Muscle Velocity (Forward)

fig1a = plt.figure()
[plt.plot(t.T,NormalizedMuscleVelocity_Forward[i].T) for i in OrderNumber]
ax1a = plt.gca()
ax1a.set_xlim(0,t_end*(1.3))
ax1a.set_ylim(-1.5,1.5)
ax1a.set_title(DescriptiveTitle+'\nNormalized Muscle Velocity')
if t_end!=1:
	ax1a.set_xlabel('Time (s)')
else:
	ax1a.set_xlabel('Normalized Time')
ax1a.set_ylabel('Normalized Muscle Velocity\nConcentric $\longleftrightarrow$ Eccentric')
[j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax1a.lines)]
ax1a.legend(OrderedMuscleList)

	# Plot Joint Angles (Forward)

fig1b = plt.figure()
plt.plot(t.T,A1.T,'g')
plt.plot(t.T,A2.T,'g:')
ax1b = plt.gca()
ax1b.set_xlim(0,t_end*(1.3))
ax1b.set_title(DescriptiveTitle+'\nAngle vs. Time')
if t_end!=1:
	ax1b.set_xlabel('Time (s)')
else:
	ax1b.set_xlabel('Normalized Time')
ax1b.set_ylabel('Joint Angles (in radians)')
ax1b.legend(["Shoulder","Elbow"])

	# Plot Configuration Space (Forward)

fig1c = plt.figure()
plt.plot(A1.T,A2.T,'g')
ax1c = plt.gca()
ax1c.set_title(DescriptiveTitle+'\nConfiguration Space')
ax1c.set_ylabel("Shouler Angle (in radians)")
ax1c.set_xlabel("Elbow Angle (in radians)")

	# Calculate Costs from Afferent-Weighted Muscle Velocities

EccentricCost_Forward = eccentric_cost(WeightedNormalizedMuscleVelocity_Forward,t_end=t_end,dt = dt,costtype='l1norm')
ConcentricCost_Forward = concentric_cost(WeightedNormalizedMuscleVelocity_Forward,t_end=t_end,dt = dt,costtype='l1norm')

# # Plot Muscle Lengths (Forward)
#
# fig1d = plt.figure()
# [plt.plot(t,MuscleLengths_Forward[i]) for i in OrderNumber]
# ax1d = plt.gca()
# ax1d.set_xlim(0,t_end*(1.3))
# ax1d.set_title(DescriptiveTitle+'\nMuscle Lengths')
# if t_end!=1:
# 	ax1d.set_xlabel('Time (s)')
# else:
# 	ax1d.set_xlabel('Normalized Time')
# ax1d.set_ylabel('Muscle Lengths (in mm)')
# [j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax1d.lines)]
# ax1d.legend(OrderedMuscleList)


###################################################################################################

	# Reverse Direction

X_Reverse = reaching_task(Xi=Xf, Xf=Xi, dt=dt, t_end=t_end, EOM=EOM)
return_MA_matrix()

	# Calculate resulting muscle velocities, lengths, etc.

NormalizedMuscleVelocity_Reverse = calculate_muscle_velocities()
WeightedPotentialTorqueVariation_shoulder_Reverse,\
WeightedPotentialTorqueVariation_elbow_Reverse = \
calculate_weighted_unscaled_potential_torque_variations_REACHING_TASK()
WeightedNormalizedMuscleVelocity_Reverse = calculate_weighted_muscle_velocities_REACHING_TASK()

	# Calculate only the lengthening Torque Variation Contributions for both the shoulder and the elbow.

EccentricTorqueVariations_Shoulder_Reverse = \
	ScalingFactor*concentric_velocities(WeightedPotentialTorqueVariation_shoulder_Reverse)
EccentricTorqueVariations_Elbow_Reverse =\
 	ScalingFactor*concentric_velocities(WeightedPotentialTorqueVariation_elbow_Reverse)
PotentialVariability_Reverse = calculate_potential_variability(X_Reverse,\
									EccentricTorqueVariations_Shoulder_Reverse,\
										EccentricTorqueVariations_Elbow_Reverse,\
											dt=dt, EOM=EOM, scheme = "Total")


	# Test to make sure that the muscle velocities are negative, time-reversed versions of each other for the forward and backwards movement.

assert np.array([np.array(abs(NormalizedMuscleVelocity_Forward[j,0,:].T-np.array(list(reversed(-NormalizedMuscleVelocity_Reverse[j,0,:]))).T)<1e-12).all() for j in range(n_muscles)]).all(), "The muscle velocities are not reversed and negative versions of each other."

	# Plot Afferent-Weighted Muscle Velocity (Reverse)

fig2e = plt.figure()
[plt.plot(t.T,WeightedNormalizedMuscleVelocity_Reverse[i].T) for i in OrderNumber]
ax2e = plt.gca()
ax2e.set_xlim(0,t_end*(1.3))
ax2e.set_ylim(-12,12)
ax2e.set_title(DescriptiveTitle+'\nAfferent-Weighted Normalized Muscle Velocity')
if t_end!=1:
	ax2e.set_xlabel('Time (s)')
else:
	ax2e.set_xlabel('Normalized Time')
ax2e.set_ylabel('Afferent-Weighted $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
[j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax2e.lines)]
ax2e.legend(OrderedMuscleList)

	# Plot Normalized Muscle Velocity (Reverse)

fig2a = plt.figure()
[plt.plot(t.T,NormalizedMuscleVelocity_Reverse[i].T) for i in OrderNumber]
ax2a = plt.gca()
ax2a.set_xlim(0,t_end*(1.3))
ax2a.set_ylim(-1.5,1.5)
ax2a.set_title(DescriptiveTitle+'\nNormalized Muscle Velocity')
if t_end!=1:
	ax2a.set_xlabel('Time (s)')
else:
	ax2a.set_xlabel('Normalized Time')
ax2a.set_ylabel('Normalized Muscle Velocity\nConcentric $\longleftrightarrow$ Eccentric')
[j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax2a.lines)]
ax2a.legend(OrderedMuscleList)

	# Plot Joint Angles (Reverse)

fig2b = plt.figure()
plt.plot(t.T,A1.T,'b')
plt.plot(t.T,A2.T,'b:')
ax2b = plt.gca()
ax2b.set_xlim(0,t_end*(1.3))
ax2b.set_title(DescriptiveTitle+'\nAngle vs. Time')
if t_end!=1:
	ax2b.set_xlabel('Time (s)')
else:
	ax2b.set_xlabel('Normalized Time')
ax2b.set_ylabel('Joint Angles (in radians)')
ax2b.legend(["Shoulder","Elbow"])

	# Plot Configuration Space (Reverse)

fig2c = plt.figure()
plt.plot(A1.T,A2.T,'b')
ax2c = plt.gca()
ax2c.set_title(DescriptiveTitle+'\nConfiguration Space')
ax2c.set_ylabel("Shouler Angle (in radians)")
ax2c.set_xlabel("Elbow Angle (in radians)")

	# Calculate Costs from Afferent-Weighted Muscle Velocities

EccentricCost_Reverse = eccentric_cost(WeightedNormalizedMuscleVelocity_Reverse,t_end=t_end,dt = dt,costtype='l1norm')
ConcentricCost_Reverse = concentric_cost(WeightedNormalizedMuscleVelocity_Reverse,t_end=t_end,dt = dt,costtype='l1norm')

# # Plot Muscle Lengths (Reverse)
#
# fig2d = plt.figure()
# [plt.plot(t,MuscleLengths_Reverse[i]) for i in OrderNumber]
# ax2d = plt.gca()
# ax2d.set_xlim(0,t_end*(1.3))
# ax2d.set_title(DescriptiveTitle+'\nMuscle Lengths')
# if t_end!=1:
# 	ax2d.set_xlabel('Time (s)')
# else:
# if t_end!=1:
# 	ax2d.set_xlabel('Time (s)')
# else:
# 	ax2d.set_xlabel('Normalized Time')
# ax2d.set_ylabel('Muscle Lengths (in mm)')
# [j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax2d.lines)]
# ax2d.legend(OrderedMuscleList)

###################################################################################################

	# Plot bar graph comparing the two directions

ConcentricCost_Forward = PotentialVariability_Forward.mean()
ConcentricCost_Reverse = PotentialVariability_Reverse.mean()
fig2 = plt.figure()
plt.bar(np.arange(2),[ConcentricCost_Forward,ConcentricCost_Reverse])
ax4 = plt.gca()
ax4.set_xticklabels(('Forward','Reverse'))
ax4.set_ylim(0,1.1*max([ConcentricCost_Forward, ConcentricCost_Reverse]))
ax4.set_yticks([0,max([ConcentricCost_Forward, ConcentricCost_Reverse])])
# ax4.set_ylim(0,20)
# ax4.set_yticks([0,10,20])
# ax4.set_yticklabels(['0','','20'])
# ax4.set_ylim(0,0.008)
# ax4.set_yticks([0,0.004,0.008])
ax4.set_yticklabels(['0','','0.008'])
ax4.set_title(DescriptiveTitle + '\nEccentric Cost for Forward and Reverse Movements')
ax4.set_ylabel('Sum of Afferent-Weighted Muscle Lengthening')
ax4.set_xticks([0,1])

###################################################################################################

	# Plot potential variability comparing the two directions
percent_cycle = np.concatenate((\
				(((X_Forward-X_Forward[:,0][:,np.newaxis])**2).sum(axis=0)**0.5)/0.35,\
					1+(((X_Reverse-X_Reverse[:,0][:,np.newaxis])**2).sum(axis=0)**0.5)/0.35),\
						axis=0)
double_t = np.array(np.arange(0,2+2*dt,dt,dtype = 'float64')*t_end, ndmin = 2) # returns a (1,2⋅(1+dt)/dt) array of ndmin = 2.
fig3 = plt.figure()
plt.plot(t[0,1:],PotentialVariability_Forward[0,1:],'k')
plt.plot(1+dt+t[0,:-1],PotentialVariability_Reverse[0,1:],'k')
# plt.plot(percent_cycle.T,np.concatenate((PotentialVariability_Forward,PotentialVariability_Reverse),\
# 						axis=1).T,'k',lw=2)
plt.plot([1,1],[0,1.1*max([max(PotentialVariability_Forward[0,1:]),max(PotentialVariability_Reverse[0,1:])])],color ='0.75',linestyle='--')
ax5 = plt.gca()
ax5.set_xticks([0,0.5,1,1.5,2])
ax5.set_xticklabels(('0','Forward','1','Reverse','2'))
# ax5.set_yticks([0,0.02])
ax5.set_yticks([0,1.1*max([max(PotentialVariability_Forward[0,1:]),max(PotentialVariability_Reverse[0,1:])])])
ax5.set_ylim([0,1.1*max([max(PotentialVariability_Forward[0,1:]),max(PotentialVariability_Reverse[0,1:])])])
ax5.set_title(DescriptiveTitle + '\nPotential For Endpoint Variability')
ax5.set_ylabel('Potential Endpoint Variability\nat Each Timestep')

	# Need to write PDF files and close the PDF files/figures

if SaveOutputFigures == True:
	pdf_forward.savefig(fig1e)
	pdf_forward.savefig(fig1a)
	pdf_forward.savefig(fig1b)
	pdf_forward.savefig(fig1c)
	# pdf_forward.savefig(fig1d)
	pdf_reverse.savefig(fig2e)
	pdf_reverse.savefig(fig2a)
	pdf_reverse.savefig(fig2b)
	pdf_reverse.savefig(fig2c)
	# pdf_reverse.savefig(fig2d)
	pdf_compare.savefig(fig2)
	pdf_compare.savefig(fig3)

	pdf_forward.close()
	pdf_reverse.close()
	pdf_compare.close()

		# Close any open plots

	plt.close('all')

plt.show()
