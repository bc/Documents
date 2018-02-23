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
def set_link_lengths(New_L1=None,New_L2=None,EOM = "Zadravec"):
	"""
	Sets the global link lengths for a 2 DOF planar reaching task. Changes the values of the link lengths. New_L1 and New_L2 must be numbers. Set EOM to Uno of Zadravec.

	Default L1 and L2 values are calculated from Winter's Anthropomorphic measurement scales for a 72 inch tall subject. L2 is currently the sum of forearm and hand measurements with some added value to reflect the reaching apparatus.
	"""
	if New_L1 != None:
		assert type(New_L1) == int or type(New_L1)==float, "New_L1 must be a number."
	if New_L2 != None:
		assert type(New_L2) == int or type(New_L2)==float, "New_L2 must be a number."
	assert EOM in [None,"Uno","Zadravec"], "EOM can be either None, 'Uno', or 'Zadravec'"
	global L1, L2
	# 72
	Height_inches = 64.5
	Height = 2.54*Height_inches/100
	if New_L1 == None and New_L2 == None:
		if EOM == "Uno":
			L1 = 0.256
			L2 = 0.315
		elif EOM == "Zadravec":
			L1 = 0.298
			L2 = 0.419
		else:
			L1 = Height*0.186
			L2 = Height*(0.146+0.108)
	elif New_L1 != None:
		L1 = New_L1
		L2 = Height*(0.146+0.108)
	elif New_L2 != None:
		L1 = Height*0.186
		L2 = New_L2
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
def c_matrix(x1,x2,x3):
	"""
	Takes in the values of x1, x2, and x3 to create the C matrix needed to find the coefficients of a clamped
	cubic spline with only one break (i.e. Cx = y, where x is an array of c coefficients for the
	piecewise polynomial equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). Returns a matrix.
	"""
	import numpy as np
	C = np.array([	[	2*(x2-x1), 		(x2-x1), 			0			],   \
					[	(x2-x1), 		2*(x3-x1), 		(x3-x2)		],   \
					[	0,				(x3-x2),		2*(x3-x2)	] 	], \
					float)
	return(C)
def y_vector(x1,x2,x3,y1,y2,y3,initial_slope,final_slope):
	"""
	Takes in the values of (x1,y1), (x2,y2), and (x3,y3) to create the y array necessary for the clamped cubic
	spline matrix manipulation for one break only (i.e. Cx = y, where x is an array of c coefficients for the
	piecewise polynomial equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). Returns an array.
	"""
	import numpy as np
	y = np.array([	3*(y2-y1)/(x2-x1) - 3*initial_slope ,  	\
					3*(y3-y2)/(x3-x2) - 3*(y2-y1)/(x2-x1),  \
					3*final_slope - 3*(y3-y2)/(x3-x2)	],  \
					float)
	return(y)
def c_coefficients(x1,x2,x3,y1,y2,y3,initial_slope,final_slope):
	"""
	Using matrix manipulations the equation Cx = y necessary for the c coefficients for a clamped cubic spline
	with only one break (i.e. Cx = y, where x is an array of c coefficients for the piecewise polynomial
	equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3) can be rearranged such that x = C.T*y. The values
	(x1,y1), (x2,y2), and (x3,y3) are the three points needed to the spline and initial_slope and final_slope
	are the endpoint conditions. Returns an array.
	"""
	import numpy as np

	C = c_matrix(x1,x2,x3)
	y = y_vector(x1,x2,x3,y1,y2,y3,initial_slope,final_slope)
	CCoefficients = (np.matrix(C)**(-1))*(np.matrix(y).T)
	return(CCoefficients)
def d_coefficients(x1,x2,x3,CCoefficients):
	"""
	Uses the c coefficients and the values of x1, x2, and x3 to find the d coefficients for the	piecewise
	polynomial equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). CCoefficients must be an array with
	three elements. Returns an array.
	"""
	import numpy as np
	DCoefficients = np.array([	(CCoefficients[1]-CCoefficients[0])/(3*(x2-x1)),  \
								(CCoefficients[2]-CCoefficients[1])/(3*(x3-x2))	],  \
								float)
	return(DCoefficients)
def b_coefficients(x1,x2,x3,y1,y2,y3,CCoefficients,DCoefficients):
	"""
	Uses the c and d coefficients and the values of (x1,y1), (x2,y2), and (x3,y3) to find the b coefficients for
	the	piecewise polynomial equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). CCoefficients must be an
	array with two or more elements and DCoefficients must be an array with two elements. Returns an array.
	"""
	import numpy as np
	BCoefficients = np.array([	((y2-y1)/(x2-x1)-CCoefficients[0]*(x2-x1) - DCoefficients[0]*((x2-x1)**2)),  \
								((y3-y2)/(x3-x2)-CCoefficients[1]*(x3-x2) - DCoefficients[1]*((x3-x2)**2)) 	]).astype(float)
	return(BCoefficients)
def test_b_coefficients(x1,x2,x3,y1,y2,y3,CCoefficients,DCoefficients,expected_slope):
	"""
	Tests to make sure that the generated b coefficients match the expected slope. Uses the c and d coefficients
	and the values of (x1,y1), (x2,y2), and (x3,y3) to find the b coefficients for the	piecewise polynomial
	equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). CCoefficients must be an array with two or more
	elements and DCoefficients must be an array with two elements. Returns TRUE if expected_slope equals b.
	"""
	import numpy as np
	B = b_coefficients(x1,x2,x3,y1,y2,y3,CCoefficients,DCoefficients)
	result = abs(B[0]-expected_slope)< 0.001
	return(result)
	assert B[0]==expected_slope, "First b coefficient (%f) does not equal initial slope (%f)." (B[0],expected_slope)
def a_coefficients(y1,y2):
	"""
	Uses the y values of (x1,y1) and (x2,y2) to find the a coefficients for the	piecewise polynomial equation
	y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3). Returns an array.
	"""
	import numpy as np
	ACoefficients = np.array([	y1,    \
								y2  ]).astype(float)
	return(ACoefficients)
def spline_coefficients(x1,x2,x3,y1,y2,y3,initial_slope,final_slope):
	"""
	Uses the values of (x1,y1), (x2,y2), and (x3,y3) to find the coefficients for the piecewise polynomial
	equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3) for a clamped cubic spline with one break only.
	Returns coefficient arrays A, B, C,and D.
	"""
	import numpy as np
	C = c_coefficients(x1,x2,x3,y1,y2,y3,initial_slope,final_slope)
	D = d_coefficients(x1,x2,x3,C)
	B = b_coefficients(x1,x2,x3,y1,y2,y3,C,D)
	A = a_coefficients(y1,y2)
	return(A,B,C[:2],D)
def generate_random_cartesian_point(xmin,xmax,ymin,ymax,distr = 'Normal'):
	"""
	Generates a random point in Cartesian space such that x is in [xmin, xmax] and y is in [ymin,ymax].
	Number is chosen at random from a uniform distribution. Returns two floats.
	"""
	assert type(distr)==str,"distr must be a string. Either 'Normal' or 'Uniform'."
	assert distr in ['Normal','Uniform'],"distr must be either 'Normal' or 'Uniform'."
	import numpy as np
	np.random.seed()
	if distr == 'Uniform':
		x_rand = np.random.uniform(xmin,xmax)
		y_rand = np.random.uniform(ymin,ymax)
	if distr == 'Normal':
		x_rand = np.random.normal(abs(xmin-xmax)/2,0.1,1000)
		y_rand = np.random.normal(ymin,ymax)
	return(x_rand,y_rand)
def test_endpoint_slope(b,c,d,x_n_minus_1,x_n,expected_slope):
	"""
	Takes in the cubic spline coefficients for the derivative of y = a + b*(x-x_n_minus_1) + c*(x-x_n_minus_1)**2 + d*(x-x_n_minus_1)**3
	(y' = b + 2*c*(x-x_n_minus_1) + 3*d*(x-x_n_minus_1)**2)	for the last piecewise polynomial and tests to see if the expected slope at
	the endpoint is equal to the actual	endpoint slope. The variable x_n_minus_1 is the initial value of the final piecewise polynomial
	and x_n is the final data point. Returns TRUE if they are equal.

	"""
	actual_slope = b + 2*c*(x_n-x_n_minus_1) + 3*d*(x_n-x_n_minus_1)**2
	result = abs(actual_slope-expected_slope)<0.001
	return(result)
def test_for_discontinuity(a_n,b_n,c_n,d_n,x_n,x_n_plus_1,y_n_plus_1):
	"""
	Takes in the coefficients for a cubic spline polynomial y = a_n + b_n*(x-x_n) + c_n*(x-x_n)**2 + d_n*(x-x_n)**3
	and tests to see if the final y value for this piecewise polynomial is equal to the initial y value of the next
	piecewise polynomial (i.e. when x = x_n_plus_1). The variable x_n is the initial x value of the preceding
	polynomial, and x_n_plus_1 is the transition value from one polynomial to the next. y_n_plus_1 is the initial y
	value for the next piecewise polynomial.
	"""
	y_n_final = a_n + b_n*(x_n_plus_1-x_n) + c_n*(x_n_plus_1-x_n)**2 + d_n*(x_n_plus_1-x_n)**3
	result = abs(y_n_final-y_n_plus_1)<0.001
	return(result)
class Spline:
	"""
	Initiate a class variable spline that has one break at x = x_break starting at x_initial and has
	the equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3.

	pp_func()
	~~~~~~~~~~~~~~~~~~~

	Takes in X array and outputs the piecewise polynomial associated with this spline.

	pp_deriv()
	~~~~~~~~~~~~~~~~~~~

	Takes in X array and outputs the piecewise polynomial associated with the spline's derivative.

	find_max_and_min()
	~~~~~~~~~~~~~~~~~~~

	Takes in the min and max values for both x and y and will find the maximum y values of the piecewise
	polynomial. To do this, first we find the extrema point (find_extrema) by inputing the x values that
	set the derivate of the piecewise polynomial equal to zero (quadratic formula). Next we ensure that
	the zero values are in fact real (is_real). We then filter out the zeros that are not in the
	appropriate domains (is_in_appropriate_domain). To see if these values are maximum or minimum, we
	plug them back into the second derivative of the appropriate piecewise polynomial (second_deriv_is_neg()
	and second_deriv_is_pos(), respectively). Finally we determine the y value of these extrema by using
	the class function self.pp_func().

	is_initial_slope_positive()
	~~~~~~~~~~~~~~~~~~~

	This takes in X and will check to make sure that for the first 2500 entries in X, that the derivative
	of the piecewise polynomial (pp_deriv()) will be positive. Make sure that X is at least 2500 in length.

	is_within_bounds()
	~~~~~~~~~~~~~~~~~~~

	This checks to see if the maximum maximum value and the minimum mininum value calculated above will
	fall between y_min and y_max. This makes use of self.find_max_and_min()

	"""
	def __init__(self,a,b,c,d,x_initial,x_break,x_final):
		self.a = a
		self.b = b
		self.c = c
		self.d = d
		self.x_initial = x_initial
		self.x_break = x_break
		self.xlim = [x_initial,x_final]
		#self.all_values = {'A': a, 'B' : b, 'C' : c, 'D' : d, 'init' : x_initial, 'break' : x_break}
	def pp_func(self,X):
		result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
			[lambda X: self.a[0] + self.b[0,0,0]*(X-self.x_initial) + self.c[0,0]*(X-self.x_initial)**2 + self.d[0,0,0]*(X-self.x_initial)**3, \
			lambda X: self.a[1] + self.b[1,0,0]*(X-self.x_break) + self.c[1,0]*(X-self.x_break)**2 + self.d[1,0,0]*(X-self.x_break)**3])
		return(result)
	def pp_deriv(self,X):
		result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
			[lambda X: self.b[0,0,0] + 2*self.c[0,0]*(X-self.x_initial) + 3*self.d[0,0,0]*(X-self.x_initial)**2, \
			lambda X: self.b[1,0,0] + 2*self.c[1,0]*(X-self.x_break) + 3*self.d[1,0,0]*(X-self.x_break)**2])
		return(result)
	def pp_2deriv(self,X):
		result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
			[lambda X: 2*self.c[0,0] + 6*self.d[0,0,0]*(X-self.x_initial), \
			lambda X: 2*self.c[1,0] + 6*self.d[1,0,0]*(X-self.x_break)])
		return(result)

	def find_max_and_min(self,x_min,x_max,y_min,y_max):
		def find_extrema():
			if (4*self.c[0,0]**2 - 12*self.b[0,0,0]*self.d[0,0,0]) >= 0:
				extrema_1 = self.x_initial + (- 2*self.c[0,0] + (4*self.c[0,0]**2 - 12*self.b[0,0,0]*self.d[0,0,0])**.5)/(6*self.d[0,0,0])
				extrema_2 = self.x_initial + (- 2*self.c[0,0] - (4*self.c[0,0]**2 - 12*self.b[0,0,0]*self.d[0,0,0])**.5)/(6*self.d[0,0,0])
			else:
				extrema_1, extrema_2 = None, None
			if (4*self.c[1,0]**2 - 12*self.b[1,0,0]*self.d[1,0,0]) >= 0:
				extrema_3 = self.x_break + (- 2*self.c[1,0] + (4*self.c[1,0]**2 - 12*self.b[1,0,0]*self.d[1,0,0])**.5)/(6*self.d[1,0,0])
				extrema_4 = self.x_break + (- 2*self.c[1,0] - np.sqrt(4*self.c[1,0]**2 - 12*self.b[1,0,0]*self.d[1,0,0]))/(6*self.d[1,0,0])
			else:
				extrema_3, extrema_4 = None, None
			return(extrema_1,extrema_2,extrema_3,extrema_4)
		def is_real(x_value):
			result = not isinstance(x_value,complex)
			return(result)
		def is_in_appropriate_domain(x_value,x_min,x_max,segment_number):
			if segment_number == 1:
				result = x_value >= x_min and x_value <= self.x_break
			elif segment_number == 2:
				result = x_value >= self.x_break and x_value <= x_max
			return(result)
		def second_deriv_is_neg(x_value,segment_number):
			if segment_number == 1:
				x_not = self.x_initial
			elif segment_number == 2:
				x_not = self.x_break
			second_deriv =  2*self.c[segment_number-1] + 6*self.d[segment_number-1]*(x_value-x_not)
			result = second_deriv<0
			return(result)
		def second_deriv_is_pos(x_value,segment_number):
			if segment_number == 1:
				x_not = self.x_initial
			elif segment_number == 2:
				x_not = self.x_break
			second_deriv =  2*self.c[segment_number-1] + 6*self.d[segment_number-1]*(x_value-x_not)
			result = second_deriv>0
			return(result)
		def determine_if_max_or_min(extrema_1,extrema_2,extrema_3,extrema_4,x_min,x_max):
			maxima = []
			minima = []
			if extrema_1 != None and is_in_appropriate_domain(extrema_1,x_min,x_max,1):
				if second_deriv_is_neg(extrema_1,1):
					maxima.append(np.float(self.pp_func(extrema_1)))
				elif second_deriv_is_pos(extrema_1,1):
					minima.append(np.float(self.pp_func(extrema_1)))
			if extrema_2 != None and is_in_appropriate_domain(extrema_2,x_min,x_max,1):
				if second_deriv_is_neg(extrema_2,1):
					maxima.append(np.float(self.pp_func(extrema_2)))
				elif second_deriv_is_pos(extrema_2,1):
					minima.append(np.float(self.pp_func(extrema_2)))
			if extrema_3 != None and is_in_appropriate_domain(extrema_3,x_min,x_max,2):
				if second_deriv_is_neg(extrema_3,2):
					maxima.append(np.float(self.pp_func(extrema_3)))
				elif second_deriv_is_pos(extrema_3,2):
					minima.append(np.float(self.pp_func(extrema_3)))
			if extrema_4 != None and is_in_appropriate_domain(extrema_4,x_min,x_max,2):
				if second_deriv_is_neg(extrema_4,2):
					maxima.append(np.float(self.pp_func(extrema_4)))
				elif second_deriv_is_pos(extrema_4,2):
					minima.append(np.float(self.pp_func(extrema_4)))
			return(maxima,minima)
		extrema_1,extrema_2,extrema_3,extrema_4 = find_extrema()
		maxima, minima = determine_if_max_or_min(extrema_1,extrema_2,extrema_3,extrema_4,x_min,x_max)
		return(maxima,minima)
	def is_initial_slope_positive(self,X,cutoff):
		result = min(self.pp_deriv(X[:cutoff]))>=0
		return(result)
	def is_within_bounds(self,x_min,x_max,y_min,y_max):
		maxima,minima = self.find_max_and_min(x_min,x_max,y_min,y_max)
		if len(maxima) == 0:
			maxima = y_max
		if len(minima) == 0:
			minima = y_min
		result = np.max(maxima) <= y_max and np.min(minima) >= y_min
		return(result)
	def print_func(self):
		from sympy import Symbol,Lambda,pprint
		x = Symbol('x')
		func_1 = Lambda(x,self.a[0] + self.b[0,0,0]*(x-self.x_initial) + self.c[0,0]*(x-self.x_initial)**2 + self.d[0,0,0]*(x-self.x_initial)**3)
		print('Function 1:\n')
		pprint(func_1)
		func_2 = Lambda(x,self.a[1] + self.b[1,0,0]*(x-self.x_break) + self.c[1,0]*(x-self.x_break)**2 + self.d[1,0,0]*(x-self.x_break)**3)
		print('Function 2:\n')
		pprint(func_2)
	def return_parameterized_X(self):
		import scipy.integrate as integrate
		dt = 1/1000
		t = np.linspace(0,1,int(1/dt) + 1)
		def ode_func(x,t):
			return(self.dS_dt(t)/np.sqrt(1 + self.pp_deriv(x)**2))
		X = integrate.odeint(ode_func,self.xlim[0],t).flatten()
		Y = np.array(list(map(lambda x: self.pp_func(x),X)))
		dS = np.array(list(map(lambda dx,dy: np.sqrt(dx**2+dy**2),\
								np.gradient(X)/dt,np.gradient(Y)/dt)))
		assert sum(abs(self.dS_dt(t)-dS))/len(dS)<1e-4, "Error in parameterizing path to dS/dt. Check ODE func."
		return(X,Y)
	def return_parameterized_dX(self,x):
		dt = 1/1000
		t = np.linspace(0,1,int(1/dt) + 1)
		dS_dt = self.dS_dt(t)
		df_dx = self.pp_deriv(x)
		dx_dt = np.array(list(map(lambda dS_dt,df_dx: dS_dt/np.sqrt(1 + df_dx**2),dS_dt,df_dx)))
		dy_dt = df_dx*dx_dt
		return(dx_dt,dy_dt)
	def return_parameterized_d2X(self,x,dx_dt):
		dt = 1/1000
		t = np.linspace(0,1,int(1/dt) + 1)

		dS_dt = self.dS_dt(t)
		d2S_dt2 = self.d2S_dt2(t)

		df_dx = self.pp_deriv(x)
		d2f_dx2 = self.pp_2deriv(x)

		d2x_dt2 = np.array(list(map(lambda dx_dt,d2S_dt2,df_dx,d2f_dx2: \
							(d2S_dt2*np.sqrt(1+df_dx**2) - df_dx*d2f_dx2*(dx_dt**2))\
									/(1 + df_dx**2), \
								dx_dt,d2S_dt2,df_dx,d2f_dx2)))
		d2y_dt2 = np.array(list(map(lambda d2f_dx2,dx_dt,df_dx,d2x_dt2:  \
							d2f_dx2*(dx_dt**2) + df_dx*(d2x_dt2),\
								d2f_dx2,dx_dt,df_dx,d2x_dt2)))
		return(d2x_dt2,d2y_dt2)
	def find_path_length(self):
		import scipy.integrate as integrate
		return(integrate.quad(lambda x: \
				np.sqrt(1 + self.pp_deriv(x)**2),self.xlim[0],self.xlim[1])[0])
	def dS_dt(self,t):
		S_initial = 0
		S_final = self.find_path_length()
		return((S_final-S_initial)*(30*t**2 - 60*t**3 + 30*t**4))
	def d2S_dt2(self,t):
		S_initial = 0
		S_final = self.find_path_length()
		return((S_final-S_initial)*(60*t - 180*t**2 + 120*t**3))
def clamped_cubic_spline(x_initial,x_final,y_initial,y_final,initial_slope,final_slope,ymin,ymax,X,**options):
	"""
	This will take in the initial and final values for both x and y, as well as the desired initial and final
	slopes and generate 1000 clamped cubic spline that produce y values that are within the bounds [ymin, ymax].
	Returns a list of arrays, each of len(X). Options allows for slope limitations on shoulder rotations such
	that the derivative of the spline is always positive to match observations (slope = "Shoulder").
	"""
	NumberOfTrials =  1000
	i = 0
	StartTime = time.time()
	while i < NumberOfTrials:
		x_rand,y_rand = generate_random_point(x_initial,x_final,ymin,ymax)
		A,B,C,D = spline_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,initial_slope,final_slope)
		assert test_b_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,C,D,initial_slope), "Initial slope does not match the expected value"
		assert test_endpoint_slope(B[1],C[1],D[1],x_rand,x_final,final_slope),"Problem with Endpoint Slope"
		assert test_for_discontinuity(A[0],B[0],C[0],D[0],x_initial,x_rand,A[1]), "Jump Discontinuity at t = %f!" %x_rand
		spline_structure = Spline(A,B,C,D,x_initial,x_rand,x_final)
		if options.get("angle") == "Shoulder":
			options_slope_condition = spline_structure.is_initial_slope_positive(X,2501)
			dof = "Shoulder: " + " "*(10-len("Shoulder: "))
		elif options.get("angle") == "Elbow":
			options_slope_condition = spline_structure.is_initial_slope_positive(X,501)
			dof = "Elbow: " + " "*(10-len("Elbow: "))
		else:
			options_slope_condition = True
			dof = "Wrist: " + " "*(10-len("Wrist: "))
		if i == 0:
			if spline_structure.is_within_bounds(x_initial,x_final, ymin, ymax) and options_slope_condition:
				Splines = spline_structure
				i+=1
				statusbar = '[' + '\u25a0'*int((i+1)/(NumberOfTrials/50)) + '\u25a1'*(50-int((i+1)/(NumberOfTrials/50))) + '] '
				print(dof + statusbar + '{0:1.1f}'.format(i/NumberOfTrials*100) + '% complete, ' + '{0:1.1f}'.format(time.time() - StartTime) + 'sec        \r', end='')
		elif i == 1:
			if spline_structure.is_within_bounds(x_initial,x_final, ymin, ymax) and options_slope_condition:
				Splines = np.concatenate(([Splines], [spline_structure]), axis = 0)
				i+=1
				statusbar = '[' + '\u25a0'*int((i+1)/(NumberOfTrials/50)) + '\u25a1'*(50-int((i+1)/(NumberOfTrials/50))) + '] '
				print(dof + statusbar + '{0:1.1f}'.format(i/NumberOfTrials*100) + '% complete, ' + '{0:1.1f}'.format(time.time() - StartTime) + 'sec        \r', end='')
		else:
			if spline_structure.is_within_bounds(x_initial,x_final, ymin, ymax) and options_slope_condition:
				Splines = np.concatenate((Splines, [spline_structure]), axis = 0)
				i+=1
				statusbar = '[' + '\u25a0'*int((i+1)/(NumberOfTrials/50)) + '\u25a1'*(50-int((i+1)/(NumberOfTrials/50))) + '] '
				print(dof + statusbar + '{0:1.1f}'.format(i/NumberOfTrials*100) + '% complete, ' + '{0:1.1f}'.format(time.time() - StartTime) + 'sec        \r', end='')
	print('\n')
	return(Splines)
def generate_default_path(PathLength = 0.35,RandomXi=False,RandomXf=False):
	ValidPath = False
	while ValidPath==False:
		x_initial_desired = 0.05 # m DISPLACEMENT TO BE TRANSLATED BACK LATER
		x_final_desired = PathLength + 0.05 # m
		y_initial_desired = 0 # m
		y_final_desired = 0 # m
		boundary_sigma = 0.0025
		allowable_error_rad = 30*(np.pi/180) #degrees
		allowable_error_y = 0.05 # m

		if RandomXi == False:
			x_initial = x_initial_desired + np.random.normal(0,boundary_sigma) # cm
			y_initial = y_initial_desired + np.random.normal(0,boundary_sigma) # cm
		else:
			x_initial = x_initial_desired # cm
			y_initial = y_initial_desired # cm

		if RandomXf == False:
			x_final = x_final_desired + np.random.normal(0,boundary_sigma) # cm
			y_final = y_final_desired + np.random.normal(0,boundary_sigma) # cm
		else:
			x_final = x_final_desired # cm
			y_final = y_final_desired # cm

		initialerror = np.random.uniform(-np.tan(allowable_error_rad),np.tan(allowable_error_rad))
		finalerror = -np.sign(initialerror)*\
				abs(np.random.uniform(-np.tan(allowable_error_rad),np.tan(allowable_error_rad)))

		if initialerror>0:
			ymax = max([y_initial,y_final])+allowable_error_y
			ymin = 0
		else:
			ymax = 0
			ymin = min([y_initial,y_final])-allowable_error_y

		xmax = x_final
		xmin = x_initial
		x_rand = np.random.normal((x_initial+x_final)/2,abs(x_initial+x_final)/4)
		y_rand = np.random.normal((ymax+ymin)/2,abs(ymax+ymin)/4)

		A,B,C,D = spline_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,initialerror,finalerror)
		assert test_b_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,C,D,initialerror), "Initial slope does not match the expected value"
		assert test_endpoint_slope(B[1],C[1],D[1],x_rand,x_final,finalerror),"Problem with Endpoint Slope"
		assert test_for_discontinuity(A[0],B[0],C[0],D[0],x_initial,x_rand,A[1]), "Jump Discontinuity at t = %f!" %x_rand
		path = Spline(A,B,C,D,x_initial,x_rand,x_final)
		if path.is_within_bounds(x_initial,x_final, ymin, ymax):
			ValidPath = True
		else:
			ValidPath = False
	return(path)
def rotate_xy(x,y,angle):
	XY = np.concatenate([x[np.newaxis,:],y[np.newaxis,:]],axis=0)
	rotation_matrix = np.matrix([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
	rotated_XY = rotation_matrix*XY
	rotated_x = np.array(rotated_XY[0,:])[0]
	rotated_y = np.array(rotated_XY[1,:])[0]
	return(rotated_x,rotated_y)
def translate_xy(x,y,px=0,py=0):
	translated_x = x + px
	translated_y = y + py
	return(translated_x,translated_y)
def return_trajectory(x,y,ReachType,IsPositionFunction):
	"""
	This takes in a string -- either 'Center','Right','Left', or 'Sideways' -- and returns the necessary initial and final positions for the movement, based on Flash/Hogan (1985).

	Parameters:
	# 0.80*(L1+L2) = 0.4586
	# 0.35/(L1+L2) = 0.5295
	# Sternum at -0.177 = -L1*(0.129/0.186) <-- Anthropomorphic Ratio
	"""
	global L1,L2
	PathLength = 0.35
	MedianPlane = -L1*(0.129/0.186) # Sternum at -L1*(0.129/0.186) <-- Anthropomorphic Ratio
	DefaultDisplacement_x = 0.05
	assert ReachType.capitalize() in ['Center','Right','Left','Sideways'], \
		"ReachType must be either 'Center','Right','Left', or 'Sideways'."
	assert type(IsPositionFunction)==bool, "IsPositionFunction is either True or False."

	if IsPositionFunction==True:
		# Side to Side, Path Length = 0.5295⋅(L1+L2) = 0.35
		if ReachType.capitalize() == 'Sideways':
			x_Trajectory,y_Trajectory = translate_xy(x,y,px=-DefaultDisplacement_x)
			# x_Trajectory, y_Trajectory = rotate_xy(x_Trajectory,y_Trajectory,0)
			x_Trajectory, y_Trajectory = translate_xy(x_Trajectory,y_Trajectory,\
												px = MedianPlane-PathLength/2,\
													py=0.10+PathLength/2)
			# Xi = [-MedianPlane-PathLength/2,0.1+PathLength/2]
			# Xf = [-MedianPlane+PathLength/2,0.1+PathLength/2]

		# Center Reach, Path Length = 0.5295⋅(L1+L2) = 0.35
		elif ReachType.capitalize() == 'Center':
			x_Trajectory,y_Trajectory = translate_xy(x,y,px=-DefaultDisplacement_x)
			x_Trajectory, y_Trajectory = rotate_xy(x_Trajectory,y_Trajectory,np.pi/2)
			x_Trajectory, y_Trajectory = translate_xy(x_Trajectory,y_Trajectory,\
												px = MedianPlane,\
													py=0.20)
			# Xi = [-MedianPlane,0.20]
			# Xf = [-MedianPlane,0.20 + PathLength]

		# # Left Diagonal Reach, Path Length = 0.5295⋅(L1+L2) = 0.35
		elif ReachType.capitalize() == 'Left':
			x_Trajectory,y_Trajectory = translate_xy(x,y,px=-DefaultDisplacement_x)
			x_Trajectory, y_Trajectory = rotate_xy(x_Trajectory,y_Trajectory,3*np.pi/4)
			x_Trajectory, y_Trajectory = translate_xy(x_Trajectory,y_Trajectory,\
												px = MedianPlane,\
													py=0.20)
			# Xi = [-MedianPlane,0.20]
			# Xf = [-MedianPlane-PathLength/(2**0.5),0.20 + PathLength/(2**0.5)]

		# # Right Diagonal Reach, Path Length = 0.5295⋅(L1+L2) = 0.35
		elif ReachType.capitalize() == 'Right':
			x_Trajectory,y_Trajectory = translate_xy(x,y,px=-DefaultDisplacement_x)
			x_Trajectory, y_Trajectory = rotate_xy(x_Trajectory,y_Trajectory,np.pi/4)
			x_Trajectory, y_Trajectory = translate_xy(x_Trajectory,y_Trajectory,\
												px = MedianPlane,\
													py=0.20)
			# Xi = [-MedianPlane,0.20]
			# Xf = [-MedianPlane+PathLength/(2**0.5),0.20 + PathLength/(2**0.5)]
	else:
		# Side to Side, Path Length = 0.5295⋅(L1+L2) = 0.35
		if ReachType.capitalize() == 'Sideways':
			x_Trajectory,y_Trajectory = x,y

		# Center Reach, Path Length = 0.5295⋅(L1+L2) = 0.35
		elif ReachType.capitalize() == 'Center':
			x_Trajectory, y_Trajectory = rotate_xy(x,y,np.pi/2)

		# # Left Diagonal Reach, Path Length = 0.5295⋅(L1+L2) = 0.35
		elif ReachType.capitalize() == 'Left':
			x_Trajectory, y_Trajectory = rotate_xy(x,y,3*np.pi/4)

		# # Right Diagonal Reach, Path Length = 0.5295⋅(L1+L2) = 0.35
		elif ReachType.capitalize() == 'Right':
			x_Trajectory, y_Trajectory = rotate_xy(x,y,np.pi/4)
	return(x_Trajectory,y_Trajectory)
def return_X_values(DefaultPath,ReachType,RandomXi=False,RandomXf=False,PathLength=0.35):
	# import ipdb; ipdb.set_trace()
	Default_x,Default_y = DefaultPath.return_parameterized_X()
	Default_ẋ,Default_ẏ = DefaultPath.return_parameterized_dX(Default_x)
	Default_ẍ,Default_ÿ = DefaultPath.return_parameterized_d2X(Default_x,Default_ẋ)

	x,y = return_trajectory(Default_x,Default_y,ReachType,IsPositionFunction=True)
	ẋ,ẏ = return_trajectory(Default_ẋ,Default_ẏ,ReachType,IsPositionFunction=False)
	ẍ,ÿ = return_trajectory(Default_ẍ,Default_ÿ,ReachType,IsPositionFunction=False)

	X = np.concatenate([x[np.newaxis,:],y[np.newaxis,:]],axis=0)
	Ẋ = np.concatenate([ẋ[np.newaxis,:],ẏ[np.newaxis,:]],axis=0)
	Ẍ = np.concatenate([ẍ[np.newaxis,:],ÿ[np.newaxis,:]],axis=0)
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
def reaching_task(ReachType='Center',RandomXi=False, RandomXf=False, \
					PathLength=0.35,EOM='Zadravec'):
	import numpy as np
	assert ReachType.capitalize() in ['Center','Right','Left','Sideways'], \
		"ReachType must be either 'Center','Right','Left', or 'Sideways'."
	set_link_lengths()
	create_angle_lists()
	# t = np.array(np.arange(0,1+dt,dt),dtype='float64',ndmin=2)
	DefaultPath = generate_default_path(PathLength=PathLength,RandomXi=RandomXi,RandomXf=RandomXf)
	X,Ẋ,Ẍ = return_X_values(DefaultPath,ReachType,RandomXi=RandomXi,\
								RandomXf=RandomXf,PathLength=PathLength)
	update_angle_lists(X,Ẋ,Ẍ)
	# calculate_torques(EOM=EOM)
	return(X,Ẋ,Ẍ)
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

	OptimalMuscleLength = np.array([AllCoefficients[key]['Optimal Muscle Length'] for key in AllCoefficients.keys()],dtype='float64')[:,np.newaxis,np.newaxis]
	NormalizedMuscleVelocity = MuscleVelocity/OptimalMuscleLength

	"""
	Next we need to multiply each muscle velocity, AT EACH TIMESTEP, by the appropriate scaling factor, which may be a function of joint angle (i.e., movement step number).
	"""

	CorrectedAfferentNumber = np.array([AllCoefficients[key]['Corrected No'] for key in AllCoefficients.keys()],dtype='float64')[:,np.newaxis] # returns an (m,1) array

	"""
	Note: (NormalizedMuscleVelocity*MomentArmMatrix.T/1000).T returns an (N,2,m) matrix. Therefore, we must construct a (2,m) matrix of Corrected Afferented Numbers that have identical rows. This is done by concatenating CorrectedAfferentNumber with itself and then transposing the result. This total product will yield a (m,2,N) array of both shoulder and elbow weighted muscle velocities. This will then be split by np.split() to give (m,1,N) arrays that are then squeezed into (m,N) arrays for the individual joint weighted muscle velocities.

	MomentArmMatrix must be rectified so that the contraction direction is maintained.
	"""

	MAWeightedMuscleVelocity = \
		((NormalizedMuscleVelocity*abs(MomentArmMatrix).T/1000).T\
			*np.concatenate((CorrectedAfferentNumber,CorrectedAfferentNumber),axis=1).T).T
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
def animate_plots(response,ReachType, Weighted=False, save_as_gif = False):
	assert type(response)==bool, "Input must be either True or False."
	assert type(Weighted)==bool, "Weighted must be either True or False."
	assert ReachType in ['Sideways','Center','Left','Right'], "ReachType must be either 'Sideways','Center','Left', or 'Right'"

	if response == True:
		import numpy as np
		import matplotlib.pyplot as plt
		from matplotlib.patches import Ellipse
		import matplotlib.animation as animation
		import matplotlib.patches as patches
		global A1,A2,L1,L2
		A1_forward = A1
		A2_forward = A2
		A1_reverse = np.array(list(reversed(A1.T))).T
		A2_reverse = np.array(list(reversed(A2.T))).T

		if Weighted == True:
			Vm_forward = WeightedNormalizedMuscleVelocity_Forward
			Vm_reverse = -np.array(list(reversed(WeightedNormalizedMuscleVelocity_Forward.T))).T
		else:
			Vm_forward = NormalizedMuscleVelocity_Forward
			Vm_reverse = -np.array(list(reversed(NormalizedMuscleVelocity_Forward.T))).T

		MedianPlane = L1*(0.129/0.186)

		fig, ((ax5,ax6),(ax1,ax2),(ax3,ax4)) = plt.subplots(3,2,figsize=(10,8))
		plt.subplots_adjust(top=0.9,hspace=0.2,bottom=0.2,left=0.2)

		if ReachType == 'Left':
			DescriptiveTitle = "45$^\circ$ Reach Left"
		elif ReachType == 'Right':
			DescriptiveTitle = "45$^\circ$ Reach Right"
		elif ReachType == 'Sideways':
			DescriptiveTitle = 'Side-to-side Reach'
		else:
			DescriptiveTitle = 'Straight Forward (Center) Reach'

		if Weighted==True:
			TypeString = "\n(Afferent-Weighted $\hat{v}_m$)\n"
		else:
			TypeString = "\n(Normalized $\hat{v}_m$)\n"
		plt.suptitle(DescriptiveTitle + TypeString,Fontsize=20,y=0.975)

		#Forward Model

		Angle1_f, = ax5.plot([0],[A1_forward.T[0]],color = '0.60')
		Angle2_f, = ax5.plot([0],[A2_forward.T[0]],color = '0.60',linestyle='--')
		ax5.set_xlim(0,t_end)
		ax5.set_xticks([0,t_end])
		ax5.set_xticklabels(['Start','Finish'])
		ax5.set_ylim(0,np.pi)
		ax5.set_yticks([0,np.pi/2,np.pi])
		ax5.set_yticklabels(['0',r'$\frac{\pi}{2}$','$\pi$'],fontsize=12)
		ax5.spines['right'].set_visible(False)
		ax5.spines['top'].set_visible(False)
		# ax5.set_ylabel('Joint Angles\n(in Radians)')
		ax5.legend(["Shoulder\nAngle","Elbow\nAngle"],loc='center right',bbox_to_anchor=(-0.1, 0.5))

		ax1.get_xaxis().set_ticks([])
		ax1.get_yaxis().set_ticks([])
		ax1.set_frame_on(True)
		RightShoulder_f = plt.Circle((0,0),radius=0.05,Color='#4682b4')
		ax1.add_patch(RightShoulder_f)
		LeftShoulder_f = plt.Circle((-MedianPlane*2,0),radius=0.05,Color='#4682b4')
		ax1.add_patch(LeftShoulder_f)
		Torso_f = plt.Rectangle((-MedianPlane*2,-0.05),MedianPlane*2,0.1,Color='#4682b4')
		ax1.add_patch(Torso_f)
		Head_f = Ellipse((-MedianPlane,0),0.2,0.225,FaceColor='w',EdgeColor='#4682b4',Linewidth=3)
		ax1.add_patch(Head_f)
		JointCoordinates_f = \
			np.concatenate([[np.cumsum([0,L1*np.cos(A1_forward[0,0]),L2*(0.146/(0.146+0.108))*np.cos(A1_forward[0,0]+A2_forward[0,0]),L2*(0.108/(0.146+0.108))*np.cos(A1_forward[0,0]+A2_forward[0,0])])],\
			[np.cumsum([0,L1*np.sin(A1_forward[0,0]),L2*(0.146/(0.146+0.108))*np.sin(A1_forward[0,0]+A2_forward[0,0]),L2*(0.108/(0.146+0.108))*np.sin(A1_forward[0,0]+A2_forward[0,0])])]],\
			axis=0)
		Elbow_f = plt.Circle((L1*np.cos(A1_forward[0,0]),L1*np.sin(A1_forward[0,0])),radius=0.03,color='#4682b4')
		Endpoint_f = plt.Circle((L1*np.cos(A1_forward[0,0])+L2*np.cos(A1_forward[0,0]+A2_forward[0,0]),\
							L1*np.sin(A1_forward[0,0])+L2*np.sin(A1_forward[0,0]+A2_forward[0,0])),\
							radius = 0.02,color='#4682b4')
		Wrist_f = plt.Circle((L1*np.cos(A1_forward[0,0])+L2*(0.146/(0.146+0.108))*np.cos(A1_forward[0,0]+A2_forward[0,0]),\
							L1*np.sin(A1_forward[0,0])+L2*(0.146/(0.146+0.108))*np.sin(A1_forward[0,0]+A2_forward[0,0])),\
							radius = 0.03,color='#4682b4')
		StickFigure_f, = ax1.plot(JointCoordinates_f[0,:],JointCoordinates_f[1,:],'ko-',LineWidth=2,MarkerFaceColor='k')
		movement_f, = ax1.plot(L1*np.cos(A1_forward[0,0])+L2*np.cos(A1_forward[0,0]+A2_forward[0,0]),\
		 					L1*np.sin(A1_forward[0,0])+L2*np.sin(A1_forward[0,0]+A2_forward[0,0]),\
							color='0.60',linestyle=':')
		UpperArm_f = plt.Rectangle((0.02*np.sin(A1_forward[0,0]),\
										-0.02*np.cos(A1_forward[0,0])),\
										L1, 0.04,\
										angle=A1_forward[0,0]*180/np.pi,color='#4682b4',animated=True)
		ax1.add_patch(UpperArm_f)
		Forearm_f = plt.Rectangle((L1*np.cos(A1_forward[0,0]) + 0.02*np.sin(A1_forward[0,0]+A2_forward[0,0]),\
										L1*np.sin(A1_forward[0,0]) -0.02*np.cos(A1_forward[0,0]+A2_forward[0,0])),\
										L2*(0.146/(0.146+0.108)), 0.04,\
										angle=(A1_forward[0,0]+A2_forward[0,0])*180/np.pi,color='#4682b4',animated=True)
		ax1.add_patch(Forearm_f)
		Hand_f = plt.Rectangle((L1*np.cos(A1_forward[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.cos(A1_forward[0,0]+A2_forward[0,0]) \
									+ 0.02*np.sin(A1_forward[0,0]+A2_forward[0,0]),\
							L1*np.sin(A1_forward[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.sin(A1_forward[0,0]+A2_forward[0,0]) \
									-0.02*np.cos(A1_forward[0,0]+A2_forward[0,0])),\
										L2*(0.108/(0.146+0.108)), 0.04,\
										angle=(A1_forward[0,0]+A2_forward[0,0])*180/np.pi,color='#4682b4',animated=True)
		ax1.add_patch(Hand_f)

		NormalizedVmPlots_forward = [ax3.plot(t.T,Vm_forward[j].T) for j in OrderNumber]
		ax3.set_xlim(0,t_end)
		ax3.set_ylim(-12,12)
		ax3.set_xticks([0,t_end])
		ax3.set_xticklabels(['Start','Finish'])
		if Weighted == True:
			ax3.set_ylabel('Afferent-Weighted $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
		else:
			ax3.set_ylabel('Normalized $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
		[k.set_color(OrderedColorsList[j]) for j,k in enumerate(ax3.lines)]
		Screen_f = plt.Rectangle((0,-13),t_end,26,color='w')
		ax3.add_patch(Screen_f)
		ax3.spines['right'].set_visible(False)
		ax3.spines['top'].set_visible(False)
		bounds = max([Vm_forward.max(),Vm_reverse.max()])
		ax3.set_ylim([-1.1*bounds,1.1*bounds])
		plt.figlegend([el[0] for el in NormalizedVmPlots_forward],OrderedMuscleList,loc='lower center',ncol=5,mode='expand')

		#Reverse Model

		Angle1_r, = ax6.plot([0],[A1_reverse.T[0]],color = '0.60')
		Angle2_r, = ax6.plot([0],[A2_reverse.T[0]],color = '0.60',linestyle='--')
		ax6.set_xlim(0,t_end)
		ax6.set_xticks([0,t_end])
		ax6.set_xticklabels(['Start','Finish'])
		ax6.set_ylim(0,np.pi)
		ax6.set_yticks([0,np.pi/2,np.pi])
		ax6.set_yticklabels(['0',r'$\frac{\pi}{2}$','$\pi$'],fontsize=12)
		ax6.spines['right'].set_visible(False)
		ax6.spines['top'].set_visible(False)

		ax2.get_xaxis().set_ticks([])
		ax2.get_yaxis().set_ticks([])
		ax2.set_frame_on(True)
		RightShoulder_r = plt.Circle((0,0),radius=0.05,Color='#4682b4')
		ax2.add_patch(RightShoulder_r)
		LeftShoulder_r = plt.Circle((-MedianPlane*2,0),radius=0.05,Color='#4682b4')
		ax2.add_patch(LeftShoulder_r)
		Torso_r = plt.Rectangle((-MedianPlane*2,-0.05),MedianPlane*2,0.1,Color='#4682b4')
		ax2.add_patch(Torso_r)
		Head_r = Ellipse((-MedianPlane,0),0.2,0.225,FaceColor='w',EdgeColor='#4682b4',Linewidth=3)
		ax2.add_patch(Head_r)
		JointCoordinates_r = \
			np.concatenate([[np.cumsum([0,L1*np.cos(A1_reverse[0,0]),L2*(0.146/(0.146+0.108))*np.cos(A1_reverse[0,0]+A2_reverse[0,0]),L2*(0.108/(0.146+0.108))*np.cos(A1_reverse[0,0]+A2_reverse[0,0])])],\
			[np.cumsum([0,L1*np.sin(A1_reverse[0,0]),L2*(0.146/(0.146+0.108))*np.sin(A1_reverse[0,0]+A2_reverse[0,0]),L2*(0.108/(0.146+0.108))*np.sin(A1_reverse[0,0]+A2_reverse[0,0])])]],\
			axis=0)
		Elbow_r = plt.Circle((L1*np.cos(A1_reverse[0,0]),L1*np.sin(A1_reverse[0,0])),radius=0.03,color='#4682b4')
		Endpoint_r = plt.Circle((L1*np.cos(A1_reverse[0,0])+L2*np.cos(A1_reverse[0,0]+A2_reverse[0,0]),\
							L1*np.sin(A1_reverse[0,0])+L2*np.sin(A1_reverse[0,0]+A2_reverse[0,0])),\
							radius = 0.02,color='#4682b4')
		Wrist_r = plt.Circle((L1*np.cos(A1_reverse[0,0])+L2*(0.146/(0.146+0.108))*np.cos(A1_reverse[0,0]+A2_reverse[0,0]),\
							L1*np.sin(A1_reverse[0,0])+L2*(0.146/(0.146+0.108))*np.sin(A1_reverse[0,0]+A2_reverse[0,0])),\
							radius = 0.03,color='#4682b4')
		StickFigure_r, = ax2.plot(JointCoordinates_r[0,:],JointCoordinates_r[1,:],'ko-',LineWidth=2,MarkerFaceColor='k')
		movement_r, = ax2.plot(L1*np.cos(A1_reverse[0,0])+L2*np.cos(A1_reverse[0,0]+A2_reverse[0,0]),\
		 					L1*np.sin(A1_reverse[0,0])+L2*np.sin(A1_reverse[0,0]+A2_reverse[0,0]),\
							color='0.60',linestyle=':')
		UpperArm_r = plt.Rectangle((0.02*np.sin(A1_reverse[0,0]),\
										-0.02*np.cos(A1_reverse[0,0])),\
										L1, 0.04,\
										angle=A1_reverse[0,0]*180/np.pi,color='#4682b4',animated=True)
		ax2.add_patch(UpperArm_r)
		Forearm_r = plt.Rectangle((L1*np.cos(A1_reverse[0,0]) \
						+ 0.02*np.sin(A1_reverse[0,0]+A2_reverse[0,0]),\
						L1*np.sin(A1_reverse[0,0]) -0.02*np.cos(A1_reverse[0,0]+A2_reverse[0,0])),\
						L2*(0.146/(0.146+0.108)), 0.04,\
						angle=(A1_reverse[0,0]+A2_reverse[0,0])*180/np.pi,color='#4682b4',animated=True)
		ax2.add_patch(Forearm_r)
		Hand_r = plt.Rectangle((L1*np.cos(A1_reverse[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.cos(A1_reverse[0,0]+A2_reverse[0,0]) \
									+ 0.02*np.sin(A1_reverse[0,0]+A2_reverse[0,0]),\
							L1*np.sin(A1_reverse[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.sin(A1_reverse[0,0]+A2_reverse[0,0]) \
									-0.02*np.cos(A1_reverse[0,0]+A2_reverse[0,0])),\
										L2*(0.108/(0.146+0.108)), 0.04,\
										angle=(A1_reverse[0,0]+A2_reverse[0,0])*180/np.pi,color='#4682b4',animated=True)
		ax2.add_patch(Hand_r)


		#Might need to add t,WeightedNormalizedMuscleVelocity_Forward,OrderNumber,etc
		NormalizedVmPlots_reverse = [ax4.plot(t.T,Vm_reverse[j].T) for j in OrderNumber]
		ax4.set_xlim(0,t_end)
		ax4.set_ylim(-12,12)
		ax4.set_xticks([0,t_end])
		ax4.set_xticklabels(['Start','Finish'])
		[k.set_color(OrderedColorsList[j]) for j,k in enumerate(ax4.lines)]
		Screen_r = plt.Rectangle((0,-13),t_end,26,color='w')
		ax4.add_patch(Screen_r)
		ax4.spines['right'].set_visible(False)
		ax4.spines['top'].set_visible(False)
		ax4.set_ylim([-1.1*bounds,1.1*bounds])

		max_x = np.concatenate(\
					[np.cumsum([0,\
						L1*np.cos(A1_reverse[0,i]),\
							L2*(0.146/(0.146+0.108))*np.cos(A1_reverse[0,i]+A2_reverse[0,i]),\
								L2*(0.108/(0.146+0.108))*np.cos(A1_reverse[0,i]+A2_reverse[0,i])]) \
									for i in range(len(t.T))],axis=0).max()
		max_y = np.concatenate(\
					[np.cumsum([0,\
						L1*np.sin(A1_reverse[0,i]),\
							L2*(0.146/(0.146+0.108))*np.sin(A1_reverse[0,i]+A2_reverse[0,i]),\
								L2*(0.108/(0.146+0.108))*np.sin(A1_reverse[0,i]+A2_reverse[0,i])]) \
									for i in range(len(t.T))],axis=0).max()
		min_x = np.concatenate(\
					[np.cumsum([0,\
						L1*np.cos(A1_reverse[0,i]),\
							L2*(0.146/(0.146+0.108))*np.cos(A1_reverse[0,i]+A2_reverse[0,i]),\
								L2*(0.108/(0.146+0.108))*np.cos(A1_reverse[0,i]+A2_reverse[0,i])]) \
									for i in range(len(t.T))],axis=0).min()
		if min_x>(-2*MedianPlane):
			min_x =	-2*MedianPlane

		ax2.set_xlim([min_x-0.1,max_x+0.1])
		ax2.set_ylim([-0.125,max_y+0.1])
		ax2.set_aspect('equal')

		ax1.set_xlim(ax2.get_xlim())
		ax1.set_ylim(ax2.get_ylim())
		ax1.set_aspect('equal')

		def animate(i):
			Angle1_f.set_xdata(t[:i])
			Angle1_f.set_ydata(A1_forward[0,:i])
			Angle2_f.set_xdata(t[:i])
			Angle2_f.set_ydata(A2_forward[0,:i])

			Angle1_r.set_xdata(t[:i])
			Angle1_r.set_ydata(A1_reverse[0,:i])
			Angle2_r.set_xdata(t[:i])
			Angle2_r.set_ydata(A2_reverse[0,:i])

			movement_f.set_xdata(list(map(lambda a1,a2: L1*np.cos(a1)+L2*np.cos(a1+a2),\
								A1_forward[0,:i],A2_forward[0,:i])))  # update the data
			movement_f.set_ydata(list(map(lambda a1,a2: L1*np.sin(a1)+L2*np.sin(a1+a2),\
								A1_forward[0,:i],A2_forward[0,:i])))
			Elbow_f.center = (L1*np.cos(A1_forward[0,i]),L1*np.sin(A1_forward[0,i]))
			Endpoint_f.center = (L1*np.cos(A1_forward[0,i])+L2*np.cos(A1_forward[0,i]+A2_forward[0,i]),\
								L1*np.sin(A1_forward[0,i])+L2*np.sin(A1_forward[0,i]+A2_forward[0,i]))
			Wrist_f.center = (L1*np.cos(A1_forward[0,i])+L2*(0.146/(0.146+0.108))*np.cos(A1_forward[0,i]+A2_forward[0,i]),\
								L1*np.sin(A1_forward[0,i])+L2*(0.146/(0.146+0.108))*np.sin(A1_forward[0,i]+A2_forward[0,i]))
			UpperArm_f._angle = A1_forward[0,i]*180/np.pi
			UpperArm_f.xy = (0.02*np.sin(A1_forward[0,i]),-0.02*np.cos(A1_forward[0,i]))
			Forearm_f._angle = (A1_forward[0,i]+A2_forward[0,i])*180/np.pi
			Forearm_f.xy = (L1*np.cos(A1_forward[0,i]) + 0.02*np.sin(A1_forward[0,i]+A2_forward[0,i]),\
							L1*np.sin(A1_forward[0,i]) -0.02*np.cos(A1_forward[0,i]+A2_forward[0,i]))
			Hand_f._angle = (A1_forward[0,i]+A2_forward[0,i])*180/np.pi
			Hand_f.xy = (L1*np.cos(A1_forward[0,i])\
							+ L2*(0.146/(0.146+0.108))*np.cos(A1_forward[0,i]+A2_forward[0,i]) \
								+ 0.02*np.sin(A1_forward[0,i]+A2_forward[0,i]),\
						L1*np.sin(A1_forward[0,i])\
							+ L2*(0.146/(0.146+0.108))*np.sin(A1_forward[0,i]+A2_forward[0,i]) \
								-0.02*np.cos(A1_forward[0,i]+A2_forward[0,i]))
			StickFigure_f.set_xdata(np.cumsum([0,L1*np.cos(A1_forward[0,i]),L2*(0.146/(0.146+0.108))*np.cos(A1_forward[0,i]+A2_forward[0,i]),L2*(0.108/(0.146+0.108))*np.cos(A1_forward[0,i]+A2_forward[0,i])]))
			StickFigure_f.set_ydata(np.cumsum([0,L1*np.sin(A1_forward[0,i]),L2*(0.146/(0.146+0.108))*np.sin(A1_forward[0,i]+A2_forward[0,i]),L2*(0.108/(0.146+0.108))*np.sin(A1_forward[0,i]+A2_forward[0,i])]))
			Screen_f.xy = (i/len(t.T),-13)
			Screen_f._width = t_end - i/len(t.T)

			movement_r.set_xdata(list(map(lambda a1,a2: L1*np.cos(a1)+L2*np.cos(a1+a2),\
								A1_reverse[0,:i],A2_reverse[0,:i])))  # update the data
			movement_r.set_ydata(list(map(lambda a1,a2: L1*np.sin(a1)+L2*np.sin(a1+a2),\
								A1_reverse[0,:i],A2_reverse[0,:i])))
			Elbow_r.center = (L1*np.cos(A1_reverse[0,i]),L1*np.sin(A1_reverse[0,i]))
			Endpoint_r.center = (L1*np.cos(A1_reverse[0,i])+L2*np.cos(A1_reverse[0,i]+A2_reverse[0,i]),\
								L1*np.sin(A1_reverse[0,i])+L2*np.sin(A1_reverse[0,i]+A2_reverse[0,i]))
			Wrist_r.center = (L1*np.cos(A1_reverse[0,i])+L2*(0.146/(0.146+0.108))*np.cos(A1_reverse[0,i]+A2_reverse[0,i]),\
								L1*np.sin(A1_reverse[0,i])+L2*(0.146/(0.146+0.108))*np.sin(A1_reverse[0,i]+A2_reverse[0,i]))
			UpperArm_r._angle = A1_reverse[0,i]*180/np.pi
			UpperArm_r.xy = (0.02*np.sin(A1_reverse[0,i]),-0.02*np.cos(A1_reverse[0,i]))
			Forearm_r._angle = (A1_reverse[0,i]+A2_reverse[0,i])*180/np.pi
			Forearm_r.xy = (L1*np.cos(A1_reverse[0,i]) + 0.02*np.sin(A1_reverse[0,i]+A2_reverse[0,i]),\
							L1*np.sin(A1_reverse[0,i]) -0.02*np.cos(A1_reverse[0,i]+A2_reverse[0,i]))
			Hand_r._angle = (A1_reverse[0,i]+A2_reverse[0,i])*180/np.pi
			Hand_r.xy = (L1*np.cos(A1_reverse[0,i])\
							+ L2*(0.146/(0.146+0.108))*np.cos(A1_reverse[0,i]+A2_reverse[0,i]) \
								+ 0.02*np.sin(A1_reverse[0,i]+A2_reverse[0,i]),\
						L1*np.sin(A1_reverse[0,i])\
							+ L2*(0.146/(0.146+0.108))*np.sin(A1_reverse[0,i]+A2_reverse[0,i]) \
								-0.02*np.cos(A1_reverse[0,i]+A2_reverse[0,i]))
			StickFigure_r.set_xdata(np.cumsum([0,L1*np.cos(A1_reverse[0,i]),L2*(0.146/(0.146+0.108))*np.cos(A1_reverse[0,i]+A2_reverse[0,i]),L2*(0.108/(0.146+0.108))*np.cos(A1_reverse[0,i]+A2_reverse[0,i])]))
			StickFigure_r.set_ydata(np.cumsum([0,L1*np.sin(A1_reverse[0,i]),L2*(0.146/(0.146+0.108))*np.sin(A1_reverse[0,i]+A2_reverse[0,i]),L2*(0.108/(0.146+0.108))*np.sin(A1_reverse[0,i]+A2_reverse[0,i])]))
			Screen_r.xy = (i/len(t.T),-13)
			Screen_r._width = t_end - i/len(t.T)
						# Arm.set_xdata(np.cumsum([0,L1*np.cos(A1_forward[0,i]),L2*(0.146/(0.146+0.108))*np.cos(A1_forward[0,i]+A2_forward[0,i]),L2*(0.108/(0.146+0.108))*np.cos(A1_forward[0,i]+A2_forward[0,i])]))
			# Arm.set_ydata(np.cumsum([0,L1*np.sin(A1_forward[0,i]),L2*(0.146/(0.146+0.108))*np.sin(A1_forward[0,i]+A2_forward[0,i]),L2*(0.108/(0.146+0.108))*np.sin(A1_forward[0,i]+A2_forward[0,i])]))
			return Elbow_f,Endpoint_f,Wrist_f,UpperArm_f,Forearm_f,Hand_f,StickFigure_f,movement_f,Screen_f,Elbow_r,Endpoint_r,Wrist_r,UpperArm_r,Forearm_r,Hand_r,StickFigure_r,movement_r,Screen_r,Angle1_f,Angle2_f,Angle1_r,Angle2_r,


		# Init only required for blitting to give a clean slate.
		def init():
		    # line.set_ydata(np.ma.array(x, mask=True))
			Elbow_f.center = (L1*np.cos(A1_forward[0,0]),L1*np.sin(A1_forward[0,0]))
			ax1.add_patch(Elbow_f)
			Endpoint_f.center = (L1*np.cos(A1_forward[0,0])+L2*np.cos(A1_forward[0,0]+A2_forward[0,0]),\
							L1*np.sin(A1_forward[0,0])+L2*np.sin(A1_forward[0,0]+A2_forward[0,0]))
			ax1.add_patch(Endpoint_f)
			Wrist_f.center = (L1*np.cos(A1_forward[0,0])+L2*(0.146/(0.146+0.108))*np.cos(A1_forward[0,0]+A2_forward[0,0]),\
							L1*np.sin(A1_forward[0,0])+L2*(0.146/(0.146+0.108))*np.sin(A1_forward[0,0]+A2_forward[0,0]))
			ax1.add_patch(Wrist_f)
			UpperArm_f = patches.Rectangle((0.02*np.sin(A1_forward[0,0]),\
		 								-0.02*np.cos(A1_forward[0,0])),\
	     								L1, 0.04,\
	     								angle=A1_forward[0,0]*180/np.pi,color='#4682b4')
			ax1.add_patch(UpperArm_f)
			Forearm_f = patches.Rectangle((L1*np.cos(A1_forward[0,0]) +\
			 			0.02*np.sin(A1_forward[0,0]+A2_forward[0,0]),\
						L1*np.sin(A1_forward[0,0]) -0.02*np.cos(A1_forward[0,0]+A2_forward[0,0])),\
						L2*(0.146/(0.146+0.108)), 0.04,\
						angle=(A1_forward[0,0]+A2_forward[0,0])*180/np.pi,color='#4682b4',animated=True)
			ax1.add_patch(Forearm_f)
			Hand_f = patches.Rectangle((L1*np.cos(A1_forward[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.cos(A1_forward[0,0]+A2_forward[0,0]) \
									+ 0.02*np.sin(A1_forward[0,0]+A2_forward[0,0]),\
							L1*np.sin(A1_forward[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.sin(A1_forward[0,0]+A2_forward[0,0]) \
									-0.02*np.cos(A1_forward[0,0]+A2_forward[0,0])),\
										L2*(0.108/(0.146+0.108)), 0.04,\
										angle=(A1_forward[0,0]+A2_forward[0,0])*180/np.pi,color='#4682b4',animated=True)
			ax1.add_patch(Hand_f)
			# Arm.set_xdata(np.ma.array(JointCoordinates[0,:], mask=True))
			# Arm.set_ydata(np.ma.array(JointCoordinates[1,:], mask=True))
			StickFigure_f.set_xdata(np.ma.array(JointCoordinates_f[0,:], mask=True))
			StickFigure_f.set_ydata(np.ma.array(JointCoordinates_f[1,:], mask=True))
			# movment.set_xdata()
			Screen_f = patches.Rectangle((0,-13),t_end,26,color='w',animated=True)
			ax3.add_patch(Screen_f)

			# line.set_ydata(np.ma.array(x, mask=True))
			Elbow_r.center = (L1*np.cos(A1_reverse[0,0]),L1*np.sin(A1_reverse[0,0]))
			ax2.add_patch(Elbow_r)
			Endpoint_r.center = (L1*np.cos(A1_reverse[0,0])+L2*np.cos(A1_reverse[0,0]+A2_reverse[0,0]),\
							L1*np.sin(A1_reverse[0,0])+L2*np.sin(A1_reverse[0,0]+A2_reverse[0,0]))
			ax2.add_patch(Endpoint_r)
			Wrist_r.center = (L1*np.cos(A1_reverse[0,0])+L2*(0.146/(0.146+0.108))*np.cos(A1_reverse[0,0]+A2_reverse[0,0]),\
							L1*np.sin(A1_reverse[0,0])+L2*(0.146/(0.146+0.108))*np.sin(A1_reverse[0,0]+A2_reverse[0,0]))
			ax2.add_patch(Wrist_r)
			UpperArm_r = patches.Rectangle((0.02*np.sin(A1_reverse[0,0]),\
		 								-0.02*np.cos(A1_reverse[0,0])),\
	     								L1, 0.04,\
	     								angle=A1_reverse[0,0]*180/np.pi,color='#4682b4')
			ax2.add_patch(UpperArm_r)
			Forearm_r = patches.Rectangle((L1*np.cos(A1_reverse[0,0]) +\
			 			0.02*np.sin(A1_reverse[0,0]+A2_reverse[0,0]),\
						L1*np.sin(A1_reverse[0,0]) -0.02*np.cos(A1_reverse[0,0]+A2_reverse[0,0])),\
						L2*(0.146/(0.146+0.108)), 0.04,\
						angle=(A1_reverse[0,0]+A2_reverse[0,0])*180/np.pi,color='#4682b4',animated=True)
			ax2.add_patch(Forearm_r)
			Hand_r = patches.Rectangle((L1*np.cos(A1_reverse[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.cos(A1_reverse[0,0]+A2_reverse[0,0]) \
									+ 0.02*np.sin(A1_reverse[0,0]+A2_reverse[0,0]),\
							L1*np.sin(A1_reverse[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.sin(A1_reverse[0,0]+A2_reverse[0,0]) \
									-0.02*np.cos(A1_reverse[0,0]+A2_reverse[0,0])),\
										L2*(0.108/(0.146+0.108)), 0.04,\
										angle=(A1_reverse[0,0]+A2_reverse[0,0])*180/np.pi,color='#4682b4',animated=True)
			ax2.add_patch(Hand_r)
			# Arm.set_xdata(np.ma.array(JointCoordinates[0,:], mask=True))
			# Arm.set_ydata(np.ma.array(JointCoordinates[1,:], mask=True))
			StickFigure_r.set_xdata(np.ma.array(JointCoordinates_r[0,:], mask=True))
			StickFigure_r.set_ydata(np.ma.array(JointCoordinates_r[1,:], mask=True))
			# movment.set_xdata()
			Screen_r = patches.Rectangle((0,-13),t_end,26,color='w',animated=True)
			ax4.add_patch(Screen_r)
			Elbow_f.set_visible(True)
			Endpoint_f.set_visible(True)
			Wrist_f.set_visible(True)
			UpperArm_f.set_visible(False)
			Forearm_f.set_visible(False)
			Hand_f.set_visible(False)
			StickFigure_f.set_visible(True)
			Screen_f.set_visible(True)
			Elbow_r.set_visible(True)
			Endpoint_r.set_visible(True)
			Wrist_r.set_visible(True)
			UpperArm_r.set_visible(False)
			Forearm_r.set_visible(False)
			Hand_r.set_visible(False)
			StickFigure_r.set_visible(True)
			Screen_r.set_visible(False)
			return Elbow_f,Endpoint_f,Wrist_f,UpperArm_f,Forearm_f,Hand_f,StickFigure_f,Screen_f,Elbow_r,Endpoint_r,Wrist_r,UpperArm_r,Forearm_r,Hand_r,StickFigure_r,Screen_r,

		ani = animation.FuncAnimation(fig, animate, np.arange(1, len(A1_forward[0,:]),10), init_func=init,interval=25, blit=True)
		if save_as_gif:
			ani.save('test.mp4',fps=50, dpi=100)
		plt.show()

def run_N_loops(NumberOfLoops):
	for LoopNumber in range(NumberOfLoops):
		XSplines = clamped_cubic_spline(0,EndTime,XInitial,XFinal,ẊInitial, \
												ẊFinal,XMin,XMax,Time,\
												angle = "Shoulder")
		Angle2Splines = clamped_cubic_spline(0,EndTime,Angle2Initial,Angle2Final,AngularVelocity2Initial, \
												AngularVelocity2Final,Angle2Bounds[0],Angle2Bounds[1],Time, \
												angle = "Elbow")
		Angle3Splines = clamped_cubic_spline(0,EndTime,Angle3Initial,Angle3Final,AngularVelocity3Initial, \
												AngularVelocity3Final,Angle3Bounds[0],Angle3Bounds[1],Time)
		if LoopNumber <= 8:
			print('-'*37 + 'End of Loop ' + str(LoopNumber+1) + '-'*37)
		else:
			print('-'*36 + 'End of Loop ' + str(LoopNumber+1) + '-'*37)
		pickle.dump([Angle1Splines, Angle2Splines, Angle3Splines], open('LoopNumber' + str(LoopNumber+1) + '.pkl','wb'),pickle.HIGHEST_PROTOCOL)

# run_N_loops(10)
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy as sp

import sys

	# Save Output Figures?

SaveOutputFigures = False

	# DescriptiveTitle should be something to identify the trial either by reaching location, (i.e., Left, Right, Center, or Sideways) or by what has changed in the most current iteration (e.g., CB_Ramsay, DELTa_Est, etc.). Spaces will be replaced by '_' symbols for the filename but kept for figure titles.

ValidResponse_1 = False
while ValidResponse_1 == False:
	ReachTypeNumber = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nPlease select reaching movement number:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n (1) - Side-to-side\n (2) - Straight (Center)\n (3) - 45° Left\n (4) - 45° Right\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nMovement Number: ")
	if ReachTypeNumber not in ['1','2','3','4']:
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInvalid Response! Please try again.')
		ValidResponse_1 = False
	else:
		ReachTypeNumber = int(ReachTypeNumber)-1
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		ValidResponse_1 = True
ReachType = ['Sideways','Center','Left','Right'][ReachTypeNumber]
DescriptiveTitle = ReachType + ' Reach'

ValidResponse_2 = False
while ValidResponse_2 == False:
	fix_starting_point_number = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nFix Starting Point?\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n (1) - True\n (2) - False\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nResponse: ")
	if fix_starting_point_number not in ['1','2']:
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInvalid Response! Please try again.')
		ValidResponse_2 = False
	else:
		fix_starting_point_number = int(fix_starting_point_number)-1
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		ValidResponse_2 = True
fix_starting_point = [True,False][fix_starting_point_number]

ValidResponse_3 = False
while ValidResponse_3 == False:
	fix_ending_point_number = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nFix Ending Point?\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n (1) - True\n (2) - False\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nResponse: ")
	if fix_ending_point_number not in ['1','2']:
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInvalid Response! Please try again.')
		ValidResponse_3 = False
	else:
		fix_ending_point_number = int(fix_ending_point_number)-1
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		ValidResponse_3 = True
fix_ending_point = [True,False][fix_ending_point_number]

x_initial_desired = 0.01 # cm
x_final_desired = 0.41 # cm
y_initial_desired = 0 # cm
y_final_desired = 0 # cm
boundary_sigma = 0.0025
allowable_error_rad = 30*(np.pi/180) #degrees
allowable_error_y = 0.05 # cm
NumberOfTrials =  100
dt = 1/1000
t = np.linspace(0,1,int(1/dt) + 1)
t_end = t[-1]

set_link_lengths()
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

###################################################################################################

	# Forward Direction

X_Forward = reaching_task(ReachType=ReachType,RandomXi=fix_starting_point,\
					RandomXf=fix_ending_point)
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

ScalingFactor = 200
# ScalingFactor is in units of N⋅s in order for its product to result in N⋅m⋅s⋅(Af#⋅lô/s) (Note: Af# and lô are unitless measures of corrected afferented number and normalized muscle length, respectively.

EccentricTorqueVariations_Shoulder_Forward = \
	ScalingFactor*eccentric_velocities(WeightedPotentialTorqueVariation_shoulder_Forward)
EccentricTorqueVariations_Elbow_Forward =\
 	ScalingFactor*eccentric_velocities(WeightedPotentialTorqueVariation_elbow_Forward)

# PotentialVariability_Forward = calculate_potential_variability(X_Forward,\
# 									EccentricTorqueVariations_Shoulder_Forward,\
# 										EccentricTorqueVariations_Elbow_Forward,\
# 											dt=dt,EOM = EOM,scheme = "Total")

	# Plot Afferent-Weighted Muscle Velocity (Forward)

fig1e = plt.figure()
[plt.plot(t.T,WeightedNormalizedMuscleVelocity_Forward[i].T) for i in OrderNumber]
ax1e = plt.gca()
ax1e.set_xlim(0,t_end*(1.3))
ax1e.set_ylim(-12,12)
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

EccentricCost_Forward = eccentric_cost(WeightedNormalizedMuscleVelocity_Forward,t_end=1,dt = dt,costtype='l1norm')
ConcentricCost_Forward = concentric_cost(WeightedNormalizedMuscleVelocity_Forward,t_end=1,dt = dt,costtype='l1norm')

# ###################################################################################################
#
# 	# Reverse Direction
#
# X_Reverse = reaching_task(Xi=Xf, Xf=Xi, dt=dt, t_end=t_end, EOM=EOM)
# return_MA_matrix()
#
# 	# Calculate resulting muscle velocities, lengths, etc.
#
# NormalizedMuscleVelocity_Reverse = calculate_muscle_velocities()
# WeightedPotentialTorqueVariation_shoulder_Reverse,\
# WeightedPotentialTorqueVariation_elbow_Reverse = \
# calculate_weighted_unscaled_potential_torque_variations_REACHING_TASK()
# WeightedNormalizedMuscleVelocity_Reverse = calculate_weighted_muscle_velocities_REACHING_TASK()
#
# 	# Calculate only the lengthening Torque Variation Contributions for both the shoulder and the elbow.
#
# EccentricTorqueVariations_Shoulder_Reverse = \
# 	ScalingFactor*eccentric_velocities(WeightedPotentialTorqueVariation_shoulder_Reverse)
# EccentricTorqueVariations_Elbow_Reverse =\
#  	ScalingFactor*eccentric_velocities(WeightedPotentialTorqueVariation_elbow_Reverse)
# PotentialVariability_Reverse = calculate_potential_variability(X_Reverse,\
# 									EccentricTorqueVariations_Shoulder_Reverse,\
# 										EccentricTorqueVariations_Elbow_Reverse,\
# 											dt=dt, EOM=EOM, scheme = "Total")
#
#
# 	# Test to make sure that the muscle velocities are negative, time-reversed versions of each other for the forward and backwards movement.
#
# assert np.array([np.array(abs(NormalizedMuscleVelocity_Forward[j,0,:].T-np.array(list(reversed(-NormalizedMuscleVelocity_Reverse[j,0,:]))).T)<1e-12).all() for j in range(n_muscles)]).all(), "The muscle velocities are not reversed and negative versions of each other."
#
# 	# Plot Afferent-Weighted Muscle Velocity (Reverse)
#
# fig2e = plt.figure()
# [plt.plot(t.T,WeightedNormalizedMuscleVelocity_Reverse[i].T) for i in OrderNumber]
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
# 	# Plot Normalized Muscle Velocity (Reverse)
#
# fig2a = plt.figure()
# [plt.plot(t.T,NormalizedMuscleVelocity_Reverse[i].T) for i in OrderNumber]
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
# 	# Plot Joint Angles (Reverse)
#
# fig2b = plt.figure()
# plt.plot(t.T,A1.T,'b')
# plt.plot(t.T,A2.T,'b:')
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
# 	# Plot Configuration Space (Reverse)
#
# fig2c = plt.figure()
# plt.plot(A1.T,A2.T,'b')
# ax2c = plt.gca()
# ax2c.set_title(DescriptiveTitle+'\nConfiguration Space')
# ax2c.set_ylabel("Shouler Angle (in radians)")
# ax2c.set_xlabel("Elbow Angle (in radians)")
#
# 	# Calculate Costs from Afferent-Weighted Muscle Velocities
#
# EccentricCost_Reverse = eccentric_cost(WeightedNormalizedMuscleVelocity_Reverse,t_end=t_end,dt = dt,costtype='l1norm')
# ConcentricCost_Reverse = concentric_cost(WeightedNormalizedMuscleVelocity_Reverse,t_end=t_end,dt = dt,costtype='l1norm')
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
# # if t_end!=1:
# # 	ax2d.set_xlabel('Time (s)')
# # else:
# # 	ax2d.set_xlabel('Normalized Time')
# # ax2d.set_ylabel('Muscle Lengths (in mm)')
# # [j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax2d.lines)]
# # ax2d.legend(OrderedMuscleList)
#
# ###################################################################################################
#
# 	# Plot bar graph comparing the two directions
#
# EccentricCost_Forward = PotentialVariability_Forward.mean()
# EccentricCost_Reverse = PotentialVariability_Reverse.mean()
# fig2 = plt.figure()
# plt.bar(np.arange(2),[EccentricCost_Forward,EccentricCost_Reverse])
# ax4 = plt.gca()
# ax4.set_xticklabels(('Forward','Reverse'))
# # ax4.set_ylim(0,20)
# # ax4.set_yticks([0,10,20])
# # ax4.set_yticklabels(['0','','20'])
# ax4.set_ylim(0,0.008)
# ax4.set_yticks([0,0.004,0.008])
# ax4.set_yticklabels(['0','','0.008'])
# ax4.set_title(DescriptiveTitle + '\nEccentric Cost for Forward and Reverse Movements')
# ax4.set_ylabel('Sum of Afferent-Weighted Muscle Lengthening')
# ax4.set_xticks([0,1])


###################################################################################################

animate_plots(True,ReachType,Weighted=False,save_as_gif=True)


##############################################################################


# fig1 = plt.figure()
# ax_rand = plt.gca()
# ax_rand.set_title(str(NumberOfTrials) + " Random Trajectories")
# plt.axes().set_aspect('equal', 'datalim')
#
# i = 0
# StartTime = time.time()
#
# while i < NumberOfTrials:
# 	if fix_starting_point == False:
# 		x_initial = x_initial_desired + np.random.normal(0,boundary_sigma) # cm
# 		y_initial = y_initial_desired + np.random.normal(0,boundary_sigma) # cm
# 	else:
# 		x_initial = x_initial_desired # cm
# 		y_initial = y_initial_desired # cm
#
# 	if fix_ending_point == False:
# 		x_final = x_final_desired + np.random.normal(0,boundary_sigma) # cm
# 		y_final = y_final_desired + np.random.normal(0,boundary_sigma) # cm
# 	else:
# 		x_final = x_final_desired # cm
# 		y_final = y_final_desired # cm
#
# 	initialerror = np.random.uniform(-np.tan(allowable_error_rad),np.tan(allowable_error_rad))
# 	finalerror = -np.sign(initialerror)*\
# 			abs(np.random.uniform(-np.tan(allowable_error_rad),np.tan(allowable_error_rad)))
# 	# finalerror = -np.random.normal(initialerror/2,abs(initialerror)/4)
# 	if initialerror>0:
# 		ymax = max([y_initial,y_final])+allowable_error_y
# 		ymin = 0
# 	else:
# 		ymax = 0
# 		ymin = min([y_initial,y_final])-allowable_error_y
#
# 	xmax = x_final
# 	xmin = x_initial
# 	x_rand = np.random.normal((x_initial+x_final)/2,abs(x_initial+x_final)/4)
# 	y_rand = np.random.normal((ymax+ymin)/2,abs(ymax+ymin)/4)
#
# 	A,B,C,D = spline_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,initialerror,finalerror)
# 	assert test_b_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,C,D,initialerror), "Initial slope does not match the expected value"
# 	assert test_endpoint_slope(B[1],C[1],D[1],x_rand,x_final,finalerror),"Problem with Endpoint Slope"
# 	assert test_for_discontinuity(A[0],B[0],C[0],D[0],x_initial,x_rand,A[1]), "Jump Discontinuity at t = %f!" %x_rand
# 	spline_structure = Spline(A,B,C,D,x_initial,x_rand,x_final)
# 	# statusbar(i,NumberOfTrials,Title="Reaching Task",StartTime=StartTime)
# 	if i == 0:
# 		if spline_structure.is_within_bounds(x_initial,x_final, ymin, ymax):
# 			Splines = spline_structure
# 			i+=1
# 			# statusbar(i,NumberOfTrials,Title="Reaching Task",StartTime=StartTime)
# 			# t = np.arange(0,1,0.001)
# 			# x = x_initial + (x_final-x_initial)*(10*t**3 - 15*t**4 + 6*t**5) # Returns an (1,N) array
# 			x = np.linspace(x_initial,x_final,1001)
# 			y = spline_structure.pp_func(x)
# 			line, = ax_rand.plot(x,y)
# 			c = line.get_color()
# 			ax_rand.scatter([x[0],x[-1]],[y[0],y[-1]],color=c,marker='o')
# 	elif i == 1:
# 		if spline_structure.is_within_bounds(x_initial,x_final, ymin, ymax):
# 			Splines = np.concatenate(([Splines], [spline_structure]), axis = 0)
# 			i+=1
# 			# statusbar(i,NumberOfTrials,Title="Reaching Task",StartTime=StartTime)
# 			# t = np.arange(0,1,0.001)
# 			# x = x_initial + (x_final-x_initial)*(10*t**3 - 15*t**4 + 6*t**5) # Returns an (1,N) array
# 			x = np.linspace(x_initial,x_final,1001)
# 			y = spline_structure.pp_func(x)
# 			line, = ax_rand.plot(x,y)
# 			c = line.get_color()
# 			ax_rand.scatter([x[0],x[-1]],[y[0],y[-1]],color=c,marker='o')
# 	else:
# 		if spline_structure.is_within_bounds(x_initial,x_final, ymin, ymax):
# 			Splines = np.concatenate((Splines, [spline_structure]), axis = 0)
# 			i+=1
# 			# t = np.arange(0,1,0.001)
# 			# x = x_initial + (x_final-x_initial)*(10*t**3 - 15*t**4 + 6*t**5) # Returns an (1,N) array
# 			x = np.linspace(x_initial,x_final,1001)
# 			y = spline_structure.pp_func(x)
# 			line, = ax_rand.plot(x,y)
# 			c = line.get_color()
# 			ax_rand.scatter([x[0],x[-1]],[y[0],y[-1]],color=c,marker='o')
# 			# if i<NumberOfTrials:
# 			# 	statusbar(i,NumberOfTrials,Title="Reaching Task",StartTime=StartTime)
#
#
#
# plt.plot([x_initial_desired,allowable_error_y/np.tan(allowable_error_rad)+x_initial_desired],[y_initial_desired, y_initial_desired + allowable_error_y],'k--',LineWidth=2)
# plt.plot([x_initial_desired,allowable_error_y/np.tan(allowable_error_rad)+x_initial_desired],[y_initial_desired, y_initial_desired - allowable_error_y],'k--',LineWidth=2)
#
# plt.plot([x_final_desired-allowable_error_y/np.tan(allowable_error_rad),x_final_desired],[y_final_desired+allowable_error_y, y_final_desired],'k--',LineWidth=2)
# plt.plot([x_final_desired-allowable_error_y/np.tan(allowable_error_rad),x_final_desired],[y_final_desired-allowable_error_y, y_final_desired],'k--',LineWidth=2)
#
# plt.plot([allowable_error_y/np.tan(allowable_error_rad)+x_initial_desired,x_final_desired-allowable_error_y/np.tan(allowable_error_rad)],[y_initial_desired+allowable_error_y,y_final_desired+allowable_error_y],'k--',LineWidth=2)
# plt.plot([allowable_error_y/np.tan(allowable_error_rad)+x_initial_desired,x_final_desired-allowable_error_y/np.tan(allowable_error_rad)],[y_initial_desired-allowable_error_y,y_final_desired-allowable_error_y],'k--',LineWidth=2)
#
# if fix_starting_point == False:
# 	x_circ_1 = np.linspace(-4*boundary_sigma + x_initial_desired + 0.00001,\
# 							4*boundary_sigma + x_initial_desired - 0.00001,\
# 								200)
# 	y_circ_1_top = np.array(list(map(lambda X: y_initial_desired+np.sqrt((4*boundary_sigma)**2-(X-x_initial_desired)**2),x_circ_1)))
# 	y_circ_1_bottom = np.array(list(map(lambda X: y_initial_desired-np.sqrt((4*boundary_sigma)**2-(X-x_initial_desired)**2),x_circ_1)))
#
# 	plt.plot(x_circ_1,y_circ_1_top,'k--',LineWidth=2)
# 	plt.plot(x_circ_1,y_circ_1_bottom,'k--',LineWidth=2)
#
# if fix_ending_point == False:
# 	x_circ_2 = np.linspace(-4*boundary_sigma + x_final_desired + 0.00001,\
# 							4*boundary_sigma + x_final_desired - 0.00001,\
# 								100)
# 	y_circ_2_top = np.array(list(map(lambda X: y_final_desired+np.sqrt((4*boundary_sigma)**2-(X-x_final_desired)**2),x_circ_2)))
# 	y_circ_2_bottom = np.array(list(map(lambda X: y_final_desired-np.sqrt((4*boundary_sigma)**2-(X-x_final_desired)**2),x_circ_2)))
#
# 	plt.plot(x_circ_2,y_circ_2_top,'k--',LineWidth=2)
# 	plt.plot(x_circ_2,y_circ_2_bottom,'k--',LineWidth=2)
#
#
# r_initial = np.sqrt(x_initial**2+y_initial**2)
# r_final = np.sqrt(x_final**2+y_final**2)
#
# # def test_ode_1(x,t):
# # 	def dr_dt(t):
# # 		return((r_final-r_initial)*(30*t**2-60*t**3+30*t**4))
# # 	def f(x):
# # 		if x<=spline_structure.x_break:
# # 			a,b,c,d = spline_structure.a[0],spline_structure.b[0,0,0],spline_structure.c[0,0],spline_structure.d[0,0,0]
# # 			return(a + b*(x-spline_structure.x_initial) + c*(x-spline_structure.x_initial)**2 + d*(x-spline_structure.x_initial)**3)
# # 		else:
# # 			a,b,c,d = spline_structure.a[1],spline_structure.b[1,0,0],spline_structure.c[1,0],spline_structure.d[1,0,0]
# # 			return(a + b*(x-spline_structure.x_break) + c*(x-spline_structure.x_break)**2 + d*(x-spline_structure.x_break)**3)
# #
# # 	def df(x):
# # 		if x<=spline_structure.x_break:
# # 			a,b,c,d = spline_structure.a[0],spline_structure.b[0,0,0],spline_structure.c[0,0],spline_structure.d[0,0,0]
# # 			return(b+ 2*c*(x-spline_structure.x_initial) + 3*d*(x-spline_structure.x_initial)**2)
# # 		else:
# # 			a,b,c,d = spline_structure.a[1],spline_structure.b[1,0,0],spline_structure.c[1,0],spline_structure.d[1,0,0]
# # 			return(b+ 2*c*(x-spline_structure.x_break) + 3*d*(x-spline_structure.x_break)**2)
# # 	return(dr_dt(t)*np.sqrt(x**2 + f(x)**2)/(x+df(x)*f(x)))
# fig1, (ax0, ax1) = plt.subplots(nrows=2, figsize=(7, 9.6))
# testx_1,testy_1 = spline_structure.return_parameterized_X()
# dS_desired_1 = spline_structure.dS_dt(t)
# testdS_1 = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),\
# 				np.gradient(testx_1)/dt,np.gradient(testy_1)/dt)))
# if sum(abs(dS_desired_1-testdS_1))/len(testdS_1)<1e-4:
# 	color_1 = 'g'
# else:
# 	color_1 = 'r'
# ax0.set_title("odeint - spline_structure piecewise function v0")
# ax0.plot(testx_1,testy_1,c=color_1,lw = 7)
# ax0.plot(np.linspace(x_initial,x_final,1001),spline_structure.pp_func(np.linspace(x_initial,x_final,1001)),'k')
#
# ax1.set_title("odeint - spline_structure piecewise function v0\nRecovered xy-plot")
# ax1.plot(t,testdS_1,c=color_1,lw=7)
# ax1.plot(t,dS_desired_1,'k')
# ax1.set_xlabel("Average Error: " + str(sum(abs(dS_desired_1-testdS_1))/len(testx_1)))


##############################################################################


# testx = [x_initial+0.01]
# for i in range(1,len(t)):
# 	testx.append(testx[-1]+(spline_structure.dS_dt(t[i])/np.sqrt(1+spline_structure.pp_deriv(testx[-1])**2))*dt)
# testx = np.array(testx)
# testy = spline_structure.pp_func(testx)
# testdS = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),\
# 						np.gradient(testx)/dt,np.gradient(testy)/dt)))
#
# fig2, (ax2, ax3) = plt.subplots(nrows=2, figsize=(7, 9.6))
# ax2.plot(testdS-dS_desired_1)
# ax2.set_title("Testing simple forward Euler againsts desired dr_dt")
# ax2.set_xlabel("Average Error: " + str(sum(abs(dS_desired_1-testdS))/len(testx)))
#
# ax3.plot(testdS_1-dS_desired_1)
# ax3.set_title("Testing odeint againsts desired dr_dt")
# ax3.set_xlabel("Average Error: " + str(sum(abs(dS_desired_1-testdS_1))/len(testx_1)))
#
# def dx_dt(t,x):
# 	r_initial = np.sqrt(x_initial**2+y_initial**2)
# 	r_final = np.sqrt(x_final**2+y_final**2)
# 	dr_dt = (r_final-r_initial)*(30*t**2 - 60*t**3 + 30*t**4)
# 	df_dx = np.piecewise(x,[x <= spline_structure.x_break, x > spline_structure.x_break], \
# 		[lambda x: spline_structure.b[0,0,0] + 	\
# 			2*spline_structure.c[0,0]*(x-spline_structure.x_initial) + \
# 			 	3*spline_structure.d[0,0,0]*(x-spline_structure.x_initial)**2, \
# 		lambda x: spline_structure.b[1,0,0] + \
# 			2*spline_structure.c[1,0]*(x-spline_structure.x_break) + \
# 				3*spline_structure.d[1,0,0]*(x-spline_structure.x_break)**2])
# 	result = dr_dt/np.sqrt(1+(df_dx)**2)
# 	return(result)
# xs=[x_initial]
# for i in range(len(t)-1):
# 	xs.append(fourth_order_runge_kutta(dx_dt,t[i],xs[i],dt))
# ys = np.array(list(map(lambda x: spline_structure.pp_func(x),xs)))
#
# fig2, (ax2, ax3) = plt.subplots(nrows=2, figsize=(7, 9.6))
#
# testr_2 = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),xs,ys)))
#
# r_desired_2 = np.array(list(map(lambda t: r_initial + (r_final-r_initial)*(10*t**3 - 15*t**4 + 6*t**5),t)))
# dr_desired_2 = np.array(list(map(lambda t: (r_final-r_initial)*(30*t**2 - 60*t**3 + 30*t**4),t)))
#
# if sum(abs(dr_desired_2-np.gradient(testr_2)/dt))/len(testr_2)<1e-4:
# 	color_2 = 'g'
# else:
# 	color_2 = 'r'
#
# ax2.set_title("RK4 - spline_structure piecewise function v1")
# ax2.plot(xs,ys,'r',lw = 7)
# ax2.plot(np.linspace(x_initial,x_final,1001),spline_structure.pp_func(np.linspace(x_initial,x_final,1001)),)
#
# ax3.set_title("RK4 - spline_structure piecewise function v1\nRecovered xy-plot")
# ax3.plot(t,np.gradient(testr_2)/dt,c=color_2,lw=7)
# ax3.plot(t,dr_desired_2,'k')
# r_break = np.sqrt(x_rand**2+y_rand**2)
# t_break = t[sum(r_break>r_desired_2)]
# ax3.plot([t_break,t_break],[0,15*(r_final-r_initial)/8],'k--')
#
#
# To show that the RK4 algorithm works...

# test_ode_1 = lambda x,t: (2*(30*t**2-60*t**3+30*t**4))/np.sqrt(4*x**2-8*x+5);
# dt = 1/1000
# t = np.linspace(0,1,int(1/dt) + 1)
# testx_3 = sp.integrate.odeint(test_ode_1,0,t).flatten()
# testy_3 = np.array(list(map(lambda x: -x**2+2*x,testx_3)))
# testr_3 = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),testx_3,testy_3)))
# testdr_3 = np.array(list(map(lambda dx,dy: np.sqrt(dx**2+dy**2),np.gradient(testx_3)\
# 												/dt,np.gradient(testy_3)/dt)))
# # assert sum(abs((2*(30*t**2-60*t**3+30*t**4))-testdr_3))/len(testr_3)<1e-4, "Error in parameterizing path to dr/dt. Check ODE func."
# r_desired_3 = np.array(list(map(lambda t: (2)*(10*t**3 - 15*t**4 + 6*t**5),t)))
# dr_desired_3 = np.array(list(map(lambda t: 2*(30*t**2-60*t**3+30*t**4),t)))
# if sum(abs(dr_desired_3-testdr_3))/len(testr_3)<1e-4:
# 	color_3 = 'g'
# else:
# 	color_3 = 'r'
# fig3, (ax4, ax5) = plt.subplots(nrows=2, figsize=(7, 9.6))
# ax4.set_title('odeint Test with $f(x) = -x^{2}+2x$')
# ax4.plot(t,testdr_3,c=color_3,lw=7)
# ax4.plot(t,dr_desired_3,'k')
#
# ax5.set_title("odeint Test with $f(x) = -x^{2}+2x$\nRecovered xy-plot")
# ax5.plot(testx_3,testy_3,c=color_3,lw=7)
# ax5.plot(np.linspace(0,2,101),np.array(list(map(lambda x: -x**2 + 2*x,np.linspace(0,2,101)))),'k')

# # def test_ode_2(t,x):
# # 	def dr_dt(t):
# # 		return((r_final-r_initial)*(30*t**2-60*t**3+30*t**4))
# # 	def f(x):
# # 		if x<=spline_structure.x_break:
# # 			a,b,c,d = spline_structure.a[0],spline_structure.b[0,0,0],spline_structure.c[0,0],spline_structure.d[0,0,0]
# # 			return(a + b*x + c*x**2 + d*x**3)
# # 		else:
# # 			a,b,c,d = spline_structure.a[1],spline_structure.b[1,0,0],spline_structure.c[1,0],spline_structure.d[1,0,0]
# # 			return(a + b*(x-spline_structure.x_break) + c*(x-spline_structure.x_break)**2 + d*(x-spline_structure.x_break)**3)
# # 	def df(x):
# # 		if x<=spline_structure.x_break:
# # 			a,b,c,d = spline_structure.a[0],spline_structure.b[0,0,0],spline_structure.c[0,0],spline_structure.d[0,0,0]
# # 			return(b+ 2*c*x + 3*d*x**2)
# # 		else:
# # 			a,b,c,d = spline_structure.a[1],spline_structure.b[1,0,0],spline_structure.c[1,0],spline_structure.d[1,0,0]
# # 			return(b+ 2*c*(x-spline_structure.x_break) + 3*d*(x-spline_structure.x_break)**2)
# # 	return(dr_dt(t)*np.sqrt(x**2 + f(x)**2)/(x+df(x)))
#
# def test_ode_2(t,x):
# 	def dr_dt(t):
# 		return((r_final-r_initial)*(30*t**2-60*t**3+30*t**4))
# 	def f(x):
# 		if x<=spline_structure.x_break:
# 			a,b,c,d = spline_structure.a[0],spline_structure.b[0,0,0],spline_structure.c[0,0],spline_structure.d[0,0,0]
# 			return(a + b*(x-spline_structure.x_initial) + c*(x-spline_structure.x_initial)**2 + d*(x-spline_structure.x_initial)**3)
# 		else:
# 			a,b,c,d = spline_structure.a[1],spline_structure.b[1,0,0],spline_structure.c[1,0],spline_structure.d[1,0,0]
# 			return(a + b*(x-spline_structure.x_break) + c*(x-spline_structure.x_break)**2 + d*(x-spline_structure.x_break)**3)
#
# 	def df(x):
# 		if x<=spline_structure.x_break:
# 			a,b,c,d = spline_structure.a[0],spline_structure.b[0,0,0],spline_structure.c[0,0],spline_structure.d[0,0,0]
# 			return(b+ 2*c*(x-spline_structure.x_initial) + 3*d*(x-spline_structure.x_initial)**2)
# 		else:
# 			a,b,c,d = spline_structure.a[1],spline_structure.b[1,0,0],spline_structure.c[1,0],spline_structure.d[1,0,0]
# 			return(b+ 2*c*(x-spline_structure.x_break) + 3*d*(x-spline_structure.x_break)**2)
# 	return(dr_dt(t)*np.sqrt(x**2 + f(x)**2)/(x+df(x)*f(x)))
#
# testx_4=[x_initial]
# for i in range(len(t)-1):
# 	testx_4.append(fourth_order_runge_kutta(test_ode_2,t[i],testx_4[i],dt))
# testy_4 = spline_structure.pp_func(np.array(testx_4))
# testr_4 = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),testx_4,testy_4)))
# dr_desired_4 = np.array(list(map(lambda t: (r_final-r_initial)*(30*t**2-60*t**3+30*t**4),t)))
# if sum(abs(dr_desired_4-np.gradient(testr_4)/dt))/len(testr_4)<1e-4:
# 	color_4 = 'g'
# else:
# 	color_4 = 'r'
# fig4, (ax6, ax7) = plt.subplots(nrows=2, figsize=(7, 9.6))
#
# ax6.set_title("RK4 - spline_structure piecewise function v2")
# ax6.plot(t,np.gradient(testr_4)/dt,c=color_4,lw=7)
# ax6.plot(t,np.array(list(map(lambda t: (r_final-r_initial)*(30*t**2-60*t**3+30*t**4),t))))
# ax6.plot([t_break,t_break],[0,15*(r_final-r_initial)/8],'k--')
#
# ax7.set_title("RK4 - spline_structure piecewise function v2\nRecovered xy-plot")
# ax7.plot(testx_4,testy_4,c=color_4,lw = 7)
# ax7.plot(np.linspace(x_initial,x_final,1001),spline_structure.pp_func(np.linspace(x_initial,x_final,1001)),)
#
#
# def test_ode_3(t,x):
# 	def dr_dt(t):
# 		return((2*np.sqrt(2))*(30*t**2-60*t**3+30*t**4))
# 	def g(x):
# 		return(x**2 + ((x-1)**3+1)**2)
# 	def dg(x):
# 		return(2*x + 2*((x-1)**3+1)*(3*(x-1)**2))
# 	return(dr_dt(t)*2*np.sqrt(g(x))/dg(x))
#
# testx_5=[0.0000001]
#
# for i in range(len(t)-1):
# 	testx_5.append(fourth_order_runge_kutta(test_ode_3,t[i],testx_5[i],dt))
# testy_5 = np.array(list(map(lambda x:(x-1)**3+1,testx_5)))
# testr_5 = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),testx_5,testy_5)))
# dr_desired_5 = np.array(list(map(lambda t: 2*np.sqrt(2)*(30*t**2-60*t**3+30*t**4),t)))
# if sum(abs(dr_desired_5-np.gradient(testr_5)/dt))/len(testr_5)<1e-4:
# 	color_5 = 'g'
# else:
# 	color_5 = 'r'
#
# fig5, (ax8, ax9) = plt.subplots(nrows=2, figsize=(7, 9.6))
# ax8.set_title('RK4 Test with $f(x) = (x-1)^{3}+1$')
#
# ax8.plot(t,np.gradient(testr_5)/dt,c=color_5,lw=7)
# ax8.plot(t,dr_desired_5,'k')
#
# ax9.set_title("RK4 Test with $f(x) = (x-1)^{3}+1$\nRecovered xy-plot")
# ax9.plot(testx_5,testy_5,c=color_5,lw=7)
# ax9.plot(np.linspace(0,2,101),np.array(list(map(lambda x: (x-1)**3 + 1,np.linspace(0,2,101)))),'k')
#
# # plt.figure()
# # plt.plot(testx,testy)
#
# # def test_ode_4(t,x):
# # 	a,b,c,d = spline_structure.a[0],spline_structure.b[0,0,0],spline_structure.c[0,0],spline_structure.d[0,0,0]
# # 	def dr_dt(t):
# # 		return((r_final-r_initial)*(30*t**2-60*t**3+30*t**4))
# # 	return(dr_dt(t)*2*np.sqrt(x**2 + spline_structure.pp_func(x))/(2*x + spline_structure.pp_deriv(x)))
# #
# # testx=[x_initial]
# #
# # for i in range(len(t)-1):
# # 	testx.append(fourth_order_runge_kutta(test_ode_4,t[i],testx[i],dt))
# # testy = np.array(list(map(lambda x:(x-1)**3+1,testx)))
# # testr = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),testx,testy)))
# # plt.figure()
# # plt.plot(t,np.gradient(testr)/dt,'r',lw=7)
# # plt.plot(t,np.array(list(map(lambda t: (r_final-r_initial)*(30*t**2-60*t**3+30*t**4),t))))
# # plt.plot([t_break,t_break],[0,15*(r_final-r_initial)/8],'k--')
#
# # plt.figure()
# # plt.plot(testx,testy)
#
# def test_ode_5(t,x):
# 	def dr_dt(t):
# 		return((2*np.sqrt(2))*(30*t**2-60*t**3+30*t**4))
# 	def g(x):
# 		if x <= 1:
# 			return(x**2 + (-x**2 + 2*x)**2)
# 		else:
# 			return(x**2 + (x**2 - 2*x + 2)**2)
# 	def dg(x):
# 		if x <= 1:
# 			return(2*x + 2*(-x**2 + 2*x)*(-2*x + 2))
# 		else:
# 			return(2*x + 2*(x**2 - 2*x + 2)*(2*x - 2))
# 	return(dr_dt(t)*2*np.sqrt(g(x))/dg(x))
#
# testx_6=[0.0000001]
#
# for i in range(len(t)-1):
# 	testx_6.append(fourth_order_runge_kutta(test_ode_5,t[i],testx_6[i],dt))
# testy_6 = np.array(list(map(lambda x:np.piecewise(x,[x<1,x>=1],[lambda x: -x**2 + 2*x, lambda x: x**2 - 2*x + 2]),testx_6)))
# testr_6 = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),testx_6,testy_6)))
# dr_desired_6 = np.array(list(map(lambda t: 2*np.sqrt(2)*(30*t**2-60*t**3+30*t**4),t)))
#
# fig6, (ax10, ax11) = plt.subplots(nrows=2, figsize=(7, 9.6))
# if sum(abs(dr_desired_6-np.gradient(testr_6)/dt))/len(testr_6)<1e-4:
# 	color_6 = 'g'
# else:
# 	color_6 = 'r'
# ax10.plot(t,np.gradient(testr_6)/dt,c=color_6,lw=7)
# ax10.plot(t,dr_desired_6,'k')
#
# ax10.set_title('RK4 Test with piecewise function\n$f(x) = -x^{2} + 2$ or $x^{2}-2x+2$;  $(x_{b} = 1)$')
#
# ax11.plot(testx_6,testy_6,c=color_6,lw=7)
# ax11.plot(np.linspace(0,2,1001),np.array(list(map(lambda x:np.piecewise(x,[x<1,x>=1],[lambda x: -x**2 + 2*x, lambda x: x**2 - 2*x + 2]),np.linspace(0,2,1001)))),'k')
# ax11.set_title("RK4 Test with piecewise function\nRecovered xy-plot")

plt.show()
