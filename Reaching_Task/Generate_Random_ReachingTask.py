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

	pp_func(X)
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
	def find_max_and_min(self,x_min,x_max,y_min,y_max):
		def find_extrema():
			extrema_1 = np.float(self.x_initial + (- 2*self.c[0,0] + (4*self.c[0,0]**2 - 12*self.b[0,0,0]*self.d[0,0,0])**.5)/(6*self.d[0,0,0]))
			extrema_2 = np.float(self.x_initial + (- 2*self.c[0,0] - (4*self.c[0,0]**2 - 12*self.b[0,0,0]*self.d[0,0,0])**.5)/(6*self.d[0,0,0]))
			extrema_3 = np.float(self.x_break + (- 2*self.c[1,0] + (4*self.c[1,0]**2 - 12*self.b[1,0,0]*self.d[1,0,0])**.5)/(6*self.d[1,0,0]))
			extrema_4 = np.float(self.x_break + (- 2*self.c[1,0] - (4*self.c[1,0]**2 - 12*self.b[1,0,0]*self.d[1,0,0])**.5)/(6*self.d[1,0,0]))
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
			if is_real(extrema_1) and is_in_appropriate_domain(extrema_1,x_min,x_max,1):
				if second_deriv_is_neg(extrema_1,1):
					maxima.append(np.float(self.pp_func(extrema_1)))
				elif second_deriv_is_pos(extrema_1,1):
					minima.append(np.float(self.pp_func(extrema_1)))
			if is_real(extrema_2) and is_in_appropriate_domain(extrema_2,x_min,x_max,1):
				if second_deriv_is_neg(extrema_2,1):
					maxima.append(np.float(self.pp_func(extrema_2)))
				elif second_deriv_is_pos(extrema_2,1):
					minima.append(np.float(self.pp_func(extrema_2)))
			if is_real(extrema_3) and is_in_appropriate_domain(extrema_3,x_min,x_max,2):
				if second_deriv_is_neg(extrema_3,2):
					maxima.append(np.float(self.pp_func(extrema_3)))
				elif second_deriv_is_pos(extrema_3,2):
					minima.append(np.float(self.pp_func(extrema_3)))
			if is_real(extrema_4) and is_in_appropriate_domain(extrema_4,x_min,x_max,2):
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
def clamped_cubic_spline(x_initial,x_final,y_initial,y_final,initial_slope,final_slope,ymin,ymax,X,**options):
	"""
	This will take in the initial and final values for both x and y, as well as the desired initial and final
	slopes and generate 1000 clamped cubic spline that produce y values that are within the bounds [ymin, ymax].
	Returns a list of arrays, each of len(X). Options allows for slope limitations on shoulder rotations such
	that the derivative of the spline is always positive to match observations (slope = "Shoulder").
	"""
	NumberOfTrials =  10000
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

x_initial_desired = 10 # cm
x_final_desired = 40 # cm
y_initial_desired = 0 # cm
y_final_desired = 0 # cm
fix_starting_point = True
fix_ending_point = True
boundary_sigma = 0.25
allowable_error_rad = 30*(np.pi/180) #degrees
allowable_error_y = 5 # cm
NumberOfTrials =  50


fig1 = plt.figure()
ax1 = plt.gca()
ax1.set_title(str(NumberOfTrials) + " Random Trajectories")
plt.axes().set_aspect('equal', 'datalim')

i = 0
StartTime = time.time()

while i < NumberOfTrials:
	if fix_starting_point == False:
		x_initial = x_initial_desired + np.random.normal(0,boundary_sigma) # cm
		y_initial = y_initial_desired + np.random.normal(0,boundary_sigma) # cm
	else:
		x_initial = x_initial_desired # cm
		y_initial = y_initial_desired # cm

	if fix_ending_point == False:
		x_final = x_final_desired + np.random.normal(0,boundary_sigma) # cm
		y_final = y_final_desired + np.random.normal(0,boundary_sigma) # cm
	else:
		x_final = x_final_desired # cm
		y_final = y_final_desired # cm

	initialerror = np.random.uniform(-np.tan(allowable_error_rad),np.tan(allowable_error_rad))
	finalerror = -np.sign(initialerror)*\
			abs(np.random.uniform(-np.tan(allowable_error_rad),np.tan(allowable_error_rad)))
	# finalerror = -np.random.normal(initialerror/2,abs(initialerror)/4)
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
	spline_structure = Spline(A,B,C,D,x_initial,x_rand,x_final)
	statusbar(i,NumberOfTrials,Title="Reaching Task",StartTime=StartTime)
	if i == 0:
		if spline_structure.is_within_bounds(x_initial,x_final, ymin, ymax):
			Splines = spline_structure
			i+=1
			statusbar(i,NumberOfTrials,Title="Reaching Task",StartTime=StartTime)
			t = np.arange(0,1,0.001)
			x = x_initial + (x_final-x_initial)*(10*t**3 - 15*t**4 + 6*t**5) # Returns an (1,N) array
			y = spline_structure.pp_func(x)
			line, = ax1.plot(x,y)
			c = line.get_color()
			ax1.scatter([x[0],x[-1]],[y[0],y[-1]],color=c,marker='o')
	elif i == 1:
		if spline_structure.is_within_bounds(x_initial,x_final, ymin, ymax):
			Splines = np.concatenate(([Splines], [spline_structure]), axis = 0)
			i+=1
			statusbar(i,NumberOfTrials,Title="Reaching Task",StartTime=StartTime)
			t = np.arange(0,1,0.001)
			x = x_initial + (x_final-x_initial)*(10*t**3 - 15*t**4 + 6*t**5) # Returns an (1,N) array
			y = spline_structure.pp_func(x)
			line, = ax1.plot(x,y)
			c = line.get_color()
			ax1.scatter([x[0],x[-1]],[y[0],y[-1]],color=c,marker='o')
	else:
		if spline_structure.is_within_bounds(x_initial,x_final, ymin, ymax):
			Splines = np.concatenate((Splines, [spline_structure]), axis = 0)
			i+=1
			t = np.arange(0,1,0.001)
			x = x_initial + (x_final-x_initial)*(10*t**3 - 15*t**4 + 6*t**5) # Returns an (1,N) array
			y = spline_structure.pp_func(x)
			line, = ax1.plot(x,y)
			c = line.get_color()
			ax1.scatter([x[0],x[-1]],[y[0],y[-1]],color=c,marker='o')
			if i<NumberOfTrials:
				statusbar(i,NumberOfTrials,Title="Reaching Task",StartTime=StartTime)



plt.plot([x_initial_desired,allowable_error_y/np.tan(allowable_error_rad)+x_initial_desired],[y_initial_desired, y_initial_desired + allowable_error_y],'k--',LineWidth=2)
plt.plot([x_initial_desired,allowable_error_y/np.tan(allowable_error_rad)+x_initial_desired],[y_initial_desired, y_initial_desired - allowable_error_y],'k--',LineWidth=2)

plt.plot([x_final_desired-allowable_error_y/np.tan(allowable_error_rad),x_final_desired],[y_final_desired+allowable_error_y, y_final_desired],'k--',LineWidth=2)
plt.plot([x_final_desired-allowable_error_y/np.tan(allowable_error_rad),x_final_desired],[y_final_desired-allowable_error_y, y_final_desired],'k--',LineWidth=2)

plt.plot([allowable_error_y/np.tan(allowable_error_rad)+x_initial_desired,x_final_desired-allowable_error_y/np.tan(allowable_error_rad)],[y_initial_desired+allowable_error_y,y_final_desired+allowable_error_y],'k--',LineWidth=2)
plt.plot([allowable_error_y/np.tan(allowable_error_rad)+x_initial_desired,x_final_desired-allowable_error_y/np.tan(allowable_error_rad)],[y_initial_desired-allowable_error_y,y_final_desired-allowable_error_y],'k--',LineWidth=2)

if fix_starting_point == False:
	x_circ_1 = np.linspace(-4*boundary_sigma + x_initial_desired,\
							4*boundary_sigma + x_initial_desired,\
								100)
	y_circ_1_top = np.array(list(map(lambda X: y_initial_desired+np.sqrt((4*boundary_sigma)**2-(X-x_initial_desired)**2),x_circ_1)))
	y_circ_1_bottom = np.array(list(map(lambda X: y_initial_desired-np.sqrt((4*boundary_sigma)**2-(X-x_initial_desired)**2),x_circ_1)))

	plt.plot(x_circ_1,y_circ_1_top,'k--',LineWidth=2)
	plt.plot(x_circ_1,y_circ_1_bottom,'k--',LineWidth=2)

if fix_ending_point == False:
	x_circ_2 = np.linspace(-4*boundary_sigma + x_final_desired,\
							4*boundary_sigma + x_final_desired,\
								100)
	y_circ_2_top = np.array(list(map(lambda X: y_final_desired+np.sqrt((4*boundary_sigma)**2-(X-x_final_desired)**2),x_circ_2)))
	y_circ_2_bottom = np.array(list(map(lambda X: y_final_desired-np.sqrt((4*boundary_sigma)**2-(X-x_final_desired)**2),x_circ_2)))

	plt.plot(x_circ_2,y_circ_2_top,'k--',LineWidth=2)
	plt.plot(x_circ_2,y_circ_2_bottom,'k--',LineWidth=2)

# plt.figure()

# def dx_dt(x,t):
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
# def ode45_step(f, x, t, dt, *args):
#     """
#     One step of 4th Order Runge-Kutta method
#     """
#     k = dt
#     k1 = k * f(t, x, *args)
#     k2 = k * f(t + 0.5*k, x + 0.5*k1, *args)
#     k3 = k * f(t + 0.5*k, x + 0.5*k2, *args)
#     k4 = k * f(t + dt, x + k3, *args)
#     return x + 1/6. * (k1 + k2 + k3 + k4)
#
# def ode45(f, t, x0, *args):
#     """
#     4th Order Runge-Kutta method
#     """
#     n = len(t)
#     x = np.zeros((n, np.size(x0)))
#     x[0] = x0
#     for i in range(n-1):
#         dt = t[i+1] - t[i]
#         x[i+1] = ode45_step(f, x[i], t[i], dt, *args)
#     return x
def fourth_order_runge_kutta(f,t,x,dt):
	def l_values(f,t,x,dt):
		l0 = f(t,x)*dt
		l1 = f(t+dt/2,x+l0/2)*dt
		l2 = f(t+dt/2,x+l1/2)*dt
		l3= f(t+dt,x+l2)*dt
		return([l0,l1,l2,l3])
	l = l_values(f,t,x,dt)
	next_x = x + (l[0] + 2*(l[1] + l[2]) + l[3])/6
	return(next_x)

dt = 1/10000
t = np.linspace(0,1,int(1/dt) + 1)
def dx_dt(t,x):
	r_initial = np.sqrt(x_initial**2+y_initial**2)
	r_final = np.sqrt(x_final**2+y_final**2)
	dr_dt = (r_final-r_initial)*(30*t**2 - 60*t**3 + 30*t**4)
	df_dx = np.piecewise(x,[x <= spline_structure.x_break, x > spline_structure.x_break], \
		[lambda x: spline_structure.b[0,0,0] + 	\
			2*spline_structure.c[0,0]*(x-spline_structure.x_initial) + \
			 	3*spline_structure.d[0,0,0]*(x-spline_structure.x_initial)**2, \
		lambda x: spline_structure.b[1,0,0] + \
			2*spline_structure.c[1,0]*(x-spline_structure.x_break) + \
				3*spline_structure.d[1,0,0]*(x-spline_structure.x_break)**2])
	result = dr_dt/np.sqrt(1+(df_dx)**2)
	return(result)
xs=[x_initial]
for i in range(len(t)-1):
	xs.append(fourth_order_runge_kutta(dx_dt,t[i],xs[i],dt))
ys = np.array(list(map(lambda x: spline_structure.pp_func(x),xs)))

# t = np.linspace(0,1,10000)
# xs = sp.integrate.odeint(dx_dt,x_initial,t).flatten()

fig2, (ax2, ax3) = plt.subplots(nrows=2, figsize=(7, 9.6))
ax2.set_title("RK4 - spline_structure piecewise function v1")
ax2.plot(xs,ys,'r',lw = 5)
ax2.plot(np.linspace(x_initial,x_final,1001),spline_structure.pp_func(np.linspace(x_initial,x_final,1001)),)

ax3.set_title("RK4 - spline_structure piecewise function v1\nRecovered xy-plot")
r_calculated = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),xs,ys)))
dr_calculated = np.gradient(r_calculated)/dt

r_initial = np.sqrt(x_initial**2+y_initial**2)
r_final = np.sqrt(x_final**2+y_final**2)
r_desired = np.array(list(map(lambda t: r_initial + (r_final-r_initial)*(10*t**3 - 15*t**4 + 6*t**5),t)))
dr_desired = np.array(list(map(lambda t: (r_final-r_initial)*(30*t**2 - 60*t**3 + 30*t**4),t)))

ax3.plot(t,dr_calculated,'r',lw=5)
ax3.plot(t,dr_desired)
r_break = np.sqrt(x_rand**2+y_rand**2)
t_break = t[sum(r_break>r_desired)]
ax3.plot([t_break,t_break],[0,15*(r_final-r_initial)/8],'k--')


# To show that the RK4 algorithm works...

test_ode_1 = lambda t,x: (2*(30*t**2-60*t**3+30*t**4))\
                    	*np.sqrt(x**2-4*x+5)/(2*x**2-6*x+5);
testx=[0]
for i in range(len(t)-1):
	testx.append(fourth_order_runge_kutta(test_ode_1,t[i],testx[i],dt))
testy = np.array(list(map(lambda x:-x**2+2*x,testx)))
testr = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),testx,testy)))
testr_desired = np.array(list(map(lambda t: (2)*(10*t**3 - 15*t**4 + 6*t**5),t)))

fig3, (ax4, ax5) = plt.subplots(nrows=2, figsize=(7, 9.6))
ax4.set_title('RK4 Test with $f(x) = -x^{2}+2x$')
ax4.plot(t,np.gradient(testr)/dt,'r',lw=5)
ax4.plot(t,np.array(list(map(lambda t: 2*(30*t**2-60*t**3+30*t**4),t))))

ax5.set_title("RK4 Test with $f(x) = -x^{2}+2x$\nRecovered xy-plot")
ax5.plot(testx,testy,'r',lw=5)
ax5.plot(np.linspace(0,2,101),np.array(list(map(lambda x: -x**2 + 2*x,np.linspace(0,2,101)))))

# def test_ode_2(t,x):
# 	def dr_dt(t):
# 		return((r_final-r_initial)*(30*t**2-60*t**3+30*t**4))
# 	def f(x):
# 		if x<=spline_structure.x_break:
# 			a,b,c,d = spline_structure.a[0],spline_structure.b[0,0,0],spline_structure.c[0,0],spline_structure.d[0,0,0]
# 			return(a + b*x + c*x**2 + d*x**3)
# 		else:
# 			a,b,c,d = spline_structure.a[1],spline_structure.b[1,0,0],spline_structure.c[1,0],spline_structure.d[1,0,0]
# 			return(a + b*(x-spline_structure.x_break) + c*(x-spline_structure.x_break)**2 + d*(x-spline_structure.x_break)**3)
# 	def df(x):
# 		if x<=spline_structure.x_break:
# 			a,b,c,d = spline_structure.a[0],spline_structure.b[0,0,0],spline_structure.c[0,0],spline_structure.d[0,0,0]
# 			return(b+ 2*c*x + 3*d*x**2)
# 		else:
# 			a,b,c,d = spline_structure.a[1],spline_structure.b[1,0,0],spline_structure.c[1,0],spline_structure.d[1,0,0]
# 			return(b+ 2*c*(x-spline_structure.x_break) + 3*d*(x-spline_structure.x_break)**2)
# 	return(dr_dt(t)*np.sqrt(x**2 + f(x)**2)/(x+df(x)))

def test_ode_2(t,x):
	def dr_dt(t):
		return((r_final-r_initial)*(30*t**2-60*t**3+30*t**4))
	def f(x):
		a,b,c,d = spline_structure.a[0],spline_structure.b[0,0,0],spline_structure.c[0,0],spline_structure.d[0,0,0]
		return(a + b*(x-spline_structure.x_initial) + c*(x-spline_structure.x_initial)**2 + d*(x-spline_structure.x_initial)**3)
	def df(x):
		a,b,c,d = spline_structure.a[0],spline_structure.b[0,0,0],spline_structure.c[0,0],spline_structure.d[0,0,0]
		return(b+ 2*c*(x-spline_structure.x_initial) + 3*d*(x-spline_structure.x_initial)**2)
	return(dr_dt(t)*np.sqrt(x**2 + f(x)**2)/(x+df(x)))

testx=[x_initial]
for i in range(len(t)-1):
	testx.append(fourth_order_runge_kutta(test_ode_2,t[i],testx[i],dt))
testy = spline_structure.pp_func(np.array(testx))
testr = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),testx,testy)))

fig4, (ax6, ax7) = plt.subplots(nrows=2, figsize=(7, 9.6))

ax6.set_title("RK4 - spline_structure piecewise function v2")
ax6.plot(t,np.gradient(testr)/dt,'r',lw=5)
ax6.plot(t,np.array(list(map(lambda t: (r_final-r_initial)*(30*t**2-60*t**3+30*t**4),t))))
ax6.plot([t_break,t_break],[0,15*(r_final-r_initial)/8],'k--')

ax7.set_title("RK4 - spline_structure piecewise function v2\nRecovered xy-plot")
ax7.plot(testx,testy,'r',lw = 5)
ax7.plot(np.linspace(x_initial,x_final,1001),spline_structure.pp_func(np.linspace(x_initial,x_final,1001)),)


def test_ode_3(t,x):
	def dr_dt(t):
		return((2*np.sqrt(2))*(30*t**2-60*t**3+30*t**4))
	def g(x):
		return(x**2 + ((x-1)**3+1)**2)
	def dg(x):
		return(2*x + 2*((x-1)**3+1)*(3*(x-1)**2))
	return(dr_dt(t)*2*np.sqrt(g(x))/dg(x))

testx=[0.0000001]

for i in range(len(t)-1):
	testx.append(fourth_order_runge_kutta(test_ode_3,t[i],testx[i],dt))
testy = np.array(list(map(lambda x:(x-1)**3+1,testx)))
testr = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),testx,testy)))

fig5, (ax8, ax9) = plt.subplots(nrows=2, figsize=(7, 9.6))
ax8.set_title('RK4 Test with $f(x) = (x-1)^{3}+1$')

ax8.plot(t,np.gradient(testr)/dt,'r',lw=5)
ax8.plot(t,np.array(list(map(lambda t: 2*np.sqrt(2)*(30*t**2-60*t**3+30*t**4),t))))

ax9.set_title("RK4 Test with $f(x) = (x-1)^{3}+1$\nRecovered xy-plot")
ax9.plot(testx,testy,'r',lw=5)
ax9.plot(np.linspace(0,2,101),np.array(list(map(lambda x: (x-1)**3 + 1,np.linspace(0,2,101)))))

# plt.figure()
# plt.plot(testx,testy)

# def test_ode_4(t,x):
# 	a,b,c,d = spline_structure.a[0],spline_structure.b[0,0,0],spline_structure.c[0,0],spline_structure.d[0,0,0]
# 	def dr_dt(t):
# 		return((r_final-r_initial)*(30*t**2-60*t**3+30*t**4))
# 	return(dr_dt(t)*2*np.sqrt(x**2 + spline_structure.pp_func(x))/(2*x + spline_structure.pp_deriv(x)))
#
# testx=[x_initial]
#
# for i in range(len(t)-1):
# 	testx.append(fourth_order_runge_kutta(test_ode_4,t[i],testx[i],dt))
# testy = np.array(list(map(lambda x:(x-1)**3+1,testx)))
# testr = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),testx,testy)))
# plt.figure()
# plt.plot(t,np.gradient(testr)/dt,'r',lw=5)
# plt.plot(t,np.array(list(map(lambda t: (r_final-r_initial)*(30*t**2-60*t**3+30*t**4),t))))
# plt.plot([t_break,t_break],[0,15*(r_final-r_initial)/8],'k--')

# plt.figure()
# plt.plot(testx,testy)

def test_ode_5(t,x):
	def dr_dt(t):
		return((2*np.sqrt(2))*(30*t**2-60*t**3+30*t**4))
	def g(x):
		if x <= 1:
			return(x**2 + (-x**2 + 2*x)**2)
		else:
			return(x**2 + (x**2 - 2*x + 2)**2)
	def dg(x):
		if x <= 1:
			return(2*x + 2*(-x**2 + 2*x)*(-2*x + 2))
		else:
			return(2*x + 2*(x**2 - 2*x + 2)*(2*x - 2))
	return(dr_dt(t)*2*np.sqrt(g(x))/dg(x))

testx=[0.0000001]

for i in range(len(t)-1):
	testx.append(fourth_order_runge_kutta(test_ode_5,t[i],testx[i],dt))
testy = np.array(list(map(lambda x:np.piecewise(x,[x<1,x>=1],[lambda x: -x**2 + 2*x, lambda x: x**2 - 2*x + 2]),testx)))
testr = np.array(list(map(lambda x,y:np.sqrt(x**2+y**2),testx,testy)))

fig6, (ax10, ax11) = plt.subplots(nrows=2, figsize=(7, 9.6))
ax10.plot(t,np.gradient(testr)/dt,'r',lw=5)
ax10.plot(t,np.array(list(map(lambda t: 2*np.sqrt(2)*(30*t**2-60*t**3+30*t**4),t))))
ax10.set_title('RK4 Test with piecewise function\n$f(x) = -x^{2} + 2$ or $x^{2}-2x+2$;  $(x_{b} = 1)$')

ax11.plot(testx,testy,'r',lw=5)
ax11.plot(np.linspace(0,2,1001),np.array(list(map(lambda x:np.piecewise(x,[x<1,x>=1],[lambda x: -x**2 + 2*x, lambda x: x**2 - 2*x + 2]),np.linspace(0,2,1001)))))
ax11.set_title("RK4 Test with piecewise function\nRecovered xy-plot")

plt.show()

# NumberOfTrials =  1000
#
# plt.figure()
# ax = plt.gca()
#
# StartTime = time.time()
#
# for i in range(NumberOfTrials):
# 	initialerror = np.random.normal(0,abs(allowable_error_rad)/2)
# 	finalerror = -np.random.normal(initialerror/2,abs(initialerror)/4)
# 	if initialerror>0:
# 		ymax = 10
# 		ymin = 0
# 	else:
# 		ymax = 0
# 		ymin = -10
# 	xmax = x_final
# 	xmin = x_initial
# 	x_rand = np.random.normal((x_initial+x_final)/2,abs(x_initial+x_final)/4)
# 	y_rand = np.random.normal((ymax+ymin)/2,abs(ymax+ymin)/4)
#
# 	A,B,C,D = spline_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,initialerror,finalerror)
#
# 	assert test_b_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,C,D,initialerror), "Initial slope does not match the expected value"
# 	assert test_endpoint_slope(B[1],C[1],D[1],x_rand,x_final,finalerror),"Problem with Endpoint Slope"
# 	assert test_for_discontinuity(A[0],B[0],C[0],D[0],x_initial,x_rand,A[1]), "Jump Discontinuity at t = %f!" %x_rand
#
# 	spline_structure = Spline(A,B,C,D,x_initial,x_rand,x_final)
#
# 	if i == 0:
# 		Splines = spline_structure
# 	elif i == 1:
# 		Splines = np.concatenate(([Splines], [spline_structure]), axis = 0)
# 	else:
# 		Splines = np.concatenate((Splines, [spline_structure]), axis = 0)
#
# 	t = np.arange(0,1,0.001)
# 	x = x_initial + (x_final-x_initial)*(10*t**3 - 15*t**4 + 6*t**5) # Returns an (1,N) array
# 	y = spline_structure.pp_func(x)
# 	ax.plot(x,y)
#
# 	statusbar(i,NumberOfTrials,Title = "Reaching Task",StartTime=StartTime)

# plt.show()
