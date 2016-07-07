import numpy as np

def is_slope_positive(x_initial,x_transition,B,C,D,X):
	"""
	Takes in the initial x values for the two piecewise polynomials made by the clamped cubic splines with one
	break and evaluates if the slope is always positive. Returns TRUE if slope is positive for the first 250 ms
	of a 550 ms shot.
	"""
	derivative = np.piecewise(X,[X <= x_transition, X > x_transition], \
										[lambda X: B[0] + 2*C[0]*(X-x_initial) + 3*D[0]*(X-x_initial)**2, \
										lambda X: B[1] + 2*C[1]*(X-x_transition) + 3*D[1]*(X-x_transition)**2])
	result = min(derivative[:2501])>=0
	return(result)

def is_within_bounds(Spline,ymin,ymax):
	"""
	This takes in a 1D Spline array and tests to see if the values are within the allowed bounds [ymin, ymax].
	Returns TRUE if all values are within bounds.

	"""
	result = max(Spline)<=ymax and min(Spline)>=ymin
	return(result)

def piecewise_for_2_intervals(x_initial,x_transition,A,B,C,D,X):
	"""
	This is to generate a piecewise polynomial array for a cubic spline that has one break ONLY. A, B, C, and D 
	must each have two elements corresponding to the coefficients of each polynomial. x_initial is the initial 
	x value for the first polynomial and x_transition is the final x value for the first polynomial AND the 
	initial x value for the second polynomial (thus being the transition x value). X is a 1D array of the
	independent variable. Returns an array with len(result)=len(X).
	"""
	result = np.piecewise(X,[X <= x_transition, X > x_transition], \
										[lambda X: A[0] + B[0]*(X-x_initial) + C[0]*(X-x_initial)**2 + D[0]*(X-x_initial)**3, \
										lambda X: A[1] + B[1]*(X-x_transition) + C[1]*(X-x_transition)**2 + D[1]*(X-x_transition)**3])
	return(result)

