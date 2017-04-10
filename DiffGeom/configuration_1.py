def SymbolicPseudoRotation(n):
	import numpy as np 
	import sympy as sp 
	assert type(n) == int and n > 0, "n (number of DOF) must be a positive integer"
	theta = [sp.Symbol('theta'+str(i+1)) for i in range(n)]
	pseudoR = sp.Matrix(sp.Identity(n+1))
	for i in range(n):
		pseudoR[0,i+1] = -theta[i]
	return(pseudoR)

def testSymbolicPseudoRotation():
	import sympy as sp
	import numpy as np
	assert SymbolicPseudoRotation(3) \
			== sp.Matrix([	[1,	-sp.Symbol('angle'),	-sp.Symbol('theta2'),	-sp.Symbol('theta3')],\
							[0,	1,						0,						0					],\
							[0,	0,						1,						0					],\
							[0,	0,						0,						1					]]),\
			"Function failed to replicate PseudoRotationMatrix"
	assert np.shape(SymbolicPseudoRotation(10)) == (10+1,10+1), \
			"Function failed to build appropriately sized PseudoRotationMatrix"

# testSymbolicPseudoRotation()
	
def SymbolicRotationMatrix_SO3(axis,angle):
	import sympy as sp
	import itertools
	assert str(type(angle)) == "<class 'sympy.core.symbol.Symbol'>", "angle must be a sympy variable."
	assert axis in ['x','y','z'], "(rotational) axis must either be 'x', 'y', or 'z'."
	R = sp.Matrix(sp.Identity(3))
	if axis == 'x': i = [1,2]
	if axis == 'y': i = [0,2]
	if axis == 'z': i = [0,1]
	trig_indeces = list(itertools.product(i,i))
	R[trig_indeces[0]] = sp.cos(angle)
	R[trig_indeces[1]] = -sp.sin(angle)
	R[trig_indeces[2]] = sp.sin(angle)
	R[trig_indeces[3]] = sp.cos(angle)	
	return(R)

def testSymbolicRotationMatrix_SO3():
	import sympy as sp 
	import numpy as np
	angle = sp.Symbol('angle')
	test1 = sp.Matrix([[1,           0,            0],\
					[0, sp.cos(angle), -sp.sin(angle)],\
					[0, sp.sin(angle),  sp.cos(angle)]])
	test2 = sp.Matrix([[sp.cos(angle), 0, -sp.sin(angle)],\
					[          0, 1,            0],\
					[sp.sin(angle), 0,  sp.cos(angle)]])
	test3 = sp.Matrix([[sp.cos(angle), -sp.sin(angle), 0],\
					[sp.sin(angle),  sp.cos(angle), 0],\
					[          0,            0, 1]])
	assert test1 == SymbolicRotationMatrix_SO3('x',angle), "Failed for rotation about 'x' by 'angle'"
	assert test2 == SymbolicRotationMatrix_SO3('y',angle), "Failed for rotation about 'y' by 'angle'"
	assert test3 == SymbolicRotationMatrix_SO3('z',angle), "Failed for rotation about 'z' by 'angle'"
	assert np.shape(SymbolicRotationMatrix_SO3('x',angle)) == (3,r), "Failed to build an SO(3) Matrix"

# testSymbolicRotationMatrix_SO3()

def SymbolicMomentArmMatrix(m,n):
	import numpy as np
	import sympy as sp
	assert type(n) == int and n > 0, "n (number of DOF) must be a positive integer"
	assert type(m) == int and m > 0, "m (number of muscles) must be a positive integer"
	MA = sp.Matrix([[sp.Symbol('r_'+str(i+1)+','+str(j+1)) for j in range(m)] for i in range(n)])
	return(MA)

def SymbolicPositionVector(x,y,z):
	import sympy as sp
	for i in [x,y,z]: assert str(type(i)) == "<class 'sympy.core.symbol.Symbol'>", str(i) + " must be a sympy variable."
	x = sp.Matrix([[x],[y],[z]])
	return(x)

def SymbolicExcursionVector(m):
	import sympy as sp
	assert type(m) == int and m > 0, "m (number of muscles) must be a positive integer."
	s = [sp.Symbol('s_'+str(i+1)) for i in range(m)]
	return(s)

def SymbolicConfiguration_WithoutPosition(m,n):
	import sympy as sp 
	import numpy as np
	import scipy.linalg 
	pseudoR = SymbolicPseudoRotation(n)
	DiagPseudoR = []
	for i in range(m):
		DiagPseudoR = scipy.linalg.block_diag(DiagPseudoR,pseudoR)
	interimConfig = np.concatenate((DiagPseudoR,np.zeros(shape = (1,(n+1)*m))),axis=0)
	MA = SymbolicMomentArmMatrix(m,n)
	S = SymbolicExcursionVector(m)
	LastColumn = np.concatenate([np.concatenate(([[S[i]]], MA[:,i].T),axis=1) for i in range(m)],axis=1).T
	LastColumn = np.concatenate((LastColumn, [[1]]),axis=0)
	ConfigurationMatrix = np.concatenate((interimConfig,LastColumn),axis=1)
	return(ConfigurationMatrix)
	

	

