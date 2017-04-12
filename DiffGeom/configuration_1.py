def SymbolicPseudoRotation(n,AngleSymbol='theta'):
	import numpy as np 
	import sympy as sp 
	assert type(n) == int and n > 0, "n (number of DOF) must be a positive integer"
	assert type(AngleSymbol) == str, "AngleSymbol must be a string."
	Angle = [sp.Symbol(AngleSymbol+str(i+1)) for i in range(n)]
	PseudoR = sp.Matrix(sp.Identity(n+1))
	for i in range(n):
		PseudoR[0,i+1] = -Angle[i]
	return(PseudoR)
def testSymbolicPseudoRotation():
	import sympy as sp
	import numpy as np
	assert SymbolicPseudoRotation(3) \
			== sp.Matrix([	[1,	-sp.Symbol('theta1'),	-sp.Symbol('theta2'),	-sp.Symbol('theta3')],\
							[0,	1,						0,						0					],\
							[0,	0,						1,						0					],\
							[0,	0,						0,						1					]]),\
			"Function failed to replicate PseudoRotationMatrix"
	assert np.shape(SymbolicPseudoRotation(10)) == (10+1,10+1), \
			"Function failed to build appropriately sized PseudoRotationMatrix"	
# testSymbolicPseudoRotation()
def SymbolicRotationMatrix_SO3(RotationalAxis,Angle):
	import sympy as sp
	import itertools
	assert str(type(Angle)) == "<class 'sympy.core.symbol.Symbol'>", "Angle must be a sympy variable."
	assert RotationalAxis in ['x','y','z'], "(rotational) RotationalAxis must either be 'x', 'y', or 'z'."
	R = sp.Matrix(sp.Identity(3))
	if RotationalAxis == 'x': i = [1,2]
	if RotationalAxis == 'y': i = [0,2]
	if RotationalAxis == 'z': i = [0,1]
	TrigonometricIndeces = list(itertools.product(i,i))
	R[TrigonometricIndeces[0]] = sp.cos(Angle)
	R[TrigonometricIndeces[1]] = -sp.sin(Angle)
	R[TrigonometricIndeces[2]] = sp.sin(Angle)
	R[TrigonometricIndeces[3]] = sp.cos(Angle)	
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
def SymbolicMomentArmMatrix(m,n,MomentArmSymbol='r'):
	import numpy as np
	import sympy as sp
	import itertools as it
	assert type(n) == int and n > 0, "n (number of DOF) must be a positive integer"
	assert type(m) == int and m > 0, "m (number of muscles) must be a positive integer"
	assert type(MomentArmSymbol) == str, "MomentArmSymbol must be a string"
	MA = sp.Matrix(np.array([sp.Symbol(MomentArmSymbol + '_' + str(el)) for el in it.product(range(1,n+1),range(1,m+1))]).reshape(n,m))
	return(MA)
def SymbolicPositionVector(StringX='x',StringY='y',StringZ='z'):
	import sympy as sp
	for i in [StringX,StringY,StringZ]: assert type(i) == str,  "Position variables must be a strings."
	Position = sp.Matrix([[sp.Symbol(StringX)],[sp.Symbol(StringY)],[sp.Symbol(StringZ)]])
	return(Position)
def SymbolicExcursionVector(m,ExcursionSymbol='s'):
	import sympy as sp
	assert type(m) == int and m > 0, "m (number of muscles) must be a positive integer."
	assert type(ExcursionSymbol) == str, "ExcursionSymbol must be a string."
	Excursion = [sp.Symbol(ExcursionSymbol+'_'+str(i+1)) for i in range(m)]
	return(Excursion)
def SymbolicConfiguration_WithoutPosition(m,n,AngleSymbol='theta',MomentArmSymbol='r',ExcursionSymbol='s'):
	import sympy as sp 
	import numpy as np
	import scipy.linalg 
	PseudoR = SymbolicPseudoRotation(n,AngleSymbol=AngleSymbol)
	DiagPseudoR = []
	for i in range(m):
		DiagPseudoR = scipy.linalg.block_diag(DiagPseudoR,PseudoR)
	InterimConfigurationMatrix = np.concatenate((DiagPseudoR,np.zeros(shape = (1,(n+1)*m))),axis=0)
	MA = SymbolicMomentArmMatrix(m,n,MomentArmSymbol=MomentArmSymbol)
	S = SymbolicExcursionVector(m,ExcursionSymbol=ExcursionSymbol)
	LastColumn = np.concatenate([np.concatenate(([[S[i]]], MA[:,i].T),axis=1) for i in range(m)],axis=1).T
	LastColumn = np.concatenate((LastColumn, [[1]]),axis=0)
	ConfigurationMatrix = sp.Matrix(np.concatenate((InterimConfigurationMatrix,LastColumn),axis=1))
	return(ConfigurationMatrix)
def ReturnSymbolicTensor(m,n,g):
	import sympy as sp
	import numpy as np
	g_tensor = sp.Matrix(np.concatenate([-g[0,1:n+1],g[range(0,(n+1)*m-1,n+1),-1].T,g[[i for i in range((n+1)*m) if i not in range(0,(n+1)*m-1,n+1)],-1].T],axis=1))
	return(g_tensor)
def SymbolicLeftLiftedAction_WithoutPosition(m,n,AngleSymbol='theta',MomentArmSymbol='r',ExcursionSymbol='s'):
	import sympy as sp
	import numpy as np
	if AngleSymbol == 'alpha': 
		Angle_g = 'theta'
	else:
		Angle_g = 'alpha'
	if MomentArmSymbol == 'ma':
		MA_g = 'r'
	else:
		MA_g = 'ma'
	if ExcursionSymbol == 'l':
		S_g = 's'
	else:
		S_g = 'l'
	g = SymbolicConfiguration_WithoutPosition(m,n,AngleSymbol=Angle_g,MomentArmSymbol=MA_g,ExcursionSymbol=S_g)
	g_tensor = ReturnSymbolicTensor(m,n,g)
	h = SymbolicConfiguration_WithoutPosition(m,n,AngleSymbol=AngleSymbol,MomentArmSymbol=MomentArmSymbol,ExcursionSymbol=ExcursionSymbol)
	hg = h*g
	hg_tensor = ReturnSymbolicTensor(m,n,hg)
	TgLh = hg_tensor.jacobian([g_tensor])
	return(TgLh)
def SymbolicRightLiftedAction_WithoutPosition(m,n,AngleSymbol='theta',MomentArmSymbol='r',ExcursionSymbol='s'):
	import sympy as sp
	import numpy as np
	if AngleSymbol == 'alpha': 
		Angle_g = 'theta'
	else:
		Angle_g = 'alpha'
	if MomentArmSymbol == 'ma':
		MA_g = 'r'
	else:
		MA_g = 'ma'
	if ExcursionSymbol == 'l':
		S_g = 's'
	else:
		S_g = 'l'
	g = SymbolicConfiguration_WithoutPosition(m,n,AngleSymbol=Angle_g,MomentArmSymbol=MA_g,ExcursionSymbol=S_g)
	g_tensor = ReturnSymbolicTensor(m,n,g)
	h = SymbolicConfiguration_WithoutPosition(m,n,AngleSymbol=AngleSymbol,MomentArmSymbol=MomentArmSymbol,ExcursionSymbol=ExcursionSymbol)
	gh = g*h
	gh_tensor = ReturnSymbolicTensor(m,n,gh)
	TgRh = gh_tensor.jacobian([g_tensor])
	return(TgRh)
def SymbolicInverseConfiguration_WithoutPosition(m,n,AngleSymbol='theta',MomentArmSymbol='r',ExcursionSymbol='s'):
	import sympy as sp
	import numpy as np
	g = SymbolicConfiguration_WithoutPosition(m,n,AngleSymbol=AngleSymbol,MomentArmSymbol=MomentArmSymbol,ExcursionSymbol=ExcursionSymbol)
	Inverse_g = g**-1
	return(Inverse_g)
def SymbolicTeLg(m,n,AngleSymbol='theta',MomentArmSymbol='r',ExcursionSymbol='s'):
	import sympy as sp
	import numpy as np
	import itertools as it
	if AngleSymbol == 'THETA': 
		Angle_h = 'theta'
	else:
		Angle_h = 'THETA'
	if MomentArmSymbol == 'R':
		MA_h = 'r'
	else:
		MA_h = 'R'
	if ExcursionSymbol == 'S':
		S_h = 's'
	else:
		S_h = 'S'
	TgLh = SymbolicLeftLiftedAction_WithoutPosition(m,n,AngleSymbol=Angle_h,MomentArmSymbol=MA_h,ExcursionSymbol=S_h)
	h = SymbolicConfiguration_WithoutPosition(m,n,AngleSymbol=Angle_h,MomentArmSymbol=MA_h,ExcursionSymbol=S_h)
	h_tensor = ReturnSymbolicTensor(m,n,h)
	Inverse_g = SymbolicInverseConfiguration_WithoutPosition(m,n,AngleSymbol=AngleSymbol,MomentArmSymbol=MomentArmSymbol,ExcursionSymbol=ExcursionSymbol)
	Inverse_g_tensor = ReturnSymbolicTensor(m,n,Inverse_g)
	TgLg_inv = TgLh.subs([(str(h_tensor[i]),str(Inverse_g_tensor[i])) for i in range(len(h_tensor))])
	TeLg = TgLg_inv**-1
	return(TeLg)
def SymbolicTeRg(m,n,AngleSymbol='theta',MomentArmSymbol='r',ExcursionSymbol='s'):
	import sympy as sp
	import numpy as np
	import itertools as it
	if AngleSymbol == 'THETA': 
		Angle_h = 'theta'
	else:
		Angle_h = 'THETA'
	if MomentArmSymbol == 'R':
		MA_h = 'r'
	else:
		MA_h = 'R'
	if ExcursionSymbol == 'S':
		S_h = 's'
	else:
		S_h = 'S'
	TgRh = SymbolicRightLiftedAction_WithoutPosition(m,n,AngleSymbol=Angle_h,MomentArmSymbol=MA_h,ExcursionSymbol=S_h)
	h = SymbolicConfiguration_WithoutPosition(m,n,AngleSymbol=Angle_h,MomentArmSymbol=MA_h,ExcursionSymbol=S_h)
	h_tensor = ReturnSymbolicTensor(m,n,h)
	Inverse_g = SymbolicInverseConfiguration_WithoutPosition(m,n,AngleSymbol=AngleSymbol,MomentArmSymbol=MomentArmSymbol,ExcursionSymbol=ExcursionSymbol)
	Inverse_g_tensor = ReturnSymbolicTensor(m,n,Inverse_g)
	TgRg_inv = TgRh.subs([(str(h_tensor[i]),str(Inverse_g_tensor[i])) for i in range(len(h_tensor))])
	TeRg = TgRg_inv**-1
	return(TeRg)
