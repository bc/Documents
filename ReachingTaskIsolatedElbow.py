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
def create_angle_lists():
	"""
	Creates global lists of angles, angular velocities, and angular accelerations.
	"""
	global A1,A2,Ȧ1,Ȧ2,Ä1,Ä2
	A1,A2=[],[]
	Ȧ1,Ȧ2=[],[]
	Ä1,Ä2=[],[]
def elbow_swing(t,Ai,Af,t_end):
	global A1,A2,Ȧ1,Ȧ2,Ä1,Ä2
	A1 = Ai[0] + (Af[0]-Ai[0])*(10*t**3 - 15*t**4 + 6*t**5)
	A2 = Ai[1] + (Af[1]-Ai[1])*(10*t**3 - 15*t**4 + 6*t**5)
	Ȧ1 = (Af[0]-Ai[0])*(30*t**2 - 60*t**3 + 30*t**4)/t_end
	Ȧ2 = (Af[1]-Ai[1])*(30*t**2 - 60*t**3 + 30*t**4)/t_end
	Ä1 = (Af[0]-Ai[0])*(60*t - 180*t**2 + 120*t**3)/(t_end**2)
	Ä2 = (Af[1]-Ai[1])*(60*t - 180*t**2 + 120*t**3)/(t_end**2)
def reaching_task_ISOLATED_ELBOW_TASK(dt = 0.001, Ai = [3.14159/2,3.14159/2], Af = [3.14159/2,3.14159/4],t_end=1):
    import numpy as np
    set_link_lengths()
    create_angle_lists()
    t = np.arange(0,1+dt,dt)
    elbow_swing(t,Ai,Af,t_end)
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

	DELTa SFE MA is listed as 33.02 mm in Pigeon and estimated as 19 mm. Using Pigeon Coefficients convention, Kuechle (1997) has the DELTp MA for [-140,90] as Pigeon_coeff_conversion([ 1.34293189,  0.20316226, -0.02339031,  0.27807828,  0.,  0.]).

	DELTp SFE MA is listed as -78.74 mm in Pigeon and estimated as -8 mm. Using Pigeon Coefficients convention, Kuechle (1997) has the DELTp MA for [-140,90] as Pigeon_coeff_conversion([ 2.28547177,  0.39721238, -0.33900829, -0.36146546,  0.,  0.]).

	PC (Clavicle attachment of Pectoralis) SFE MA is listed as 50.80 mm in Pigeon.
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
			'MA' : Pigeon_coeff_conversion([ 1.34293189,  0.20316226, -0.02339031,  0.27807828,  0.,  0.]),\
			'src' : 'Pigeon', 'eq' : None, 'threshold' : np.pi/2, \
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
			'MA' : Pigeon_coeff_conversion([ -2.28547177,  -0.39721238, 0.33900829, 0.36146546,  0.,  0.]), \
			'src' : 'Pigeon', 'eq' : None, 'threshold' : np.pi/2, \
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
def calculate_weighted_muscle_velocities_ISOLATED_ROTATION(A1,A2,Ȧ1,Ȧ2,dof='Elbow'):
    import numpy as np
    import sympy as sp
    global q1,q2,q_PS,AllCoefficients
    R = return_MA_matrix(A1,A2)
    if dof =='Elbow':
        MomentArmList = [1]
    elif dof == 'Shoulder':
        MomentArmList = [0]
    else:
        MomentArmList = [0,1]
    #dR = np.concatenate([sp.diff(RMatrix_Transpose[:,0],q1).subs(q1,A1),sp.diff(RMatrix_Transpose[:,1],q2).subs([(q2,A2),(q_PS,np.pi/2)])],axis=1)
    MuscleVelocity = -R.T*(np.matrix([Ȧ1,Ȧ2]).T)#-dR*(np.matrix(np.multiply([A1,A2],[Ȧ1,Ȧ2])).T)#
    MuscleList = AllCoefficients.keys()
    # WeightedNormalizedMuscleVelocity = np.array([[float(MuscleVelocity[i,0])/AllCoefficients[list(MuscleList)[i]]['Optimal Muscle Length'] for i in range(len(MuscleVelocity))]])

    WeightedNormalizedMuscleVelocity = np.array([[float(MuscleVelocity[i,0])*abs(sum([R.T[i,col] for col in MomentArmList])/1000)*AllCoefficients[list(MuscleList)[i]]['Corrected No']/AllCoefficients[list(MuscleList)[i]]['Optimal Muscle Length'] for i in range(len(AllCoefficients))]])
    return(WeightedNormalizedMuscleVelocity)
def calculate_muscle_velocities(A1,A2,Ȧ1,Ȧ2):
    import numpy as np
    import sympy as sp
    global q1,q2,q_PS,AllCoefficients
    R = return_MA_matrix(A1,A2)
    #dR = np.concatenate([sp.diff(RMatrix_Transpose[:,0],q1).subs(q1,A1),sp.diff(RMatrix_Transpose[:,1],q2).subs([(q2,A2),(q_PS,np.pi/2)])],axis=1)
    MuscleVelocity = -R.T*(np.matrix([Ȧ1,Ȧ2]).T)#-dR*(np.matrix(np.multiply([A1,A2],[Ȧ1,Ȧ2])).T)#
    MuscleList = list(AllCoefficients.keys())
    # WeightedNormalizedMuscleVelocity = np.array([[float(MuscleVelocity[i,0])/AllCoefficients[list(MuscleList)[i]]['Optimal Muscle Length'] for i in range(len(MuscleVelocity))]])

    NormalizedMuscleVelocity = np.array([[float(MuscleVelocity[i,0])/AllCoefficients[MuscleList[i]]['Optimal Muscle Length'] for i in range(len(AllCoefficients))]])
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

import numpy as np
import matplotlib.pyplot as plt
import math
import time

dof = 'Shoulder'
t_end = 1
dt = 0.001
t = np.arange(0,1+dt,dt)*t_end
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

# Forward Direction

if dof == 'Elbow':
    # Elbow Rotation
    Ai = [np.pi/2,np.pi/2]
    Af = [np.pi/2,np.pi/4]
elif dof == 'Shoulder':
    # Shoulder Rotation
    Ai = [np.pi/2-np.pi/8,0]
    Af = [np.pi/2+np.pi/8,0]

reaching_task_ISOLATED_ELBOW_TASK(dt = dt, Ai=Ai, Af=Af, t_end=t_end)
WeightedNormalizedMuscleVelocity_Forward = calculate_weighted_muscle_velocities_ISOLATED_ROTATION(A1[0],A2[0],Ȧ1[0],Ȧ2[0],dof=dof)
NormalizedMuscleVelocity_Forward = calculate_muscle_velocities(A1[0],A2[0],Ȧ1[0],Ȧ2[0])
StartTime = time.time()
statusbar(0,len(A1),Title = 'Forward Vm',StartTime=StartTime)
for i in range(1,len(A1)):
    statusbar(i,len(A1),Title = 'Forward Vm',StartTime=StartTime)
    WeightedNormalizedMuscleVelocity_Forward = np.concatenate((WeightedNormalizedMuscleVelocity_Forward,calculate_weighted_muscle_velocities_ISOLATED_ROTATION(A1[i],A2[i],Ȧ1[i],Ȧ2[i],dof=dof)),axis=0)
    NormalizedMuscleVelocity_Forward = np.concatenate((NormalizedMuscleVelocity_Forward,calculate_muscle_velocities(A1[i],A2[i],Ȧ1[i],Ȧ2[i])),axis=0)
print('\n')

plt.figure()
[plt.plot(t,WeightedNormalizedMuscleVelocity_Forward[:,i]) for i in OrderNumber]
ax1a = plt.gca()
ax1a.set_xlim(0,t_end*(1.3))
ax1a.set_title('Forward Direction')
if t_end!=1:
	ax1a.set_xlabel('Time (s)')
else:
	ax1a.set_xlabel('Normalized Time')
ax1a.set_ylabel('Normalized Muscle Velocity\nConcentric $\longleftrightarrow$ Eccentric')
[j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax1a.lines)]
ax1a.legend(OrderedMuscleList)

plt.figure()
[plt.plot(t,NormalizedMuscleVelocity_Forward[:,i]) for i in OrderNumber]
ax1b = plt.gca()
ax1b.set_xlim(0,t_end*(1.3))
ax1b.set_title('Forward Direction')
if t_end!=1:
	ax1b.set_xlabel('Time (s)')
else:
	ax1b.set_xlabel('Normalized Time')
ax1b.set_ylabel('Normalized Muscle Velocity\nConcentric $\longleftrightarrow$ Eccentric')
[j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax1b.lines)]
ax1b.legend(OrderedMuscleList)

EccentricCost_Forward = eccentric_cost(WeightedNormalizedMuscleVelocity_Forward,t_end=t_end,dt = dt,type='l1norm')
ConcentricCost_Forward = concentric_cost(WeightedNormalizedMuscleVelocity_Forward,t_end=t_end,dt = dt,type='l1norm')

# Reverse Direction

Ai_temp,Af_temp = Ai,Af
Af = Ai_temp
Ai = Af_temp

reaching_task_ISOLATED_ELBOW_TASK(dt = dt, Ai=Ai, Af=Af, t_end=t_end)
WeightedNormalizedMuscleVelocity_Reverse = calculate_weighted_muscle_velocities_ISOLATED_ROTATION(A1[0],A2[0],Ȧ1[0],Ȧ2[0],dof=dof)
NormalizedMuscleVelocity_Reverse = calculate_muscle_velocities(A1[0],A2[0],Ȧ1[0],Ȧ2[0])
StartTime = time.time()
statusbar(0,len(A1),Title = 'Reverse Vm',StartTime=StartTime)
for i in range(1,len(A1)):
    statusbar(i,len(A1),Title = 'Reverse Vm',StartTime=StartTime)
    WeightedNormalizedMuscleVelocity_Reverse = np.concatenate((WeightedNormalizedMuscleVelocity_Reverse,calculate_weighted_muscle_velocities_ISOLATED_ROTATION(A1[i],A2[i],Ȧ1[i],Ȧ2[i],dof=dof)),axis=0)
    NormalizedMuscleVelocity_Reverse = np.concatenate((NormalizedMuscleVelocity_Reverse,calculate_muscle_velocities(A1[i],A2[i],Ȧ1[i],Ȧ2[i])),axis=0)
print('\n')

plt.figure()
[plt.plot(t,WeightedNormalizedMuscleVelocity_Reverse[:,i]) for i in OrderNumber]
ax2a = plt.gca()
ax2a.set_xlim(0,t_end*(1.3))
ax2a.set_title('Reverse Direction')
if t_end!=1:
	ax2a.set_xlabel('Time (s)')
else:
	ax2a.set_xlabel('Normalized Time')
ax2a.set_ylabel('Normalized Muscle Velocity\nConcentric $\longleftrightarrow$ Eccentric')
[j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax2a.lines)]
ax2a.legend(OrderedMuscleList)

plt.figure()
[plt.plot(t,NormalizedMuscleVelocity_Reverse[:,i]) for i in OrderNumber]
ax2b = plt.gca()
ax2b.set_xlim(0,t_end*(1.3))
ax2b.set_title('Reverse Direction')
if t_end!=1:
	ax2b.set_xlabel('Time (s)')
else:
	ax2b.set_xlabel('Normalized Time')
ax2b.set_ylabel('Normalized Muscle Velocity\nConcentric $\longleftrightarrow$ Eccentric')
[j.set_color(OrderedColorsList[i]) for i,j in enumerate(ax2b.lines)]
ax2b.legend(OrderedMuscleList)

EccentricCost_Reverse = eccentric_cost(WeightedNormalizedMuscleVelocity_Reverse,t_end=t_end,dt = dt,type = 'l1norm')
ConcentricCost_Reverse = concentric_cost(WeightedNormalizedMuscleVelocity_Reverse,t_end=t_end,dt = dt,type = 'l1norm')

# Compare Costs

# plt.figure()
# plt.plot(np.arange(0,2*math.ceil(max([EccentricCost_Forward,ConcentricCost_Forward]))+1,1), np.arange(0,2*math.ceil(max([EccentricCost_Forward,ConcentricCost_Forward]))+1,1),'0.75',linestyle='--')
# plt.scatter([EccentricCost_Forward],[ConcentricCost_Forward],facecolor ='g', edgecolor = 'k', marker='o',s=100)
# plt.scatter([EccentricCost_Reverse],[ConcentricCost_Reverse],facecolor ='b', edgecolor = 'k', marker='o',s=100)
# ax3 = plt.gca()
# ax3.set_xlim(0,2*math.ceil(max([EccentricCost_Forward,EccentricCost_Reverse])))
# ax3.set_ylim(0,2*math.ceil(max([EccentricCost_Forward,EccentricCost_Reverse])))
# ax3.set_aspect('equal')
# ax3.set_title('Eccentric vs. Concentric Cost\nFor Forward and Reverse Movments')
# ax3.set_xlabel('Eccentric Cost (in $\hat{l}_{o}\cdot m$)')
# ax3.set_ylabel('Concentric Cost (in $\hat{l}_{o}\cdot m$)')
# ax3.legend(['Equal Cost Line','Forward Direction','Reversed Direction'])
plt.figure()
plt.bar(np.arange(2),[EccentricCost_Reverse,EccentricCost_Forward])
ax4 = plt.gca()
ax4.set_xticklabels(('Reverse','Forward'))
ax4.set_ylim(0,7.1)
ax4.set_yticks([0,3,6])
ax4.set_yticklabels(['0','','6'])
ax4.set_title('Eccentric Cost\nFor Forward and Reverse Movements')
ax4.set_ylabel('Sum of Afferent-Weighted Muscle Lengthening')
ax4.set_xticks([0,1])

plt.show()
