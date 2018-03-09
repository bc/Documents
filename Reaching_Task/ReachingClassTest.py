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
	import numpy as np
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
	TimeBreak = abs
	if StartTime != False:
		if i==0:
			time_array = []
			TimeLeft = '--'
		elif i==int(0.02*N):
			time_array.append(time.time()-StartTime)
			TimeLeft = '{0:1.1f}'.format(time_array[-1]*(N/(i+1)))
		elif i%int(0.02*N)==0:
			time_array.append(time.time()-StartTime)
			TimeLeft = '{0:1.1f}'.format(float(interpolate.interp1d(np.arange(len(time_array)),time_array,fill_value='extrapolate')(49))-time_array[-1])
		print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete, ' + '{0:1.1f}'.format(time.time() - StartTime) \
			+ 'sec, (est. ' + TimeLeft,' sec left)		\r', end='')
	else:
		print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete           \r',end = '')
def return_muscle_settings():
	"""
	Notes:
	Coefficients from observation, Ramsay, FVC, Holtzbaur, Pigeon, Kuechle, or Banks.

	BRA (Brachialis) EFE MA for Ramsay has R² = 0.990 whereas Pigeon has R² = 0.9988. Curve appears to be a better fit, as it experiences its smallest MA when Elbow angle = 0. Coefficients and equation number/type are listed below to test either implementation.

	src = 'Ramsay', eq = 1, Coefficients = [16.1991,-16.1463,24.5512,-6.3335,0], threshold = None
	src = 'Pigeon', dof = 'elbow', eq = None, Coefficients = Pigeon_coeff_conversion([5.5492,2.3080,2.3425,-2.0530,0,0]), threshold = None

	BRD (Brachioradialis) for Ramsay has R² = 0.988 whereas Pigeon has R² = 0.9989. Pigeon, however, only takes elbow angle into account, whereas Ramsay takes in variable PS angles. Coefficients and equation number/type are listed below to test either implementation.

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

	DELTp SFE MA is listed as -78.74 mm in Pigeon. Using Pigeon Coefficients convention, Kuechle (1997) has the DELTp MA for [-140,90] as Pigeon_coeff_conversion([ 2.28547177,  0.39721238, -0.33900829, -0.36146546,  0.,  0.]). This will yield a piecewise function that creates jumps with the new velocity formulation. Instead, we are going to try Pigeon_coeff_conversion([-2.38165173, -0.4486164 ,  0.58655808,  0.65003255, -0.82736695,0.20812998]) so that the function is within range and continuous during the ROM. Threshold (pi/2) has been removed for this new MA function.

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

	def Pigeon_coeff_conversion(Coefficients):
		"""
		Takes in Coefficient values from Pigeon (1996) -- which take in angles in degrees -- and coverts them into the properly scaled coefficients for radians, additionally scaled by the magnitude listed in the paper.

		Note that the coefficients listed in Pigeon (1996) are given in decending order (i.e., c₅,c₄,c₃,c₂,c₁,c₀). However to maintain continuity with the equations given in Ramsay (2009), we list coefficients in order of increasing power (i.e., c₀,c₁,c₂,c₃,c₄,c₅).
		"""
		import numpy as np
		assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
		assert type(Coefficients)==list, 'Coefficients must be a 6 element list.'
		Rad_Conversion = np.multiply(Coefficients,\
				np.array([1,(180/np.pi),(180/np.pi)**2,(180/np.pi)**3,(180/np.pi)**4,(180/np.pi)**5],dtype = 'float64'))
		new_Coefficients =\
			np.multiply(Rad_Conversion,np.array([1,1e-1,1e-3,1e-5,1e-7,1e-9],dtype='float64'))
		return(new_Coefficients)

	PC_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : [50.80,0,0,0,0,0],\
			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' :  295.6, \
		'Actual No' : 450, \
		'Corrected No' : 389.7, \
		'Relative Abundance' : 1.2,\
		'Optimal Muscle Length' : 150, \
		'Group' : 'flexor'}
	DELTa_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : Pigeon_coeff_conversion([ 12.7928795,  2.0480346,  0.8917734,  3.2207214, -2.3928223,  0.        ]),\
			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 355.7/3, \
		'Actual No' : 182/3, \
		'Corrected No' : 426.3/3, \
		'Relative Abundance' : 0.43,\
		'Optimal Muscle Length' : 98,\
		'Group' : 'flexor'}
	CB_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 20, \
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 39.8, \
		'Actual No' : 123, \
		'Corrected No' : 147.3, \
		'Relative Abundance' : 0.83,\
		'Optimal Muscle Length' : 93,\
		'Group' : 'flexor'}
	DELTp_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : Pigeon_coeff_conversion([-23.8165173, -4.486164 ,  5.8655808,  6.5003255, -8.2736695,2.0812998]), \
			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 355.7/3, \
		'Actual No' : 182/3, \
		'Corrected No' : 426.3/3, \
		'Relative Abundance' : 0.43,\
		'Optimal Muscle Length' : 137,\
		'Group' : 'extensor'}
	BIC_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : [29.21,0,0,0,0,0],\
			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : [8.4533,36.6147,2.4777,-19.432,2.0571,0,13.6502,0,0,-5.6172,0,-2.0854,0,0,0,0],\
			'Source' : 'Ramsay', 'Equation Number' : 2, 'Threshold' : 3.021, \
			'dof' : 'Elbow'}, \
		'Mass' : 163.8,\
		'Actual No' : 320,\
		'Corrected No' : 292.6,\
		'Relative Abundance' : 1.1,\
		'Optimal Muscle Length' : 116,\
		'Group' : 'flexor'}
	TRI_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : [-25.40,0,0,0,0,0], \
			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : [-24.5454,-8.8691,9.3509,-1.7518,0],\
			'Source' : 'Ramsay', 'Equation Number' : 1, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : (94.2+138.4+92.5), \
		'Actual No' : (200+222+98),\
		'Corrected No' : (223.7+269.6+221.8),\
		'Relative Abundance' : (0.89+0.82+0.44)/3,\
		'Optimal Muscle Length' : 134,\
		'Group' : 'extensor'}
	BRA_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : [16.1991,-16.1463,24.5512,-6.3335,0],\
			'Source' : 'Ramsay', 'Equation Number' : 1, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 141,\
		'Actual No' : 256,\
		'Corrected No' : 272.1,\
		'Relative Abundance' : 0.94,\
		'Optimal Muscle Length' : 86,\
		'Group' : 'flexor'}
	BRD_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0, \
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : 	[15.2564,-11.8355,2.8129,-5.7781,44.8143,0,2.9032,0,0,-13.4956,0,-0.3940,0,0,0,0],\
			'Source' : 'Ramsay', 'Equation Number' : 2, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 64.7,\
		'Actual No' : 70,\
		'Corrected No' : 190.2,\
		'Relative Abundance' : 0.37,\
		'Optimal Muscle Length' : 173,\
		'Group' : 'flexor'}
	PRO_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : 	[11.0405,-1.0079,0.3933,-10.4824,-12.1639,-0.4369,36.9174,3.5232,-10.4223,21.2604,-37.2444,10.2666,-11.0060,14.5974,-3.9919,1.7526,-2.0089,0.5460],\
			'Source' : 'Ramsay', 'Equation Number' : 3,'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 38.8, \
		'Actual No' : 187.6, \
		'Corrected No' : 185.5, \
		'Relative Abundance' : 1.3,\
		'Optimal Muscle Length' : 49,\
		'Group' : 'flexor'}
	FCR_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0, \
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : Pigeon_coeff_conversion([0.9351,0.5375,-0.3627,0,0,0]),\
			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : 2.86, \
			'dof' : 'Elbow'}, \
		'Mass' : 28.7, \
		'Actual No' : 129, \
		'Corrected No' : 125.7, \
		'Relative Abundance' : 1.0,\
		'Optimal Muscle Length' : 63,\
		'Group' : 'flexor'}
	ECRB_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : [-11.256,17.8548,1.6398,-0.5073,-2.8827,0,-0.0942,0,0,0,0,0,0,0,0,0],\
			'Source' : 'Ramsay', 'Equation Number' : 2, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 32.1, \
		'Actual No' : 102, \
		'Corrected No' : 132.7, \
		'Relative Abundance' : 0.77,\
		'Optimal Muscle Length' : 59,\
		'Group' : 'flexor'}
	ECRL_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : Pigeon_coeff_conversion([4.7304,1.2590,4.4347,-3.0229,0,0]),\
			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 44.3, \
		'Actual No' : 74, \
		'Corrected No' : 155.2, \
		'Relative Abundance' : 0.48,\
		'Optimal Muscle Length' : 81,\
		'Group' : 'flexor'}
	FCU_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : 5,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 36.5,\
		'Actual No' : 175,\
		'Corrected No' : 141.2,\
		'Relative Abundance' : 1.2,\
		'Optimal Muscle Length' : 51,\
		'Group' : 'flexor'}
	FDS_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : 5,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 95.2,\
		'Actual No' : 356,\
		'Corrected No' : 224.9,\
		'Relative Abundance' : 1.6,\
		'Optimal Muscle Length' : 84,\
		'Group' : 'flexor'}
	PL_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \

			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : 10,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : None, \
		'Actual No' : None, \
		'Corrected No' : None, \
		'Relative Abundance' : None,\
		'Optimal Muscle Length' : 64,\
		'Group' : 'flexor'}
	ECU_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : Pigeon_coeff_conversion([-2.1826,-1.7386,1.1491,0,0,0]),\
			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 25.2,\
		'Actual No' : 157,\
		'Corrected No' : 118,\
		'Relative Abundance' : 1.3,\
		'Optimal Muscle Length' : 62,\
		'Group' : 'extensor'}
	EDM_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0, \
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : -5, \
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 6.2, \
		'Actual No' : 53, \
		'Corrected No' : 59.8, \
		'Relative Abundance' : 0.89,\
		'Optimal Muscle Length' : 68,\
		'Group' : 'extensor'}
	EDC_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0, \
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : -5,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : 42.8, \
		'Actual No' : 219, \
		'Corrected No' : 152.6, \
		'Relative Abundance' : 1.4,\
		'Optimal Muscle Length' : 70,\
		'Group' : 'extensor'}
	AN_Settings = {\
		'Shoulder' : {\
			'MA Coefficients' : 0,\
			'Source' : 'Est', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Shoulder'}, \
		'Elbow' : {\
			'MA Coefficients' : Pigeon_coeff_conversion([-5.3450,-2.2841e-1,8.4297e-3,-14.329e-5,10.448e-7,-2.736e-9]),\
			'Source' : 'Pigeon', 'Equation Number' : None, 'Threshold' : None, \
			'dof' : 'Elbow'}, \
		'Mass' : None, \
		'Actual No' : None, \
		'Corrected No' : None, \
		'Relative Abundance' : None,\
		'Optimal Muscle Length' : None,\
		'Group' : 'extensor'}

	AllAvailableMuscles =[	"PC", "DELTa", "CB", "DELTp", "BIC", \
							"TRI", "BRA", "BRD", "PRO", "FCR",\
	 						"ECRB", "ECRL", "FCU", "FDS", "PL",\
	  						"ECU", "EDM", "EDC", "AN"]
	AllMuscleSettings = {	'PC': PC_Settings, 'DELTa' : DELTa_Settings, \
							'CB' : CB_Settings, 'DELTp' : DELTp_Settings,\
							'BIC' : BIC_Settings, 'TRI' : TRI_Settings, \
							'BRA' : BRA_Settings, 'BRD' : BRD_Settings, \
							'PRO' : PRO_Settings, 'FCR' : FCR_Settings, \
							'ECRB' : ECRB_Settings, 'ECRL' : ECRL_Settings, \
							'FCU' : FCU_Settings, 'FDS' : FDS_Settings, \
							'PL' : PL_Settings,'ECU' : ECU_Settings, \
							'EDM' : EDM_Settings, 'EDC' : EDC_Settings,\
							'AN' : AN_Settings}
	ValidResponse_1 = False
	while ValidResponse_1 == False:
		MuscleSelectionType = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nMuscle Selection:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n (1) - Default\n (2) - Custom\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nResponse: ")
		# import ipdb; ipdb.set_trace()
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
		if MuscleSelectionType not in ['1','2','']:
			print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInvalid Response! Please try again.')
			ValidResponse_1 = False
		elif MuscleSelectionType == '' or MuscleSelectionType == '1':
			for Muscle in ["PRO","AN"]:
				del(AllMuscleSettings[Muscle])
			ValidResponse_1 = True
		elif MuscleSelectionType == '2':
			ValidResponse_2 = False
			while ValidResponse_2 == False:
				MuscleListString = ""
				for i in range(len(AllAvailableMuscles)):
					MuscleListString += " (" + str(i+1) + ") - " + AllAvailableMuscles[i] + "\n"
				MuscleSelectionNumbers = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nSelect Muscle Number(s)\n(separated by commas & groups with hyphens):\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" + MuscleListString + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nMuscle Number(s): ")
				print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
				MuscleSelectionNumbers = [el.strip() for el in MuscleSelectionNumbers.split(",")]
				for el in MuscleSelectionNumbers:
					if "-" in el:
						temp = el.split("-")
						MuscleSelectionNumbers.remove(el)
						[MuscleSelectionNumbers.append(str(i)) \
										for i in range(int(temp[0]),int(temp[1])+1)]
				if np.array([el in [str(i+1) for i in range(len(AllAvailableMuscles))] \
									for el in MuscleSelectionNumbers]).all() == False:
					print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInvalid Response! Please check muscle numbers and try again.')
					ValidResponse_2 = False
				else:
					SelectedMuscles = [AllAvailableMuscles[int(el)-1] \
											for el in MuscleSelectionNumbers]
					MusclesToBeDeleted = [Muscle for Muscle in AllAvailableMuscles \
												if Muscle not in SelectedMuscles]
					for Muscle in MusclesToBeDeleted:
						del(AllMuscleSettings[Muscle])
					ValidResponse_2 = True
			ValidResponse_1 = True

	# MuscleList = AllMuscleSettings.keys()
	#
	# Rᵀ_symbolic = sp.Matrix([[MA_function(AllMuscleSettings[muscle][dof]) for \
	# 				dof in ['Shoulder','Elbow']] for muscle in MuscleList])
	# Ṙᵀ_symbolic = sp.Matrix(np.concatenate((sp.diff(Rᵀ_symbolic[:,0],q1),\
	# 											sp.diff(Rᵀ_symbolic[:,1],q2)),axis=1))
	#
	# Rᵀ_func = lambdify([q1,q2,q_PS],Rᵀ_symbolic)
	# Ṙᵀ_func = lambdify([q1,q2,q_PS],Ṙᵀ_symbolic)
	# # returns an (n,m) matrix when n is the number of muscles and m is the number of DOFS. We chose to return R.T because this is commonly utilized in muscle velocity calculations.
	return(AllMuscleSettings)
def return_MA_matrix_functions(AllMuscleSettings):
	import numpy as np
	import sympy as sp
	from sympy.utilities import lambdify
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

		src = Parameters['Source']
		Coefficients = Parameters['MA Coefficients']
		eq = Parameters['Equation Number']
		dof = Parameters['dof']
		threshold = Parameters['Threshold']

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

	MuscleList = AllMuscleSettings.keys()

	Rᵀ_symbolic = sp.Matrix([[MA_function(AllMuscleSettings[muscle][dof]) for \
					dof in ['Shoulder','Elbow']] for muscle in MuscleList])
	Ṙᵀ_symbolic = sp.Matrix(np.concatenate((sp.diff(Rᵀ_symbolic[:,0],q1),\
												sp.diff(Rᵀ_symbolic[:,1],q2)),axis=1))

	Rᵀ_func = lambdify([q1,q2,q_PS],Rᵀ_symbolic)
	Ṙᵀ_func = lambdify([q1,q2,q_PS],Ṙᵀ_symbolic)
	# returns an (n,m) matrix when n is the number of muscles and m is the number of DOFS. We chose to return R.T because this is commonly utilized in muscle velocity calculations.
	return(Rᵀ_func,Ṙᵀ_func)
def plot_MA_ranges(TrialData,ReturnFig=False):
	import numpy as np
	import matplotlib.pyplot as plt
	import sympy as sp
	from sympy.utilities import lambdify
	PS = np.pi/2
	A1 = np.linspace(-10*(np.pi/180),120*(np.pi/180),1001)
	A2 = np.linspace(0,160*(np.pi/180),1001)
	Rᵀ_func,_ = return_MA_matrix_functions(TrialData["All Muscle Settings"])
	Shoulder_MA_ranges = np.array(list(map(lambda a1: Rᵀ_func(a1,0,PS),A1)))[:,:,0]
	Elbow_MA_ranges = np.array(list(map(lambda a2: Rᵀ_func(np.pi/2,a2,PS),A2)))[:,:,1]
	fig, [ax1,ax2] = plt.subplots(1,2,figsize=(10,8))
	plt.subplots_adjust(hspace=0.2,bottom=0.2)
	for i in TrialData["Ordered Muscle Numbers"]:
		Shoulder_MA_plots = [ax1.plot(A1,\
								Shoulder_MA_ranges[:,TrialData["Ordered Muscle Numbers"][i]],\
									c=TrialData["Ordered Muscle Colors List"][i])[0] \
										for i in range(len(TrialData["All Muscle Settings"]))]
		ax1.set_title("Shoulder Moment Arms")
		ax1.set_ylabel("Moment Arm (mm)")
		ax1.set_xlabel("Shoulder Angle")
		Elbow_MA_plots = [ax2.plot(A2,\
								Elbow_MA_ranges[:,TrialData["Ordered Muscle Numbers"][i]],\
									c=TrialData["Ordered Muscle Colors List"][i])[0] \
										for i in range(len(TrialData["All Muscle Settings"]))]
		ax2.set_title("Elbow Moment Arms")
		ax2.set_xlabel("Elbow Angle")
	plt.figlegend(Shoulder_MA_plots,TrialData["Ordered Muscle List"],\
						loc='lower center',ncol=5,mode='expand')
	if ReturnFig == True: return(fig)
def input_notes():
	"""
	This function is designed to allow for a string of specific notes to be placed in the overall trial dict so that specific information regarding muscle parameters, MA functions, reaching directions, etc. can be documented so the trials are repeatable.
	"""
	Notes = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInput Notes for Specific Trial? (⏎ : None)\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
	if Notes == "":
		return(None)
	else:
		return(Notes)
def reach_type_prompt():
	import numpy as np
	import random
	DefaultSettings = False
	ValidResponse_1 = False
	while ValidResponse_1 == False:
		ReachTypeNumber = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nPlease select reaching movement number:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n (1) - Side-to-side\n (2) - Straight (Center)\n (3) - 45° Left\n (4) - 45° Right\n (5) - Fixed-target\n  ⏎  - Random\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nMovement Type: ")
		if ReachTypeNumber not in ['1','2','3','4','5','']:
			print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInvalid Response! Please try again.')
			ValidResponse_1 = False
		elif ReachTypeNumber == '':
			ReachTypeNumber = random.randint(0,3)
			DefaultSettings = True
			ValidResponse_1 = True
		else:
			ReachTypeNumber = int(ReachTypeNumber)-1
			print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
			ValidResponse_1 = True
	ReachType = ['Sideways','Center','Left','Right','Fixed-target'][ReachTypeNumber]
	DescriptiveTitle = ReachType + ' Reach'

	if DefaultSettings == False:
		ValidResponse_2 = False
		while ValidResponse_2 == False:
			RandomXiResponse = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nFix Starting Point?\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n (1) - True\n (2) - False\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nResponse: ")
			if RandomXiResponse not in ['1','2']:
				print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInvalid Response! Please try again.')
				ValidResponse_2 = False
			else:
				RandomXiResponse = int(RandomXiResponse)-1
				print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
				ValidResponse_2 = True
		RandomXiBool = [True,False][RandomXiResponse]

		ValidResponse_3 = False
		while ValidResponse_3 == False:
			RandomXfResponse = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nFix Ending Point?\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n (1) - True\n (2) - False\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nResponse: ")
			if RandomXfResponse not in ['1','2']:
				print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInvalid Response! Please try again.')
				ValidResponse_3 = False
			else:
				RandomXfResponse = int(RandomXfResponse)-1
				print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
				ValidResponse_3 = True
		RandomXfBool = [True,False][RandomXfResponse]
	else:
		RandomXiBool = True
		RandomXfBool = True

	if DescriptiveTitle == "Fixed-target Reach":
		"""
		This will allow for an input to be passed along to the reaching_movement.return_X_values(TrialData) in the DescriptiveTitle that will denote where along the arc of radius TrialData["Target Amplitude"] to start the movement. By convention, 0 radians will be associated with the 3 o'clock position and follow a counterclockwise rotation through 2π. Realistic starting positions will likely be in [π,2π].
		"""
		ValidResponse_4 = False
		while ValidResponse_4 == False:
			StartingPositionInRadians = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInput starting position along the arc with\nradius equal to the target amplitude as a \nmultiple of π. (Note: 0 radians corresp. to \n3 o'clock from  target.)\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nResponse (in rads): π ⨉ ")
			if any(ch.isalpha() for ch in StartingPositionInRadians) == True:
				print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInvalid Response (Numbers only)! Please try again.')
				ValidResponse_4 = False
			else:
				DescriptiveTitle = DescriptiveTitle + "_" + StartingPositionInRadians
				ValidResponse_4 = True
	return(DescriptiveTitle,ReachType,RandomXiBool,RandomXfBool)
def return_ordered_muscle_list_with_colors(TrialData):
	import numpy as np
	import matplotlib.pyplot as plt

	NumExtensors = sum([TrialData["All Muscle Settings"][key]['Group']=='extensor' \
							for key in TrialData["All Muscle Settings"].keys()])
	NumFlexors = sum([TrialData["All Muscle Settings"][key]['Group']=='flexor' \
							for key in TrialData["All Muscle Settings"].keys()])
	FlexorCMap=plt.get_cmap('autumn')
	ExtensorCMap = plt.get_cmap('YlGnBu')

	if NumFlexors == 1:
		FlexorColors = FlexorCMap(0.5)
	elif NumFlexors != 0:
		FlexorColors = iter(FlexorCMap(np.linspace(0,0.75,NumFlexors)))

	if NumExtensors == 1:
		ExtensorColors = ExtensorCMap(0.5)
	elif NumExtensors != 0:
		ExtensorColors = iter(list(reversed(ExtensorCMap(np.linspace(0.5,1,NumExtensors)))))

	FlexorColorsList = []
	ExtensorColorsList = []
	FlexorOrderedMuscleList = []
	ExtensorOrderedMuscleList = []
	for Muscle in TrialData["All Muscle Settings"]:
		if TrialData["All Muscle Settings"][Muscle]['Group'] == 'flexor':
			if NumFlexors == 1:
				FlexorColorsList.append(FlexorColors)
			else:
				FlexorColorsList.append(next(FlexorColors))
			FlexorOrderedMuscleList.append(Muscle)
		else:
			if NumExtensors == 1:
				ExtensorColorsList.append(ExtensorColors)
			else:
				ExtensorColorsList.append(next(ExtensorColors))
			ExtensorOrderedMuscleList.append(Muscle)


	OrderedColorsList = FlexorColorsList + list(reversed(ExtensorColorsList))
	OrderedMuscleList = FlexorOrderedMuscleList +\
	 											list(reversed(ExtensorOrderedMuscleList))
	OrderNumber = [list(TrialData["All Muscle Settings"].keys()).index(el)\
	 													for el in OrderedMuscleList]

	return(OrderNumber,OrderedMuscleList,OrderedColorsList)
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
	return(L1,L2)
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

	pp_2deriv()
	~~~~~~~~~~~~~~~~~~~

	Takes in X array and outputs the piecewise polynomial associated with the spline's second derivative.

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

	This checks to see if the maximum maximum value and the minimum mininum value calculated above will fall between y_min and y_max. This makes use of self.find_max_and_min()

	print_func()
	~~~~~~~~~~~~~~~~~~~

	This function uses pprint() to return a printout of the piecewise polynomial f(x).

	return_parameterized_X()
	~~~~~~~~~~~~~~~~~~~

	This function will return the parameterized x(t) and y(t) that follows path S along the curve y=f(x) subject to some tangential velocity profile dS/dt. This utilizes scipy.integrate.odeint to solve the time derivative to the path length function S = ∫(dx_dt)√(1 + f'(x)²)dt.

	return_parameterized_dX()
	~~~~~~~~~~~~~~~~~~~

	This function takes in x(t) and returns dx/dt and dy/dt derived from the relationship of arc length S, its time derivative dS/dt, and the path f(x(t)) and its derivative with respect to x.

	return_parameterized_d2X()
	~~~~~~~~~~~~~~~~~~~

	This function takes in x(t) and dx/dt and returns d²x/dt² and d²y/dt² derived from the relationship of arc length S, its first and second time derivatives, dS/dt and d²S/dt², respectively, and the path f(x(t)) and its first and second derivatives with respect to x.

	find_path_length()
	~~~~~~~~~~~~~~~~~~~

	Calculates the path length of the curve y=f(x) from S = ∫(dx_dt)√(1 + f'(x)²)dt. This is needed in order to describe the minimum jerk criterion tangential velocity equation, dS/dt.

	dS_dt()
	~~~~~~~~~~~~~~~~~~~

	Returns the minimum jerk criterion tangential velocity equation, dS/dt, given by:

							dS/dt = S*(30t² - 60t³ + 30t⁴)

	Where S is found from S = self.find_path_length()

	d2S_dt2()
	~~~~~~~~~~~~~~~~~~~

	Returns the minimum jerk criterion for acceleration along the path S, d²S/dt², given by:

							d²S/dt² = S*(60t - 180t² + 120t³)

	Where S is found from S = self.find_path_length()

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
		import numpy as np
		result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
			[lambda X: self.a[0] + self.b[0,0,0]*(X-self.x_initial) + self.c[0,0]*(X-self.x_initial)**2 + self.d[0,0,0]*(X-self.x_initial)**3, \
			lambda X: self.a[1] + self.b[1,0,0]*(X-self.x_break) + self.c[1,0]*(X-self.x_break)**2 + self.d[1,0,0]*(X-self.x_break)**3])
		return(result)
	def pp_deriv(self,X):
		import numpy as np
		result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
			[lambda X: self.b[0,0,0] + 2*self.c[0,0]*(X-self.x_initial) + 3*self.d[0,0,0]*(X-self.x_initial)**2, \
			lambda X: self.b[1,0,0] + 2*self.c[1,0]*(X-self.x_break) + 3*self.d[1,0,0]*(X-self.x_break)**2])
		return(result)
	def pp_2deriv(self,X):
		import numpy as np
		result = np.piecewise(X,[X <= self.x_break, X > self.x_break], \
			[lambda X: 2*self.c[0,0] + 6*self.d[0,0,0]*(X-self.x_initial), \
			lambda X: 2*self.c[1,0] + 6*self.d[1,0,0]*(X-self.x_break)])
		return(result)
	def find_max_and_min(self,x_min,x_max,y_min,y_max):
		def find_extrema():
			import numpy as np
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
			import numpy as np
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
		import numpy as np
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
	def return_parameterized_X(self,t_end=1):
		import scipy.integrate as integrate
		import numpy as np
		N = 1000
		t = np.linspace(0,t_end, N + 1)
		dt = t[1]-t[0]
		def ode_func(x,t,t_end):
			return(self.dS_dt(t,t_end=t_end)/np.sqrt(1 + self.pp_deriv(x)**2))
		X = integrate.odeint(lambda x,t: ode_func(x,t,t_end),self.xlim[0],t).flatten()
		Y = np.array(list(map(lambda x: self.pp_func(x),X)))
		dS = np.array(list(map(lambda dx,dy: np.sqrt(dx**2+dy**2),\
								np.gradient(X)/dt,np.gradient(Y)/dt)))
		assert sum(abs(self.dS_dt(t,t_end=t_end)-dS))/len(dS)<1e-4, "Error in parameterizing path to dS/dt. Check ODE func."
		return(X,Y)
	def return_parameterized_dX(self,x,t_end=1):
		import numpy as np
		N = 1000
		t = np.linspace(0,t_end,N + 1)
		# dt = t[1]-t[0]
		dS_dt = self.dS_dt(t,t_end=t_end)
		df_dx = self.pp_deriv(x)
		dx_dt = np.array(list(map(lambda dS_dt,df_dx: dS_dt/np.sqrt(1 + df_dx**2),dS_dt,df_dx)))
		dy_dt = df_dx*dx_dt
		return(dx_dt,dy_dt)
	def return_parameterized_d2X(self,x,dx_dt,t_end=1):
		import numpy as np
		N = 1000
		t = np.linspace(0,t_end,N + 1)
		# dt = t[1]-t[0]

		dS_dt = self.dS_dt(t,t_end=t_end)
		d2S_dt2 = self.d2S_dt2(t,t_end=t_end)

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
		import numpy as np
		import scipy.integrate as integrate
		return(integrate.quad(lambda x: \
				np.sqrt(1 + self.pp_deriv(x)**2),self.xlim[0],self.xlim[1])[0])
	def dS_dt(self,t,t_end=1):
		S_initial = 0
		S_final = self.find_path_length()
		return((S_final-S_initial)*(30*(t/t_end)**2 - 60*(t/t_end)**3 + 30*(t/t_end)**4)/t_end)
	def d2S_dt2(self,t,t_end=1):
		S_initial = 0
		S_final = self.find_path_length()
		return((S_final-S_initial)*(60*(t/t_end) - 180*(t/t_end)**2 + 120*(t/t_end)**3)/(t_end**2))
def generate_default_path(TrialData):
	import numpy as np
	def spline_coefficients(x1,x2,x3,y1,y2,y3,initial_slope,final_slope):
		"""
		Uses the values of (x1,y1), (x2,y2), and (x3,y3) to find the coefficients for the piecewise polynomial
		equation y = a + b*(x-x_o) + c*(x-x_o)**2 + d*(x-x_o)**3) for a clamped cubic spline with one break only.
		Returns coefficient arrays A, B, C,and D.
		"""
		import numpy as np
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

		C = c_coefficients(x1,x2,x3,y1,y2,y3,initial_slope,final_slope)
		D = d_coefficients(x1,x2,x3,C)
		B = b_coefficients(x1,x2,x3,y1,y2,y3,C,D)
		A = a_coefficients(y1,y2)

		assert test_b_coefficients(x1,x2,x3,y1,y2,y3,C,D,initial_slope), "Initial slope does not match the expected value"
		assert test_endpoint_slope(B[1,0,0],C[1,0],D[1,0,0],x2,x3,final_slope),"Problem with Endpoint Slope"
		assert test_for_discontinuity(A[0],B[0,0,0],C[0,0],D[0,0,0],x1,x2,A[1]), "Jump Discontinuity at t = %f!" %x2


		return(A,B,C[:2],D)
	ValidPath = False
	while ValidPath==False:
		DefaultXi = 0.05 # m DISPLACEMENT TO BE TRANSLATED BACK LATER
		DefaultXf = TrialData["Target Amplitude"] + 0.05 # m
		DefaultYi = 0 # m
		DefaultYf = 0 # m
		EndpointErrorSigma = 0.0025
		EndpointErrorTheta = 30*(np.pi/180) # Allowable error in initial/final slope, tan(ϑ)
		MaximumDeviationInY = 0.02 # m

		if TrialData["Randomize Boundary Positions"][0] == False:
			x_initial = DefaultXi + np.random.normal(0,EndpointErrorSigma) # cm
			y_initial = DefaultYi + np.random.normal(0,EndpointErrorSigma) # cm
		else:
			x_initial = DefaultXi # cm
			y_initial = DefaultYi # cm

		if TrialData["Randomize Boundary Positions"][1] == False:
			x_final = DefaultXf + np.random.normal(0,EndpointErrorSigma) # cm
			y_final = DefaultYf + np.random.normal(0,EndpointErrorSigma) # cm
		else:
			x_final = DefaultXf # cm
			y_final = DefaultYf # cm

		initialerror = np.random.uniform(-np.tan(EndpointErrorTheta),np.tan(EndpointErrorTheta))
		finalerror = -np.sign(initialerror)*\
				abs(np.random.uniform(-np.tan(EndpointErrorTheta),np.tan(EndpointErrorTheta)))

		if initialerror>0:
			ymax = max([y_initial,y_final])+MaximumDeviationInY
			ymin = 0
		else:
			ymax = 0
			ymin = min([y_initial,y_final])-MaximumDeviationInY

		xmax = x_final
		xmin = x_initial
		x_rand = np.random.normal((x_initial+x_final)/2,abs(x_initial+x_final)/4)
		y_rand = np.random.normal((ymax+ymin)/2,abs(ymax+ymin)/4)

		A,B,C,D = spline_coefficients(x_initial,x_rand,x_final,y_initial,y_rand,y_final,initialerror,finalerror)
		path = Spline(A,B,C,D,x_initial,x_rand,x_final)
		if path.is_within_bounds(x_initial,x_final, ymin, ymax):
			ValidPath = True
		else:
			ValidPath = False
	return(path)
class reaching_movement:
	def __init__(self,DefaultPath):
		self.DefaultPath = DefaultPath
	def return_X_values(self,TrialData):
		"""
		This takes in a string -- either 'Center','Right','Left', or 'Sideways' -- and returns the necessary initial and final positions for the movement, based on Flash/Hogan (1985).

		Parameters:
		# 0.80*(L1+L2) = 0.4586
		# 0.35/(L1+L2) = 0.5295
		# Sternum at -0.177 = -L1*(0.129/0.186) <-- Anthropomorphic Ratio
		"""
		import numpy as np
		t_end = TrialData["Movement Duration"]
		TargetAmplitude = TrialData["Target Amplitude"]
		L1,_ = TrialData["Limb Lengths"]
		x,y = self.DefaultPath.return_parameterized_X(t_end=t_end)
		ẋ,ẏ = self.DefaultPath.return_parameterized_dX(x,t_end=t_end)
		ẍ,ÿ = self.DefaultPath.return_parameterized_d2X(x,ẋ,t_end=t_end)

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

		MedianPlane = -L1*(0.129/0.186) # Sternum at -L1*(0.129/0.186) <-- Anthropomorphic Ratio
		DefaultDisplacement_x = 0.05

		assert TrialData["Reach Type"][:-6].capitalize() in\
		 			['Center','Right','Left','Sideways','Fixed-target'], \
						"ReachType must be either 'Center','Right','Left', or 'Sideways'."

		if TrialData["Reach Type"][:-6].capitalize() == 'Sideways':
			"""
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			~~~~~~~~~~~~~~~~~~    Side to Side   ~~~~~~~~~~~~~~~~~~
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Xi = [-MedianPlane-TargetAmplitude/2,0.1+TargetAmplitude/2]
			Xf = [-MedianPlane+TargetAmplitude/2,0.1+TargetAmplitude/2]
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			"""
			x,y = translate_xy(x,y,px=-DefaultDisplacement_x)
			# x,y = rotate_xy(x,y,0)
			x,y = translate_xy(x,y, px=MedianPlane-TargetAmplitude/2, py=0.10+TargetAmplitude/2)

			# ẋ,ẏ = ẋ,ẏ

			# ẍ,ÿ = ẍ,ÿ

		elif TrialData["Reach Type"][:-6].capitalize() == 'Center':
			"""
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			~~~~~~~~~~~~~~~~~    Center Reach     ~~~~~~~~~~~~~~~~~
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Xi = [-MedianPlane,0.20]
			Xf = [-MedianPlane,0.20 + TargetAmplitude]
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			"""
			x,y = translate_xy(x,y,px=-DefaultDisplacement_x)
			x,y = rotate_xy(x,y,np.pi/2)
			x,y = translate_xy(x,y, px=MedianPlane, py=0.20)

			ẋ,ẏ = rotate_xy(ẋ,ẏ,np.pi/2)

			ẍ,ÿ = rotate_xy(ẍ,ÿ,np.pi/2)

		elif TrialData["Reach Type"][:-6].capitalize() == 'Left':
			"""
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			~~~~~~~~~~~~~~~~~ Left Diagonal Reach ~~~~~~~~~~~~~~~~~
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Xi = [-MedianPlane,0.20]
			Xf = [-MedianPlane-TargetAmplitude/(2**0.5), 0.20 + TargetAmplitude/(2**0.5)]
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			"""
			x,y = translate_xy(x,y,px=-DefaultDisplacement_x)
			x,y = rotate_xy(x,y,3*np.pi/4)
			x,y = translate_xy(x,y, px=MedianPlane, py=0.20)

			ẋ,ẏ = rotate_xy(ẋ,ẏ,3*np.pi/4)

			ẍ,ÿ = rotate_xy(ẍ,ÿ,3*np.pi/4)

		elif TrialData["Reach Type"][:-6].capitalize() == 'Right':
			"""
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			~~~~~~~~~~~~~~~~ Right Diagonal Reach ~~~~~~~~~~~~~~~~~
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Xi = [-MedianPlane,0.20]
			Xf = [-MedianPlane+TargetAmplitude/(2**0.5), 0.20 + TargetAmplitude/(2**0.5)]
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			"""
			x,y = translate_xy(x,y,px=-DefaultDisplacement_x)
			x,y = rotate_xy(x,y,np.pi/4)
			x,y = translate_xy(x,y,px = MedianPlane, py=0.20)

			ẋ,ẏ = rotate_xy(ẋ,ẏ,np.pi/4)

			ẍ,ÿ = rotate_xy(ẍ,ÿ,np.pi/4)

		elif TrialData["Reach Type"][:-6].capitalize() == 'Fixed-target':
			"""
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			~~~~~~~~~~~~~~~~~ Fixed Target Reach ~~~~~~~~~~~~~~~~~~
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Xi = [-MedianPlane,0.20]
			Xf = [-MedianPlane,0.20 + TargetAmplitude]
			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			"""
			ϕ = TrialData["Reach Angle"] - np.pi
			TargetAmplitude = TrialData["Target Amplitude"]
			x,y = translate_xy(x,y,px=-DefaultDisplacement_x)
			x,y = rotate_xy(x,y,ϕ)
			px = MedianPlane - TargetAmplitude*np.cos(ϕ)
			py = 0.20 + TargetAmplitude - TargetAmplitude*np.sin(ϕ)
			x,y = translate_xy(x,y,px=px,py=py)

			ẋ,ẏ = rotate_xy(ẋ,ẏ,ϕ)

			ẍ,ÿ = rotate_xy(ẍ,ÿ,ϕ)

		X = np.concatenate([x[np.newaxis,:],y[np.newaxis,:]],axis=0)
		Ẋ = np.concatenate([ẋ[np.newaxis,:],ẏ[np.newaxis,:]],axis=0)
		Ẍ = np.concatenate([ẍ[np.newaxis,:],ÿ[np.newaxis,:]],axis=0)
		return(X,Ẋ,Ẍ)
	def inverse_kinematics(self,X,TrialData):
		"""
		Takes in a (2,N) list/array with values for x and y (endpoint) and maps to the  A1 and A2 with the current angles (in radians) from the inverse kinematics.
		"""
		import numpy as np
		from math import acos, sin, cos, atan, atan2
		from numpy import pi
		assert np.shape(X)[0]==2, "X must be either a (2,len(X)) list/array"
		x,y = np.split(X,2,axis=0) # Returns 2 (1,N) arrays
		L1,L2 = TrialData["Limb Lengths"]
		"""
		Math Logic:
		x² + y² > (L₁+L₂)² = L₁² + 2L₁L₂ + L₂² > L₁² + L₂²
		x² + y² - L₁² - L₂² > 0
		a₂ = cos⁻¹((x² + y² - L₁² - L₂²)/2L₁L₂) ∊ (0,π/2)

		Map functions take in (N,1) list/arrays -- i.e., only length-1 arrays can be converted to Python scalars needed for the map function. Therefore, the transpose is needed to change (1,N) to (N,1)
		"""
		a2 = lambda x,y: acos((x**2 + y**2 - L1**2 - L2**2)/(2*L1*L2))
		a1 = lambda x,y,a2: atan2(y,x) - atan2(L2*sin(a2),(L1+L2*cos(a2)))
		self.A2 = np.array(list(map(a2,x.T,y.T)), dtype='float64', ndmin=2) # Returns a (1,N) array
		self.A1 = np.array(list(map(a1,x.T,y.T,self.A2.T)), dtype='float64', ndmin=2) # Returns a (1,N) array
	def update_angular_velocity(self,X,Ẋ,TrialData):
		import numpy as np
		if hasattr(self,'A1') == False:
			self.inverse_kinematics(X,TrialData)
		def inverse_jacobian(a1,a2):
			import numpy as np
			from math import cos, sin
			L1,L2 = TrialData["Limb Lengths"]
			det_J = L1*L2*cos(a1 + a2)*sin(a1) - L1*L2*sin(a1 + a2)*cos(a1)
			J_inv = (1/det_J)*np.matrix([[- L2*cos(a1 + a2),- L2*sin(a1 + a2)],\
							[L1*cos(a1)+L2*cos(a1 + a2),	L2*sin(a1 + a2) + L1*sin(a1)]], \
										dtype = 'float64')
			return(J_inv)
		"""
		Note: numpy.float64 does not support inverse matrix operators. Luckily, this is an always invertible 2x2 matrix (J is nonsingular within the ROM because a₂ ∊ (0,π/2)), so the inverse Jacobian function has been mapped instead.
		"""
		# J = list(map(jacobian,A1.T,A2.T)) # Returns a list of shape (N,2,2)
		J_inv = list(map(inverse_jacobian,self.A1.T,self.A2.T)) # Returns a list of shape (N,2,2)
		"""
		In order to map properly, X (shape (2,N)) must be be first split into N (2,1) arrays (length-1). Therefore, the third argument to map() does not need to be transposed. Similar logic follows for why J_inv is not transposed, as it is has N (2,2) matrix arrays. Changing the output to an array aids in creating arrays Ȧ1 and Ȧ2.
		"""
		Ȧ = list(map(lambda j_inv,ẋ: np.array(j_inv*ẋ), \
							J_inv, np.split(Ẋ,np.shape(Ẋ)[1],axis=1)))
		Ȧ1,Ȧ2 = np.split(np.concatenate(Ȧ, axis=1), 2, axis=0) # both are (1,N) arrays
		return(Ȧ1,Ȧ2)
	def update_angular_acceleration(self,Ȧ1,Ȧ2,Ẋ,Ẍ,TrialData):
		from math import cos,sin
		import numpy as np

		assert np.shape(Ẋ)[0]==2, "Ẋ must be a (2,len(Ẋ)) list/array"
		assert np.shape(Ẍ)[0]==2, "Ẍ must be a (2,len(Ẍ)) list/array"
		ẍ,ÿ = np.split(Ẍ,2,axis=0)
		L1,L2 = TrialData["Limb Lengths"]

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

		Ä1 = np.array(list(map(ä1,self.A1.T,self.A2.T,Ȧ1.T,Ȧ2.T,Ẋ[0].T,Ẋ[1].T,Ẍ[0].T,Ẍ[1].T))\
						,dtype = 'float64',ndmin=2).T # returns a (1,N) array
		Ä2 = np.array(list(map(ä2,self.A1.T,self.A2.T,Ȧ1.T,Ȧ2.T,Ẋ[0].T,Ẋ[1].T,Ẍ[0].T,Ẍ[1].T))\
						,dtype = 'float64',ndmin=2).T # returns a (1,N) array
		return(Ä1,Ä2)
	def update_angle_lists(self,X,Ẋ,Ẍ,TrialData):
		"""
		Takes in three (2,N) endpoint arrays and returns global lists for angles 1 and 2 of shape (1,N).
		"""
		import numpy as np
		self.inverse_kinematics(X,TrialData)
		Ȧ1,Ȧ2 = self.update_angular_velocity(X,Ẋ,TrialData)
		Ä1,Ä2 = self.update_angular_acceleration(Ȧ1,Ȧ2,Ẋ,Ẍ,TrialData)
		return(Ȧ1,Ȧ2,Ä1,Ä2)
	def reaching_task_kinematics(self,TrialData,EOM='Zadravec'):
		# def set_link_lengths(New_L1=None,New_L2=None,EOM = "Zadravec"):
		# 	"""
		# 	Sets the global link lengths for a 2 DOF planar reaching task. Changes the values of the link lengths. New_L1 and New_L2 must be numbers. Set EOM to Uno of Zadravec.
        #
		# 	Default L1 and L2 values are calculated from Winter's Anthropomorphic measurement scales for a 72 inch tall subject. L2 is currently the sum of forearm and hand measurements with some added value to reflect the reaching apparatus.
		# 	"""
		# 	if New_L1 != None:
		# 		assert type(New_L1) == int or type(New_L1)==float, "New_L1 must be a number."
		# 	if New_L2 != None:
		# 		assert type(New_L2) == int or type(New_L2)==float, "New_L2 must be a number."
		# 	assert EOM in [None,"Uno","Zadravec"], "EOM can be either None, 'Uno', or 'Zadravec'"
		# 	global L1, L2
		# 	# 72
		# 	Height_inches = 64.5
		# 	Height = 2.54*Height_inches/100
		# 	if New_L1 == None and New_L2 == None:
		# 		if EOM == "Uno":
		# 			L1 = 0.256
		# 			L2 = 0.315
		# 		elif EOM == "Zadravec":
		# 			L1 = 0.298
		# 			L2 = 0.419
		# 		else:
		# 			L1 = Height*0.186
		# 			L2 = Height*(0.146+0.108)
		# 	elif New_L1 != None:
		# 		L1 = New_L1
		# 		L2 = Height*(0.146+0.108)
		# 	elif New_L2 != None:
		# 		L1 = Height*0.186
		# 		L2 = New_L2
		import numpy as np
		assert TrialData["Reach Type"][:-6].capitalize() in \
			['Center','Right','Left','Sideways','Fixed-target'], \
				"ReachType must be either 'Center','Right','Left', or 'Sideways'."
		set_link_lengths()
		X,Ẋ,Ẍ = self.return_X_values(TrialData)
		Ȧ1,Ȧ2,Ä1,Ä2=self.update_angle_lists(X,Ẋ,Ẍ,TrialData)
		# calculate_torques(EOM=EOM)
		return(Ȧ1,Ȧ2,Ä1,Ä2,X,Ẋ,Ẍ)
	def return_MA_matrix(self,TrialData):
		"""
		Notes:

		The angle of pronation is fixed at pi/2 for this reaching paradigm.
		The angle of radial/ulnar deviation is set to zero for the fixed wrist apparatus of the model paradigm.
		These functions have been verified to match the previous posture dependent MAs from Ramsey (2010) - SEE ERRATUM
		"""
		import numpy as np
		from numpy import pi
		import sympy as sp
		from sympy.utilities import lambdify
		def return_MA_matrix_functions(AllMuscleSettings):
			import numpy as np
			import sympy as sp
			from sympy.utilities import lambdify
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

				src = Parameters['Source']
				Coefficients = Parameters['MA Coefficients']
				eq = Parameters['Equation Number']
				dof = Parameters['dof']
				threshold = Parameters['Threshold']

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

			MuscleList = AllMuscleSettings.keys()

			Rᵀ_symbolic = sp.Matrix([[MA_function(AllMuscleSettings[muscle][dof]) for \
							dof in ['Shoulder','Elbow']] for muscle in MuscleList])
			Ṙᵀ_symbolic = sp.Matrix(np.concatenate((sp.diff(Rᵀ_symbolic[:,0],q1),\
														sp.diff(Rᵀ_symbolic[:,1],q2)),axis=1))

			Rᵀ_func = lambdify([q1,q2,q_PS],Rᵀ_symbolic)
			Ṙᵀ_func = lambdify([q1,q2,q_PS],Ṙᵀ_symbolic)
			# returns an (n,m) matrix when n is the number of muscles and m is the number of DOFS. We chose to return R.T because this is commonly utilized in muscle velocity calculations.
			return(Rᵀ_func,Ṙᵀ_func)

		if hasattr(self,'A1') == False:
			X,_,_ = self.return_X_values(TrialData)
			self.inverse_kinematics(X,TrialData)

		Rᵀ_func,Ṙᵀ_func = return_MA_matrix_functions(TrialData["All Muscle Settings"])
		# Rᵀ_func,Ṙᵀ_func = [TrialData[key] for key in \
		# 						["Moment Arm Matrix Function","Moment Arm Matrix Derivative"]]

		Rᵀ = np.array(list(map(lambda A1,A2: \
							np.float64(Rᵀ_func(A1,A2,pi/2).T),\
							self.A1.T,self.A2.T)))
		Ṙᵀ = np.array(list(map(lambda A1,A2:\
							np.float64(Ṙᵀ_func(A1,A2,pi/2).T),\
							self.A1.T,self.A2.T)))
		# returns two matrices of size (N,2,m)
		return(Rᵀ,Ṙᵀ)
	def calculate_muscle_velocities(self,Ȧ1,Ȧ2,Rᵀ,Ṙᵀ,TrialData):
		"""
		v = (-dR/dϑ)⋅(dϑ/dt)⋅ϑ + (-R)⋅(dϑ/dt)
		(dsⱼ/dϑ₁)⋅(dϑ₁/dt) = √((dϑ₁/dt)²⋅(dr₁ⱼ/dϑ₁)² + (dϑ₁/dt)²⋅(r₁ⱼ)²)
		(dsⱼ/dt) = sum([√(el) for el in diag(dϑ/dt)(J²_{Rⱼ})(dϑ/dt) + diag(Rⱼ)(diag(dϑ/dt)²)(dϑ₁/dt)])

		Note 02/27/18: The last line for the above equation is correct, but it is not utilized. Instead it is a calculated in one long matrix multiplication.

		Need to hard code whether the rotation causes lengthening or shortening. This is determined by the sign of -r₁ⱼ⋅dϑ₁/dt for each d.o.f. CORRECTION: by removing the angular velocity from the square root, the sign of the rotation is preserved -- regaining the sign convention.

		Note 02/27/18: The sign still needs to hardcoded to reflect the line of action of the moment arm relative to the angular velocity (positive or negative torque production).
		"""

		import numpy as np
		# global Ȧ1,Ȧ2,AllMuscleSettings,MomentArmMatrix,dMomentArmMatrix

		# MuscleList = list(AllMuscleSettings.keys())
		Ȧ = np.swapaxes(np.array(np.concatenate((Ȧ1,Ȧ2),axis=0),ndmin=3),0,2)

		"""
		Note: Ȧ is a 3D array of size (N,2,1). Thus, Ȧ.T has shape (1,2,N). Numpy.ndarray multiplication of arrays with IDENTICAL shapes will result in an element by element multiplication similar to np.multiply. Thus we multiply -Ȧ.T[0] by the (2,N) array that is the result of the norm of dMomentArmMatrix.T[j] and MomentArmMatrix.T[j] (each of shape (2,N)). This is equivalent to -ȧ₁*√(ṙ₁ⱼ² + r₁ⱼ²) for all t.

		Normalizing these muscle velocities can similarly be done by dividing the (m,1,N) MuscleVelocity array by the OptimalMuscleLength array of shape (m,1,1).
		"""

		MuscleVelocity = \
						(abs(Ȧ.T)*(((Ṙᵀ.T**2)+(Rᵀ.T**2))**0.5)\
							*np.sign(-Ȧ.T*Rᵀ.T))\
								.sum(axis=1)[:,np.newaxis]

		OptimalMuscleLength = np.array(\
		[TrialData["All Muscle Settings"][key]['Optimal Muscle Length'] \
			for key in TrialData["All Muscle Settings"].keys()],\
				dtype='float64')[:,np.newaxis,np.newaxis]
		NormalizedMuscleVelocity = MuscleVelocity/OptimalMuscleLength
		return(NormalizedMuscleVelocity)
	def calculate_weighted_muscle_velocities(self,Rᵀ,NormalizedMuscleVelocity,TrialData):

		import numpy as np

		"""
		Next we need to multiply each muscle velocity, AT EACH TIMESTEP, by the appropriate scaling factor, which may be a function of joint angle (i.e., movement step number).
		"""

		CorrectedAfferentNumber = np.array(\
			[TrialData["All Muscle Settings"][key]['Corrected No'] \
				for key in TrialData["All Muscle Settings"].keys()],\
					dtype='float64')[:,np.newaxis]
		WeightedMuscleVelocity = (NormalizedMuscleVelocity*abs(Rᵀ).T/1000).sum(axis=1)  \
									/((abs(Rᵀ)>0).T.sum(axis=1))  \
										*CorrectedAfferentNumber

		"""
		If you wanted to change the weighting scheme, use the above notation whereby the numpy.ndarrays have the same primary shape (i.e., length). Removing the .sum(axis=1) for the first line of the WeightedMuscleVelocity equation will return the individual weighted muscle velocity components per joint.
		"""
		return(WeightedMuscleVelocity)
	def return_muscle_velocities(self,TrialData,Weighted=False):
		"""
		Returns an (n,N) array where n is the number of muscles and N is the length of the array. Weighted can be set to True to return the weighted muscle velocities.
		"""
		Ȧ1,Ȧ2,Ä1,Ä2,X,Ẋ,Ẍ = self.reaching_task_kinematics(TrialData)
		Rᵀ,Ṙᵀ = self.return_MA_matrix(TrialData)
		MuscleVelocities = self.calculate_muscle_velocities(Ȧ1,Ȧ2,Rᵀ,Ṙᵀ,TrialData)
		if Weighted == True:
			MuscleVelocities =\
			 	self.calculate_weighted_muscle_velocities(Rᵀ,MuscleVelocities,TrialData)
		return(MuscleVelocities.squeeze())
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
def animate_plots(response,Movement,TrialData,Weighted=False, save_as_gif = False):
	assert type(response)==bool, "Input must be either True or False."
	assert type(Weighted)==bool, "Weighted must be either True or False."
	assert TrialData["Reach Type"][:-6] in ['Sideways','Center','Left','Right','Fixed-target'], "ReachType must be either 'Sideways','Center','Left', or 'Right'"

	if response == True:
		import numpy as np
		import matplotlib.pyplot as plt
		from matplotlib.patches import Ellipse
		import matplotlib.animation as animation
		import matplotlib.patches as patches

		t_end = TrialData["Movement Duration"]
		N = 1000
		t = np.linspace(0,t_end, N + 1)
		L1,L2 = TrialData["Limb Lengths"]
		MuscleVelocities = Movement.return_muscle_velocities(TrialData,Weighted=Weighted)

		A1_forward = Movement.A1
		A2_forward = Movement.A2
		A1_reverse = np.array(list(reversed(Movement.A1.T))).T
		A2_reverse = np.array(list(reversed(Movement.A2.T))).T

		Vm_forward = MuscleVelocities
		Vm_reverse = -np.array(list(reversed(Vm_forward.T))).T

		MedianPlane = L1*(0.129/0.186)

		fig, ((ax5,ax6),(ax1,ax2),(ax3,ax4)) = plt.subplots(3,2,figsize=(10,8))
		plt.subplots_adjust(top=0.9,hspace=0.2,bottom=0.2,left=0.2)

		if TrialData["Reach Type"][:-6] == 'Left':
			DescriptiveTitle = "45$^\circ$ Reach Left\n"
		elif TrialData["Reach Type"][:-6] == 'Right':
			DescriptiveTitle = "45$^\circ$ Reach Right\n"
		elif TrialData["Reach Type"][:-6] == 'Sideways':
			DescriptiveTitle = "Side-to-side Reach\n"
		elif TrialData["Reach Type"][:-6] == 'Center':
			DescriptiveTitle = "Straight Forward (Center) Reach\n"
		elif TrialData["Reach Type"][:-6] == 'Fixed-target':
			DescriptiveTitle = "Fixed Target Reach\n"

		# if Weighted==True:
		# 	TypeString = "\n(Afferent-Weighted $\hat{v}_m$)\n"
		# else:
		# 	TypeString = "\n(Normalized $\hat{v}_m$)\n"
		if t_end == 1:
			MovementDurationString = "Movement Duration : " + str(t_end) \
									+ " sec\n"
		else:
			MovementDurationString = "Movement Duration : " + str(t_end) \
									+ " secs\n"
		plt.suptitle(DescriptiveTitle + MovementDurationString,Fontsize=20,y=0.975)

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
		# import ipdb; ipdb.set_trace()
		ax3.set_xlim(0,t_end)
		ax3.set_xticks([0,t_end])
		ax3.set_xticklabels(['Start','Finish'])
		if Weighted == True:
			ax3.set_ylabel('Afferent-Weighted $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
		else:
			ax3.set_ylabel('Normalized $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
		if np.shape(Vm_forward) == (len(Vm_forward),):
			NormalizedVmPlots_forward = ax3.plot(t,Vm_forward,c=TrialData["Ordered Muscle Colors List"][0])
			plt.figlegend(NormalizedVmPlots_forward,TrialData["Ordered Muscle List"],\
							loc='lower center',ncol=5,mode='expand')
			bounds = max([Vm_forward.max(),Vm_reverse.max()])
		else:
			NormalizedVmPlots_forward = [ax3.plot(t.T,Vm_forward[j].T) \
												for j in TrialData["Ordered Muscle Numbers"]]
			[k.set_color(TrialData["Ordered Muscle Colors List"][j]) \
												for j,k in enumerate(ax3.lines)]

			bounds = max([	max([max(Vm_forward[i]) \
								for i in range(len(TrialData["Ordered Muscle List"]))]),\
							max([max(Vm_reverse[i]) \
								for i in range(len(TrialData["Ordered Muscle List"]))]) ] )
			plt.figlegend([el[0] for el in NormalizedVmPlots_forward],TrialData["Ordered Muscle List"],loc='lower center',ncol=5,mode='expand')

		Screen_f = plt.Rectangle((0,-13),t_end,26,color='w')
		ax3.add_patch(Screen_f)
		ax3.spines['right'].set_visible(False)
		ax3.spines['top'].set_visible(False)



		ax3.set_ylim([-1.1*bounds,1.1*bounds])
		# plt.figlegend([el[0] for el in NormalizedVmPlots_forward],TrialData["Ordered Muscle List"],loc='lower center',ncol=5,mode='expand')

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


		#Might need to add t,WeightedNormalizedMuscleVelocity_Forward,TrialData["Ordered Muscle Numbers"],etc
		ax4.set_xlim(0,t_end)
		ax4.set_ylim(-12,12)
		ax4.set_xticks([0,t_end])
		ax4.set_xticklabels(['Start','Finish'])

		if np.shape(Vm_reverse) == (len(Vm_reverse),):
			NormalizedVmPlots_reverse = ax4.plot(t,Vm_reverse,c=TrialData["Ordered Muscle Colors List"][0])
		else:
			NormalizedVmPlots_reverse = [ax4.plot(t.T,Vm_reverse[j].T) \
												for j in TrialData["Ordered Muscle Numbers"]]
			[k.set_color(TrialData["Ordered Muscle Colors List"][j]) \
												for j,k in enumerate(ax4.lines)]
        #
		# NormalizedVmPlots_reverse = [ax4.plot(t.T,Vm_reverse[j].T) \
		# 									for j in TrialData["Ordered Muscle Numbers"]]

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
			Screen_f.xy = (t_end*i/len(t.T),-13)
			Screen_f._width = t_end*(1 - i/len(t.T))

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
			Screen_r.xy = (t_end*i/len(t.T),-13)
			Screen_r._width = t_end*(1 - i/len(t.T))
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
		# if save_as_gif:
		# 	ani.save('test.gif', writer='imagemagick', fps=30)
		plt.show()
def create_trial_data(NumberOfTrials=100):
	"""
	Returns dict TrialData with NumberOfTrials properly oriented reaching movements.
	"""
	import numpy as np
	import matplotlib.pyplot as plt

	DescriptiveTitle,ReachType,RandomXiBool,RandomXfBool = reach_type_prompt()
	AllMuscleSettings = return_muscle_settings()
	Notes = input_notes()

	TrialData = {	"All Muscle Settings" : AllMuscleSettings, \
					"Notes" : Notes, \
					"Equations of Motion" : "Zadravec", \
					"Target Amplitude" : 0.35,\
					"Movement Duration" : 1,\
					"Randomize Boundary Positions" : [RandomXiBool,RandomXfBool],\
					"Default Paths" : []}
	if DescriptiveTitle[:5]=="Fixed":
		TrialData["Reach Type"] = DescriptiveTitle[:18]
		TrialData["Reach Angle"] = np.pi*eval(DescriptiveTitle[19:])
	else:
		TrialData["Reach Type"] = DescriptiveTitle

	OrderNumber, OrderedMuscleList, OrderedColorsList = \
	 												return_ordered_muscle_list_with_colors(TrialData)
	TrialData["Ordered Muscle List"] = OrderedMuscleList
	TrialData["Ordered Muscle Colors List"] = OrderedColorsList
	TrialData["Ordered Muscle Numbers"] = OrderNumber
	TrialData["Limb Lengths"] = list(set_link_lengths(EOM=TrialData["Equations of Motion"]))

	for i in range(NumberOfTrials):
		DefaultPath = generate_default_path(TrialData)
		TrialData["Default Paths"].append(DefaultPath)

	return(TrialData)
def animate_random_trajectory(TrialData):
	import numpy as np
	import matplotlib.pyplot as plt
	np.random.seed()

	RandomPath = TrialData["Default Paths"][np.random.randint(0,len(TrialData["Default Paths"]))]
	RandomMovement = reaching_movement(RandomPath)
	animate_plots(True,RandomMovement,TrialData)
def plot_all_on_same_axes(TrialData,Weighted=False,Statusbar = False,ReturnFig=False):
	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.patches import Ellipse
	import matplotlib.patches as patches
	import time

	t_end = TrialData["Movement Duration"]
	N = 1000
	t = np.linspace(0,t_end, N + 1)
	L1,L2 = TrialData["Limb Lengths"]
	MedianPlane = L1*(0.129/0.186)

	if TrialData["Reach Type"][:-6] == 'Left':
		DescriptiveTitle = "45$^\circ$ Reach Left\n"
	elif TrialData["Reach Type"][:-6] == 'Right':
		DescriptiveTitle = "45$^\circ$ Reach Right\n"
	elif TrialData["Reach Type"][:-6] == 'Sideways':
		DescriptiveTitle = "Side-to-side Reach\n"
	elif TrialData["Reach Type"][:-6] == 'Center':
		DescriptiveTitle = "Straight Forward (Center) Reach\n"
	elif TrialData["Reach Type"][:-6] == 'Fixed-target':
		DescriptiveTitle = "Fixed Target Reach\n"
	else:
		DescriptiveTitle = "Error in Title"

	fig, ((ax5,ax6),(ax1,ax2),(ax3,ax4)) = plt.subplots(3,2,figsize=(11,8))

	if t_end == 1:
		MovementDurationString = "Movement Duration : " + str(t_end) \
								+ " sec\n"
	else:
		MovementDurationString = "Movement Duration : " + str(t_end) \
								+ " secs\n"
	plt.suptitle(DescriptiveTitle + MovementDurationString,Fontsize=20,y=0.975)

	StartTime = time.time()
	for i in range(len(TrialData["Default Paths"])):
		Path = TrialData["Default Paths"][i]
		Movement = reaching_movement(Path)

		MuscleVelocities = Movement.return_muscle_velocities(TrialData,Weighted=Weighted)
		Vm_forward = MuscleVelocities
		Vm_reverse = -np.array(list(reversed(Vm_forward.T))).T

		A1_forward = Movement.A1
		A2_forward = Movement.A2
		A1_reverse = np.array(list(reversed(Movement.A1.T))).T
		A2_reverse = np.array(list(reversed(Movement.A2.T))).T

		#Forward Model

		Angle1_f, = ax5.plot(t.T,A1_forward.T,color = 'c')
		Angle2_f, = ax5.plot(t.T,A2_forward.T,color = 'r')
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
		movement_f, = ax1.plot(L1*np.cos(A1_forward[0,:])+L2*np.cos(A1_forward[0,:]+A2_forward[0,:]),\
							L1*np.sin(A1_forward[0,:])+L2*np.sin(A1_forward[0,:]+A2_forward[0,:]),\
							color='0.60')
		UpperArm_f = plt.Rectangle((0.02*np.sin(A1_forward[0,0]),\
										-0.02*np.cos(A1_forward[0,0])),\
										L1, 0.04,\
										angle=A1_forward[0,0]*180/np.pi,color='#4682b4')
		ax1.add_patch(UpperArm_f)
		Forearm_f = plt.Rectangle((L1*np.cos(A1_forward[0,0]) + 0.02*np.sin(A1_forward[0,0]+A2_forward[0,0]),\
										L1*np.sin(A1_forward[0,0]) -0.02*np.cos(A1_forward[0,0]+A2_forward[0,0])),\
										L2*(0.146/(0.146+0.108)), 0.04,\
										angle=(A1_forward[0,0]+A2_forward[0,0])*180/np.pi,color='#4682b4')
		ax1.add_patch(Forearm_f)
		Hand_f = plt.Rectangle((L1*np.cos(A1_forward[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.cos(A1_forward[0,0]+A2_forward[0,0]) \
									+ 0.02*np.sin(A1_forward[0,0]+A2_forward[0,0]),\
							L1*np.sin(A1_forward[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.sin(A1_forward[0,0]+A2_forward[0,0]) \
									-0.02*np.cos(A1_forward[0,0]+A2_forward[0,0])),\
										L2*(0.108/(0.146+0.108)), 0.04,\
										angle=(A1_forward[0,0]+A2_forward[0,0])*180/np.pi,color='#4682b4')
		ax1.add_patch(Hand_f)

		ax3.set_xlim(0,t_end)
		ax3.set_xticks([0,t_end])
		ax3.set_xticklabels(['Start','Finish'])
		if Weighted == True:
			ax3.set_ylabel('Afferent-Weighted $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
		else:
			ax3.set_ylabel('Normalized $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
		if np.shape(Vm_forward) == (len(Vm_forward),):
			NormalizedVmPlots_forward = ax3.plot(t,Vm_forward,c=TrialData["Ordered Muscle Colors List"][0])
			plt.figlegend(NormalizedVmPlots_forward,TrialData["Ordered Muscle List"],\
							loc='lower center',ncol=5,mode='expand')
			bounds = max([Vm_forward.max(),Vm_reverse.max()])
		else:
			NormalizedVmPlots_forward = [ax3.plot(t.T,Vm_forward[j].T) \
												for j in TrialData["Ordered Muscle Numbers"]]
			TotalColorsList = TrialData["Ordered Muscle Colors List"]\
								*len(TrialData["Default Paths"])
			[k.set_color(TotalColorsList[j]) for j,k in enumerate(ax3.lines)]

			bounds = max([	max([max(Vm_forward[i]) \
								for i in range(len(TrialData["Ordered Muscle List"]))]),\
							max([max(Vm_reverse[i]) \
								for i in range(len(TrialData["Ordered Muscle List"]))]) ] )
			plt.figlegend([el[0] for el in NormalizedVmPlots_forward],TrialData["Ordered Muscle List"],loc='lower center',ncol=5,mode='expand')

		ax3.spines['right'].set_visible(False)
		ax3.spines['top'].set_visible(False)
		ax3.set_ylim([-1.1*bounds,1.1*bounds])

		#Reverse Model

		Angle1_r, = ax6.plot(t.T,A1_reverse.T,color = 'c')
		Angle2_r, = ax6.plot(t.T,A2_reverse.T,color = 'r')
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
		movement_r, = ax2.plot(L1*np.cos(A1_reverse[0,:])+L2*np.cos(A1_reverse[0,:]+A2_reverse[0,:]),\
							L1*np.sin(A1_reverse[0,:])+L2*np.sin(A1_reverse[0,:]+A2_reverse[0,:]),\
							color='0.60')
		UpperArm_r = plt.Rectangle((0.02*np.sin(A1_reverse[0,0]),\
										-0.02*np.cos(A1_reverse[0,0])),\
										L1, 0.04,\
										angle=A1_reverse[0,0]*180/np.pi,color='#4682b4')
		ax2.add_patch(UpperArm_r)
		Forearm_r = plt.Rectangle((L1*np.cos(A1_reverse[0,0]) \
						+ 0.02*np.sin(A1_reverse[0,0]+A2_reverse[0,0]),\
						L1*np.sin(A1_reverse[0,0]) -0.02*np.cos(A1_reverse[0,0]+A2_reverse[0,0])),\
						L2*(0.146/(0.146+0.108)), 0.04,\
						angle=(A1_reverse[0,0]+A2_reverse[0,0])*180/np.pi,color='#4682b4')
		ax2.add_patch(Forearm_r)
		Hand_r = plt.Rectangle((L1*np.cos(A1_reverse[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.cos(A1_reverse[0,0]+A2_reverse[0,0]) \
									+ 0.02*np.sin(A1_reverse[0,0]+A2_reverse[0,0]),\
							L1*np.sin(A1_reverse[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.sin(A1_reverse[0,0]+A2_reverse[0,0]) \
									-0.02*np.cos(A1_reverse[0,0]+A2_reverse[0,0])),\
										L2*(0.108/(0.146+0.108)), 0.04,\
										angle=(A1_reverse[0,0]+A2_reverse[0,0])*180/np.pi,color='#4682b4')
		ax2.add_patch(Hand_r)


		#Might need to add t,WeightedNormalizedMuscleVelocity_Forward,TrialData["Ordered Muscle Numbers"],etc
		ax4.set_xlim(0,t_end)
		ax4.set_ylim(-12,12)
		ax4.set_xticks([0,t_end])
		ax4.set_xticklabels(['Start','Finish'])

		if np.shape(Vm_reverse) == (len(Vm_reverse),):
			NormalizedVmPlots_reverse = ax4.plot(t,Vm_reverse,c=TrialData["Ordered Muscle Colors List"][0])
		else:
			NormalizedVmPlots_reverse = [ax4.plot(t.T,Vm_reverse[j].T) \
												for j in TrialData["Ordered Muscle Numbers"]]
			TotalColorsList = TrialData["Ordered Muscle Colors List"]\
								*len(TrialData["Default Paths"])
			[k.set_color(TotalColorsList[j]) for j,k in enumerate(ax4.lines)]
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

		if Statusbar == True:
			statusbar(i,len(TrialData["Default Paths"]),StartTime=StartTime,\
							Title = TrialData["Reach Type"])
	print('\n')
	if ReturnFig == True: return(fig)
def plot_all_trajectories(TrialData,Statusbar = False,ReturnFig=False):
	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.patches import Ellipse
	import matplotlib.patches as patches
	import time

	t_end = TrialData["Movement Duration"]
	N = 1000
	t = np.linspace(0,t_end, N + 1)
	L1,L2 = TrialData["Limb Lengths"]
	MedianPlane = L1*(0.129/0.186)

	if TrialData["Reach Type"][:-6] == 'Left':
		DescriptiveTitle = "45$^\circ$ Reach Left\n"
	elif TrialData["Reach Type"][:-6] == 'Right':
		DescriptiveTitle = "45$^\circ$ Reach Right\n"
	elif TrialData["Reach Type"][:-6] == 'Sideways':
		DescriptiveTitle = "Side-to-side Reach\n"
	elif TrialData["Reach Type"][:-6] == 'Straight':
		DescriptiveTitle = "Straight Forward (Center) Reach\n"
	elif TrialData["Reach Type"][:-6] == 'Fixed-target':
		DescriptiveTitle = "Fixed-Target Reach\n"
	else:
		DescriptiveTitle = "Error in Title"

	fig, (ax1,ax2) = plt.subplots(2,1,figsize=(11,8))

	if t_end == 1:
		MovementDurationString = "Movement Duration : " + str(t_end) \
								+ " sec\n"
	else:
		MovementDurationString = "Movement Duration : " + str(t_end) \
								+ " secs\n"
	plt.suptitle(DescriptiveTitle + MovementDurationString,Fontsize=20,y=0.975)

	StartTime = time.time()
	for i in range(len(TrialData["Default Paths"])):
		Path = TrialData["Default Paths"][i]
		Movement = reaching_movement(Path)
		X,_,_ = Movement.return_X_values(TrialData)
		Movement.inverse_kinematics(X,TrialData)
		A1 = Movement.A1
		A2 = Movement.A2

		Angle1, = ax1.plot(t.T,A1.T,color = 'c')
		Angle2, = ax1.plot(t.T,A2.T,color = 'r')
		ax1.set_xlim(0,t_end)
		ax1.set_xticks([0,t_end])
		ax1.set_xticklabels(['Start','Finish'])
		ax1.set_ylim(-np.pi/2,np.pi)
		ax1.set_yticks([-np.pi/2,0,np.pi/2,np.pi])
		ax1.set_yticklabels([r'$\frac{-\pi}{2}$','0',r'$\frac{\pi}{2}$','$\pi$'],fontsize=12)
		ax1.spines['right'].set_visible(False)
		ax1.spines['top'].set_visible(False)
		# ax1.set_ylabel('Joint Angles\n(in Radians)')
		ax1.legend(["Shoulder\nAngle","Elbow\nAngle"],loc='best')

		ax2.get_xaxis().set_ticks([])
		ax2.get_yaxis().set_ticks([])
		ax2.set_frame_on(True)
		RightShoulder = plt.Circle((0,0),radius=0.05,Color='#4682b4')
		ax2.add_patch(RightShoulder)
		LeftShoulder = plt.Circle((-MedianPlane*2,0),radius=0.05,Color='#4682b4')
		ax2.add_patch(LeftShoulder)
		Torso = plt.Rectangle((-MedianPlane*2,-0.05),MedianPlane*2,0.1,Color='#4682b4')
		ax2.add_patch(Torso)
		Head = Ellipse((-MedianPlane,0),0.2,0.225,FaceColor='w',EdgeColor='#4682b4',Linewidth=3)
		ax2.add_patch(Head)
		JointCoordinates = \
			np.concatenate([[np.cumsum([0,L1*np.cos(A1[0,0]),L2*(0.146/(0.146+0.108))*np.cos(A1[0,0]+A2[0,0]),L2*(0.108/(0.146+0.108))*np.cos(A1[0,0]+A2[0,0])])],\
			[np.cumsum([0,L1*np.sin(A1[0,0]),L2*(0.146/(0.146+0.108))*np.sin(A1[0,0]+A2[0,0]),L2*(0.108/(0.146+0.108))*np.sin(A1[0,0]+A2[0,0])])]],\
			axis=0)
		Elbow = plt.Circle((L1*np.cos(A1[0,0]),L1*np.sin(A1[0,0])),radius=0.03,color='#4682b4')
		ax2.add_patch(Elbow)
		Endpoint = plt.Circle((L1*np.cos(A1[0,0])+L2*np.cos(A1[0,0]+A2[0,0]),\
							L1*np.sin(A1[0,0])+L2*np.sin(A1[0,0]+A2[0,0])),\
							radius = 0.02,color='#4682b4')
		ax2.add_patch(Endpoint)
		Wrist = plt.Circle((L1*np.cos(A1[0,0])+L2*(0.146/(0.146+0.108))*np.cos(A1[0,0]+A2[0,0]),\
							L1*np.sin(A1[0,0])+L2*(0.146/(0.146+0.108))*np.sin(A1[0,0]+A2[0,0])),\
							radius = 0.03,color='#4682b4')
		ax2.add_patch(Wrist)
		StickFigure, = ax2.plot(JointCoordinates[0,:],JointCoordinates[1,:],'ko-',LineWidth=2,MarkerFaceColor='k')
		movement, = ax2.plot(L1*np.cos(A1[0,:])+L2*np.cos(A1[0,:]+A2[0,:]),\
							L1*np.sin(A1[0,:])+L2*np.sin(A1[0,:]+A2[0,:]),\
							color='0.60')
		UpperArm = plt.Rectangle((0.02*np.sin(A1[0,0]),\
										-0.02*np.cos(A1[0,0])),\
										L1, 0.04,\
										angle=A1[0,0]*180/np.pi,color='#4682b4')
		ax2.add_patch(UpperArm)
		Forearm = plt.Rectangle((L1*np.cos(A1[0,0]) + 0.02*np.sin(A1[0,0]+A2[0,0]),\
										L1*np.sin(A1[0,0]) -0.02*np.cos(A1[0,0]+A2[0,0])),\
										L2*(0.146/(0.146+0.108)), 0.04,\
										angle=(A1[0,0]+A2[0,0])*180/np.pi,color='#4682b4')
		ax2.add_patch(Forearm)
		Hand = plt.Rectangle((L1*np.cos(A1[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.cos(A1[0,0]+A2[0,0]) \
									+ 0.02*np.sin(A1[0,0]+A2[0,0]),\
							L1*np.sin(A1[0,0])\
								+ L2*(0.146/(0.146+0.108))*np.sin(A1[0,0]+A2[0,0]) \
									-0.02*np.cos(A1[0,0]+A2[0,0])),\
										L2*(0.108/(0.146+0.108)), 0.04,\
										angle=(A1[0,0]+A2[0,0])*180/np.pi,color='#4682b4')
		ax2.add_patch(Hand)
		ax2.set_aspect('equal')
		if Statusbar == True:
			statusbar(i,len(TrialData["Default Paths"]),StartTime=StartTime,\
							Title = TrialData["Reach Type"])
	print('\n')
	if ReturnFig == True: return(fig)
def plot_individual_muscles(TrialData,Weighted=False,Statusbar=False,\
								SameScale=False,Scale=[],ReturnFig=False):
	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.patches import Ellipse
	import matplotlib.patches as patches
	import time

	TotalColorsList = TrialData["Ordered Muscle Colors List"]*len(TrialData["Default Paths"])

	NumMuscles = len(TrialData["All Muscle Settings"])
	NumRows = int(np.ceil(NumMuscles/5))
	if NumMuscles < 5:
		NumColumns = NumMuscles
	else:
		NumColumns = 5

	ColumnNumber = [el%5 for el in np.arange(0,NumMuscles,1)]
	RowNumber = [int(el/5) for el in np.arange(0,NumMuscles,1)]

	fig, axes = plt.subplots(NumRows,NumColumns,figsize=(3*NumColumns,2*NumRows + 2))

	t_end = TrialData["Movement Duration"]
	N = 1000
	t = np.linspace(0,t_end, N + 1)
	L1,L2 = TrialData["Limb Lengths"]
	MedianPlane = L1*(0.129/0.186)
	bounds = [0]*NumMuscles

	if TrialData["Reach Type"][:-6] == 'Left':
		DescriptiveTitle = "45$^\circ$ Reach Left\n"
	elif TrialData["Reach Type"][:-6] == 'Right':
		DescriptiveTitle = "45$^\circ$ Reach Right\n"
	elif TrialData["Reach Type"][:-6] == 'Sideways':
		DescriptiveTitle = "Side-to-side Reach\n"
	elif TrialData["Reach Type"][:-6] == 'Center':
		DescriptiveTitle = "Straight Forward (Center) Reach\n"
	elif TrialData["Reach Type"][:-6] == 'Fixed-target':
		DescriptiveTitle = "Fixed Target Reach\n"
	else:
		DescriptiveTitle = "Error in Title\n"

	if t_end == 1:
		MovementDurationString = "Movement Duration : " + str(t_end) \
								+ " sec\n"
	else:
		MovementDurationString = "Movement Duration : " + str(t_end) \
								+ " secs\n"
	plt.suptitle(DescriptiveTitle + MovementDurationString,Fontsize=20,y=0.975)
	if NumMuscles == 1:
		axes.set_xlim(0,t_end)
		axes.set_xticks([0,t_end])
		axes.set_title(list(TrialData["All Muscle Settings"].keys())[0],Color = TrialData["Ordered Muscle Colors List"][0])
		axes.spines['right'].set_visible(False)
		axes.spines['top'].set_visible(False)
		if Weighted == True:
			axes.set_ylabel('Afferent-Weighted $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
		else:
			axes.set_ylabel('Normalized $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
		# bounds = max([Vm_forward.max(),Vm_reverse.max()])
		StartTime = time.time()
		for i in range(len(TrialData["Default Paths"])):
			Path = TrialData["Default Paths"][i]
			Movement = reaching_movement(Path)

			MuscleVelocities = Movement.return_muscle_velocities(TrialData,Weighted=Weighted)
			Vm_forward = MuscleVelocities
			Vm_reverse = -np.array(list(reversed(Vm_forward.T))).T
			axes.plot(t.T,Vm_forward.T,c=TrialData["Ordered Muscle Colors List"][0])
			bounds = max([ max(Vm_forward), bounds ])
			if Statusbar == True:
				statusbar(i,len(TrialData["Default Paths"]),StartTime=StartTime,\
								Title = TrialData["Reach Type"])
		axes.set_ylim([-1.1*bounds,1.1*bounds])
	else:
		for j in range(NumMuscles):
			axes[RowNumber[j],ColumnNumber[j]].set_xlim(0,t_end)
			axes[RowNumber[j],ColumnNumber[j]].set_xticks([0,t_end])
			if RowNumber[j] == RowNumber[-1] and ColumnNumber[j]==0:
				axes[RowNumber[j],ColumnNumber[j]].set_xticklabels(['Start','Finish'])
				if Weighted == True:
					axes[RowNumber[j],ColumnNumber[j]].set_ylabel('Afferent-Weighted $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
				else:
					axes[RowNumber[j],ColumnNumber[j]].set_ylabel('Normalized $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
			else:
				axes[RowNumber[j],ColumnNumber[j]].set_xticklabels(['',''])
			axes[RowNumber[j],ColumnNumber[j]].set_title(TrialData["Ordered Muscle List"][j],\
									Color = TrialData["Ordered Muscle Colors List"][j])
			# if ColumnNumber[j] == 0:
			# 	if Weighted == True:
			# 		axes[RowNumber[j],ColumnNumber[j]].set_ylabel('Afferent-Weighted $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
			# 	else:
			# 		axes[RowNumber[j],ColumnNumber[j]].set_ylabel('Normalized $\hat{v}_m$\nConcentric $\longleftrightarrow$ Eccentric')
		if NumMuscles%5!=0:
			[fig.delaxes(axes[RowNumber[-1],el]) for el in range(ColumnNumber[-1]+1,5)]
		StartTime = time.time()
		for i in range(len(TrialData["Default Paths"])):
			Path = TrialData["Default Paths"][i]
			Movement = reaching_movement(Path)

			MuscleVelocities = Movement.return_muscle_velocities(TrialData,Weighted=Weighted)
			Vm_forward = MuscleVelocities
			Vm_reverse = -np.array(list(reversed(Vm_forward.T))).T

			for j in range(NumMuscles):
				MuscleNumber = TrialData["Ordered Muscle Numbers"][j]
				axes[RowNumber[j],ColumnNumber[j]].plot(t.T,Vm_forward[MuscleNumber].T,\
						c=TrialData["Ordered Muscle Colors List"][j])
				bounds[j] = \
					max([ max(abs(Vm_forward[MuscleNumber])), bounds[j] ])
			if Statusbar == True:
				statusbar(i,len(TrialData["Default Paths"]),StartTime=StartTime,\
								Title = TrialData["Reach Type"])
		print('\n')
		if SameScale == False:
			for j in range(NumMuscles):
				axes[RowNumber[j],ColumnNumber[j]].spines['right'].set_visible(False)
				axes[RowNumber[j],ColumnNumber[j]].spines['top'].set_visible(False)
				axes[RowNumber[j],ColumnNumber[j]].set_ylim([-1.1*bounds[j],1.1*bounds[j]])
		else:
			if Scale == []:
				MaxBounds = max(bounds)
				for j in range(NumMuscles):
					axes[RowNumber[j],ColumnNumber[j]].spines['right'].set_visible(False)
					axes[RowNumber[j],ColumnNumber[j]].spines['top'].set_visible(False)
					axes[RowNumber[j],ColumnNumber[j]].set_ylim([-1.1*MaxBounds,1.1*MaxBounds])
			else:
				assert type(Scale)==list and len(Scale)==2, "Scale must be list of length 2"
				for j in range(NumMuscles):
					axes[RowNumber[j],ColumnNumber[j]].spines['right'].set_visible(False)
					axes[RowNumber[j],ColumnNumber[j]].spines['top'].set_visible(False)
					axes[RowNumber[j],ColumnNumber[j]].set_ylim(Scale)
	if ReturnFig == True: return(fig)
def save_trial_data_prompt(TrialData):
	import numpy as np
	import pickle
	import os.path
	ValidResponse_1 = False
	while ValidResponse_1 == False:
		Response = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nSave Trial Data:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n (1) - Yes (Default)\n (2) - No\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nResponse: ")
		if Response not in ['1','2','']:
			print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInvalid Response! Please try again.')
			ValidResponse_1 = False
		elif Response == '' or Response == '1':
			ValidResponse_2 = False
			while ValidResponse_2 == False:
				FileName = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nFileName: ")
				if os.path.exists(FileName.capitalize() + ".pkl") == True:
					ValidResponse_3 = False
					while ValidResponse_3 == False:
						Replace = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" + FileName +".pkl already exists. Replace file?\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n (1) - Yes (Default)\n (2) - No\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nResponse: ")
						if Replace not in ['1','2','']:
							print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInvalid Response! Please try again.')
							ValidResponse_3 = False
						elif Replace == "" or Replace == "1":
							print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
							pickle.dump(TrialData,open(FileName.capitalize()+'.pkl','wb'),\
								pickle.HIGHEST_PROTOCOL)
							ValidResponse_2 = True
							ValidResponse_3 = True
						elif Replace == "2":
							print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
							ValidResponse_2 = False
							ValidResponse_3 = True
				else:
					print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
					pickle.dump(TrialData,open(FileName.capitalize()+'.pkl','wb'),\
						pickle.HIGHEST_PROTOCOL)
					ValidResponse_2 = True
			ValidResponse_1 = True
		else:
			print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
			ValidResponse_1 = True
def load_trial_data_prompt():
	import numpy as np
	import pickle
	import os.path
	ValidResponse_1 = False
	while ValidResponse_1 == False:
		Response = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nLoad Trial Data:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n (1) - Yes (Default)\n (2) - No\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nResponse: ")
		if Response not in ['1','2','']:
			print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nInvalid Response! Please try again.')
			ValidResponse_1 = False
		elif Response == '' or Response == '1':
			ValidResponse_2 = False
			while ValidResponse_2 == False:
				FileName = input("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nFileName: ")
				if FileName.capitalize() == 'Cancel' or FileName.capitalize() == 'Exit':
					ValidResponse_2 = True
				elif os.path.exists(FileName.capitalize() + ".pkl") == False:
					print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n' + FileName.capitalize() + '.pkl not found. Please try again.\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
					ValidResponse_2 = False
				else:
					print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
					TrialData = pickle.load(open(FileName+'.pkl','rb'))
					ValidResponse_2 = True
			ValidResponse_1 = True
			return(TrialData)
		else:
			print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
			ValidResponse_1 = True
			pass
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
def save_figures(TrialData,figs):
	import os.path
	from matplotlib.backends.backend_pdf import PdfPages
	i = 1
	FileName = TrialData["Reach Type"].replace(' ','_').capitalize() + "_" \
					+ "{:0>2d}".format(i) +".pdf"
	if os.path.exists(FileName) == True:
		while os.path.exists(FileName) == True:
			i += 1
			FileName = TrialData["Reach Type"].replace(' ','_').capitalize()\
			 				+ "_" + "{:0>2d}".format(i) +".pdf"
	PDFFile = PdfPages(FileName)
	if len(figs)==1:
		PDFFile.savefig(figs)
	else:
		[PDFFile.savefig(fig) for fig in figs]
	PDFFile.close()

# NumberOfTrials should be greater than 50 if using statusbar())
NumberOfTrials = 100
if NumberOfTrials<50:
	Statusbar_bool = False
else:
	Statusbar_bool = True
TrialData = create_trial_data(NumberOfTrials=NumberOfTrials)
# animate_random_trajectory(TrialData)
# plot_all_on_same_axes(TrialData,Statusbar=Statusbar_bool)
fig1 = plot_all_trajectories(TrialData,Statusbar=Statusbar_bool,ReturnFig=True)
fig2 = plot_individual_muscles(TrialData,Statusbar=Statusbar_bool,\
								SameScale=True,Scale=[-1.2,1.2],ReturnFig=True)
save_figures(TrialData,[fig1,fig2])
save_trial_data_prompt(TrialData)
# load_trial_data_prompt()
