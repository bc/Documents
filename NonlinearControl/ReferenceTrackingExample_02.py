import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import sympy as sy
from sympy.utilities import lambdify
import time
from collections import namedtuple

def return_primary_source(Settings):
	import numpy as np
	assert Settings["Primary Source"]!=None, "No sources were found for this setting."
	TotalSources = Settings["Sources"]
	PrimarySource = Settings["Primary Source"]
	assert PrimarySource in [settings.Source for settings in TotalSources], "Error! Primary Source is not referenced."
	return(TotalSources[np.where([settings.Source == PrimarySource for settings in TotalSources])[0][0]])
def return_muscle_settings(PreselectedMuscles=None):
	"""
	Notes:
	Coefficients from observation, Ramsay; 2009, FVC, Holtzbaur, Pigeon, Kuechle, or Banks. Optimal Muscle Length given in mm. Optimal tendon/muscle lengths and PCSA were taken from Garner and Pandy (2003)
	"""
	import sympy as sp
	from sympy.utilities import lambdify
	import numpy as np
	from numpy import pi

	global θ_SFE,θ_EFE,θ_PS
	θ_SFE,θ_EFE,θ_PS = sp.symbols('θ_SFE'),sp.symbols('θ_EFE'),sp.symbols('θ_PS')

	# Coefficients from observation, Ramsay; 2009, Pigeon, FVC, Holtzbaur, or Banks.

	MA_Settings = namedtuple("MA_Settings",["Values","Source","Units","Equation_Number","Threshold","DOF"])
	Spindle_Settings = namedtuple("Spindle_Settings",["ActualNumber",'CorrectedNumber','RelativeAbundance',"Source"])
	Input_Source = namedtuple("Source_Settings"["Values","Source","Units"])
	def Pigeon_coeff_conversion(Coefficients):
		"""
		Takes in Coefficient values from Pigeon (1996) -- which take in angles in degrees -- and coverts them into the properly scaled coefficients for radians, additionally scaled by the magnitude listed in the paper.

		Note that the coefficients listed in Pigeon (1996) are given in decending order (i.e., c₅,c₄,c₃,c₂,c₁,c₀). However to maintain continuity with the equations given in Ramsay; 2009 (2009), we list coefficients in order of increasing power (i.e., c₀,c₁,c₂,c₃,c₄,c₅).
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
		'Notes' : [\
						'This is the *clavicular* portion of the pectoralis major.',\
						'Banks and Garner & Pandy are parameter values for the entire muscle. Pigeon and Holzbaur have the values for the clavicular portion only.',\
						'Holzbaur parameters for shoulder are for frontal plane (ABD/ADD) only! This explains the relatively small values for some measures.'\
					],\
		'Shoulder MA' : {	"Primary Source" : "Pigeon; 1996",\
		 					"Sources" : \
								[\
									MA_Settings([50.80,0,0,0,0,0], 'Pigeon; 1996', "mm", None, None, "Shoulder"),\
									MA_Settings([2,0,0,0,0,0], 'Holzbaur; 2005', "mm", None, None, "Shoulder")\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Est", \
		 				"Sources" : \
							[\
								MA_Settings(0, "m", None, None, 'Elbow', "Est")\
							]}, \
		'Spindle' : Spindle_Settings(450,389.7,1.2,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006", \
					"Sources" : \
						[\
							Input_Source(295.6, "Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(14.4, 'Holzbaur; 2005', 'cm'),\
											Input_Source(150, 'Est', 'mm')\
										]}, \
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(0.3, 'Holzbaur; 2005', 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(17, 'Holzbaur; 2005', 'degrees')\
									]}, \
		'PCSA' : {	"Primary Source" : "Garner & Pandy; 2003", \
					'Sources' : \
						[\
							Input_Source(36.20,'Garner & Pandy; 2003','sq cm'),\
							Input_Source(2.6,'Holzbaur; 2005','sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Garner & Pandy; 2003", \
										'Sources': \
											[\
												Input_Source(1175.01,'Garner & Pandy; 2003','N'),\
												Input_Source(364.4,'Holzbaur; 2005','N')\
											]}\
		}

	DELTa_Settings = {\
		'Notes' : [\
					"SFE MA is listed as 33.02 mm in Pigeon and estimated as 19 mm. Using Pigeon Coefficients convention, Kuechle (1997) has the DELTp MA for [-140,90] as Pigeon_coeff_conversion([ 13.4293189,  2.0316226, -0.2339031,  2.7807828,  0.,  0.]). This will yield a piecewise function that creates jumps with the new velocity formulation. Instead, we are going to try Pigeon_coeff_conversion([ 12.7928795,  2.0480346,  0.8917734,  3.2207214, -2.3928223,  0.]) so that the function is within range and continuous during the ROM. Threshold (pi/2) has been removed for this new MA function.",\
					"Garner & Pandy have much larger PCSA and Peak Force Values but only consider the entire Deltoid.",\
					"Holzbaur parameters for shoulder are for frontal plane (ABD/ADD) only! This explains the relatively small values for some measures.",\
					"Banks only had mass and spindle settings for the entire deltoid. As a result the parameters are divided by 3 as an estimate of the individual muscles. Will need to do sensitivity analysis as a result."\
					], \
		'Shoulder MA' : {	"Primary Source" : "Kuechle; 1997",\
		 					"Sources" : \
								[\
									MA_Settings(Pigeon_coeff_conversion([12.7928795,  2.0480346,  0.8917734,  3.2207214, -2.3928223,  0.]), 'Kuechle; 1997', "mm", None, None, 'Shoulder'),\
									MA_Settings([33.02,0,0,0,0,0], 'Pigeon; 1996', "mm", None, None, 'Shoulder'),\
									MA_Settings([1.9,0,0,0,0,0], 'Holzbaur; 2005', "cm", None, None, 'Shoulder'),\
									MA_Settings(19, "Est", "mm", None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Est", \
		 				"Sources" : \
							[\
								MA_Settings(0, "Est", "m", None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(182/3,426.3/3,0.43,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006", \
					"Sources" : \
						[\
							Input_Source(355.7/3, "Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(9.8, 'Holzbaur; 2005', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(9.3, 'Holzbaur; 2005', 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(22, 'Holzbaur; 2005', 'degrees')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(82.98,'Garner & Pandy; 2003','sq cm'),\
							Input_Source(8.2,'Holzbaur; 2005','sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(2044.65,'Garner & Pandy; 2003','N'),\
												Input_Source(1142.6,'Holzbaur; 2005','N')\
											]}\
		}

	CB_Settings = {\
		'Notes' : [\
						"Holzbaur parameters for shoulder are for frontal plane (ABD/ADD) only! This explains the relatively small values for some measures. MA is negative as a result.",\
						"Garner & Pandy values for muscle length, PCSA, and peak force are very different from those reported in Wood (1989), Veeger (1991), Bassett (1990), Chen (1988), Keating (1993), Veeger (1997), An (1981), and Cutts (1991)." \
					],\
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(20, "Est", "mm", None, None, "Shoulder"),\
									MA_Settings([-20,0,0,0,0,0], "Holzbaur; 2005", "mm", None, None, "Shoulder")\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Est", \
		 				"Sources" : \
							[\
								MA_Settings(0, "Est", "m", None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(123,147.3,0.83,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : 'Banks; 2006',\
					'Sources' : \
						[\
							Input_Source(39.8, "Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(9.3, 'Holzbaur; 2005', 'cm'),\
											Input_Source(17.60, 'Garner & Pandy; 2003', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : 'Holzbaur; 2005', \
									"Sources" : \
										[\
											Input_Source(9.7, 'Holzbaur; 2005', 'cm'),\
											Input_Source(4.23, 'Garner & Pandy; 2003', 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(27, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(1.7,"Holzbaur; 2005","sq cm"),\
							Input_Source(4.55,"Garner & Pandy; 2003","sq cm")\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(242.5, "Holzbaur; 2005", "N"),\
												Input_Source(150.02, "Garner & Pandy; 2003", "N")\
											]}\
		}

	DELTp_Settings = {\
		'Notes' : [\
						"DELTp SFE MA is listed as -78.74 mm in Pigeon. Using Pigeon Coefficients convention, Kuechle (1997) has the DELTp MA for [-140,90] as Pigeon_coeff_conversion([ 22.8547177,  3.9721238, -3.3900829, -3.6146546,  0.,  0.]). This will yield a piecewise function that creates jumps with the new velocity formulation. Instead, we are going to try Pigeon_coeff_conversion([-23.8165173, -4.486164 ,  5.8655808,  6.5003255, -8.2736695,2.0812998]) so that the function is within range and continuous during the ROM. Threshold (pi/2) has been removed for this new MA function.",\
						"Holzbaur parameters for shoulder are for frontal plane (ABD/ADD) only! This explains the relatively small values for some measures. MA is negative as a result.",\
						"Garner & Pandy values for muscle length, PCSA, and peak force are very different from those reported in Wood (1989), Veeger (1991), Bassett (1990), Chen (1988), Keating (1993), Veeger (1997), An (1981), and Cutts (1991). Also, they do no distinguish between ant, mid, post.",\
						"Holzbaur parameters for shoulder are for frontal plane (ABD/ADD) only! This explains the relatively small values for some measures.",\
						"Banks only had mass and spindle settings for the entire deltoid. As a result the parameters are divided by 3 as an estimate of the individual muscles. Will need to do sensitivity analysis as a result."\
					],\
		'Shoulder MA' : {	"Primary Source" : "Kuechle; 1997",\
		 					"Sources" : \
								[\
									MA_Settings(Pigeon_coeff_conversion([-23.8165173, -4.486164 ,  5.8655808,  6.5003255, -8.2736695,2.0812998]), 'Kuechle; 1997', "mm", None, None, "Shoulder"),\
									MA_Settings([-78.74,0,0,0,0,0], 'Pigeon; 1996', "mm", None, None, "Shoulder"),\
									MA_Settings([-8,0,0,0,0,0], 'Holzbaur; 2005', "mm", None, None, "Shoulder")
								]}, \
		'Elbow MA' : {	"Primary Source" : "Est", \
		 				"Sources" : \
							[\
								MA_Settings(0, "Est", "m", None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(182/3,426.3/3,0.43,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006", \
					"Sources" : \
						[\
							Input_Source(355.7/3, "Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(13.7, 'Holzbaur; 2005', 'cm'),\
											Input_Source(12.8, 'Garner & Pandy; 2003', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(3.8, "Holzbaur; 2005", "cm"),\
											Input_Source(5.38, 'Garner & Pandy; 2003', 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(18, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(1.9, "Holzbaur; 2005", 'sq cm'),\
							Input_Source(81.98,"Garner & Pandy; 2003","sq cm")\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(259.9,"Holzbaur; 2005","N"),\
												Input_Source(2044.65,"Garner & Pandy; 2003","N")\
											]}\
		}}

	BIC_Settings = {\
		'Notes' : [\
					"BIC EFE MA for Ramsay; 2009 has R² = 0.985 whereas Pigeon has R² = 0.9918. Pigeon, however, only takes elbow angle into account, whereas Ramsay; 2009 takes in variable PS angles. It appears that because Pigeon uses an average of fully pronated and fully supinated MAs, the BIC moment arm is similar but subject to variation as the level of PS is changed. (NOTE: BIC becomes slightly negative when q2 > 3.021. If trajectory has elbow angles exceding this value, enter a threshold of 3.021 into the model.)",\
					"Note: Only using the long head for optimal length, see Holzbaur (2005) for additional head parameters. Adding when logical."
					],\
		'Shoulder MA' : {	"Primary Source" : "Pigeon; 1996",\
		 					"Sources" : \
								[\
									MA_Settings([29.21,0,0,0,0,0], 'Pigeon; 1996', "mm", None, None, "Shoulder"),\
									MA_Settings(15, 'Est', "mm", None, None, "Shoulder")\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Ramsay; 2009", \
	 					"Sources" : \
							[\
								MA_Settings([8.4533,36.6147,2.4777,-19.432,2.0571,0,13.6502,0,0,-5.6172,0,-2.0854,0,0,0,0], 'Ramsay; 2009', "mm", 2, 3.021, 'Elbow'),\
								MA_Settings(Pigeon_coeff_conversion([14.660,4.5322,1.8047,-2.9883,0,0]), 'Pigeon; 1996', "mm", None, 2.9326, 'Elbow'),\
								MA_Settings([36,0,0,0,0,0], 'Holzbaur; 2005', "mm", None, None, "Elbow")
							]}, \
		'Spindle' : Spindle_Settings(320,292.6,1.1,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006", \
					"Sources" : \
						[\
							Input_Source(163.8,"Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(11.6, 'Holzbaur; 2005', 'cm'),\
											Input_Source(14.22, "Garner & Pandy; 2003", "cm")
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(27.2, "Holzbaur; 2005", "cm"),\
											Input_Source(22.98, "Garner & Pandy; 2003", "cm")
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(0, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source((4.5+3.1), "Holzbaur; 2005", "sq cm"),\
							Input_Source(25.90, "Garner & Pandy; 2003", "sq cm")
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source((624.3+435.6), "Holzbaur; 2005", "N"),\
												Input_Source(849.29, "Garner & Pandy; 2003", "N")
											]}\
		}}

	TRI_Settings = {\
		'Notes' : [\
					"TRI EFE MA for Ramsay; 2009 has R² = 0.997 whereas Pigeon has R² = 0.9904. Pigeon appears to really fail when the elbow angle is greater than 140°. For this reason, Ramsay; 2009 should be used. However the approach of fixing the MA for values greater than 140° can be adopted for completeness. Coefficients and equation number/type are listed below to test either implementation.",\
					"Note: Only using the long head for optimal length, see Holzbaur (2005) for additional head parameters.",\
					"Banks had the parameters for each head of the triceps, values were added.",\
					"Holzbaur settings only utilizes the long head of the TRI."\
					],
		'Shoulder MA' : {	"Primary Source" : "Pigeon; 1996",\
		 					"Sources" : \
								[\
									MA_Settings([-25.40,0,0,0,0,0], 'Pigeon; 1996', "mm", None, None, "Shoulder"),\
									MA_Settings(-15, 'Est', "mm", None, None, "Shoulder")\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Ramsay; 2009", \
	 					"Sources" : \
							[\
								MA_Settings([-24.5454,-8.8691,9.3509,-1.7518,0], 'Ramsay; 2009', 'mm', 1, None, 'Elbow'),\
								MA_Settings(Pigeon_coeff_conversion([-23.287,-3.0284,12.886,-19.092,13.277,-3.5171]), 'Pigeon; 1996', 'mm', None, None, 'Elbow'),\
								MA_Settings([-21,0,0,0,0,0], 'Holzbaur; 2005', 'mm', None, None, 'Elbow')
							]}, \
		'Spindle' : Spindle_Settings((200+222+98),(223.7+269.6+221.8),(0.89+0.82+0.44)/3,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006",\
					"Sources" : \
						[\
							Input_Source((94.2+138.4+92.5), "Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005",\
		 							"Sources" : \
										[\
											Input_Source(13.4, 'Holzbaur; 2005', 'cm'),\
											Input_Source(8.77, 'Garner & Pandy; 2003', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(14.3, "Holzbaur; 2005", 'cm'),\
											Input_Source(19.05, 'Garner & Pandy; 2003', 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(12, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source((5.7+4.5+4.5),"Holzbaur; 2005",'sq cm'),\
							Input_Source(76.30, 'Garner & Pandy; 2003', 'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source((798.5+624.3+624.3),"Holzbaur; 2005",'N'),\
												Input_Source(2332.92, 'Garner & Pandy; 2003', 'N')\
											]}\
		}}

	BRA_Settings = {\
		"Notes" : [\
					"BRA (Brachialis) EFE MA for Ramsay; 2009 has R² = 0.990 whereas Pigeon has R² = 0.9988. Curve appears to be a better fit, as it experiences its smallest MA when Elbow angle = 0. Coefficients and equation number/type are listed below to test either implementation."\
					],
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Ramsay; 2009", \
	 					"Sources" : \
							[\
								MA_Settings([16.1991,-16.1463,24.5512,-6.3335,0], 'Ramsay; 2009', 'mm', 1, None, 'Elbow'),\
								MA_Settings(Pigeon_coeff_conversion([5.5492,2.3080,2.3425,-2.0530,0,0]), 'Pigeon; 1996', 'mm', None, None, 'Elbow'),\
								MA_Settings([18,0,0,0,0,0], 'Holzbaur; 2005', 'mm', None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(256,272.1,0.94,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006", \
					"Sources" : \
						[\
							Input_Source(141, "Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(8.6, 'Holzbaur; 2005', 'cm'),\
											Input_Source(10.28, 'Holzbaur; 2005', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(5.4, "Holzbaur; 2005", 'cm'),\
											Input_Source(1.75, "Holzbaur; 2005", 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(0, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(7.1,"Holzbaur; 2005",'sq cm'),\
							Input_Source(25.88,"Garner & Pandy; 2003",'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(987.3,"Holzbaur; 2005","N"),\
												Input_Source(583.76,"Garner & Pandy; 2003","N")\
											]}\
		}}

	BRD_Settings = {\
		"Notes" : [\
					"BRD (Brachioradialis) for Ramsay; 2009 has R² = 0.988 whereas Pigeon has R² = 0.9989. Pigeon, however, only takes elbow angle into account, whereas Ramsay; 2009 takes in variable PS angles. Coefficients and equation number/type are listed below to test either implementation."\
					],
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Ramsay; 2009", \
	 					"Sources" : \
							[\
								MA_Settings(	[15.2564,-11.8355,2.8129,-5.7781,44.8143,0,2.9032,0,0,-13.4956,0,-0.3940,0,0,0,0], 'Ramsay; 2009', 'mm', 2, None, 'Elbow'), \
								MA_Settings(	Pigeon_coeff_conversion([19.490,1.6681,10.084,-6.5171,0,0]), 'Pigeon; 1996', 'mm', None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(70,190.2,0.37,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006",\
					"Sources" : \
						[\
							Input_Source(64.7,"Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(17.3, 'Holzbaur; 2005', 'cm'),\
											Input_Source(27.03, 'Garner & Pandy; 2003', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(13.3, "Holzbaur; 2005", 'cm'),\
											Input_Source(6.04, "Garner & Pandy; 2003", 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(0, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(1.9,"Holzbaur; 2005",'sq cm'),\
							Input_Source(3.08,"Garner & Pandy; 2003",'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(261.3,"Holzbaur; 2005",'N'),\
												Input_Source(101.56,"Garner & Pandy; 2003",'N')\
											]}\
		}}

	PRO_Settings = {\
		'Notes' : [\
					""\
					]
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Ramsay; 2009", \
	 					"Sources" : \
							[\
								MA_Settings(	[11.0405,-1.0079,0.3933,-10.4824,-12.1639,-0.4369,36.9174,3.5232,-10.4223,21.2604,-37.2444,10.2666,-11.0060,14.5974,-3.9919,1.7526,-2.0089,0.5460], 'Ramsay; 2009', 'mm', 3, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(187.6,185.5,1.3,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006",\
					"Sources" : \
						[\
							Input_Source(38.8, "Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(4.9, 'Holzbaur; 2005', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(9.8, "Holzbaur; 2005", 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(10, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(4.0,"Holzbaur; 2005",'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(566.2,"Holzbaur; 2005",'N')\
											]}\
		}}
		\\\\\\\\\\\\\\\\
	FCR_Settings = {\
		'Notes' : [\
					"FCR EFE MA is not listed in Ramsay; 2009 but Pigeon has a quadratic function with R² = 0.9975. Pigeon only takes elbow angle into account. Coefficients and equation number/type are listed below to test either implementation. EFE MA was estimated to be constant and 10 mm for this muscle. If you use Pigeon, make sure to only accept positive moment arm values, as this model fails outside the ROM. One option is to set the MA to the constant value (i.e., MA[theta>140°] = MA[140°])"\
					]
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Pigeon; 1996", \
	 					"Sources" : \
							[\
								MA_Settings(Pigeon_coeff_conversion([0.9351,0.5375,-0.3627,0,0,0]), 'Pigeon; 1996', 'mm', None, 2.86, 'Elbow'),\
								MA_Settings(1, 'Est', 'cm', None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(129,125.7,1.0,"Banks; 2006"),\
		'Mass' : (28.7, "Banks; 2006", 'g'),\
		'Optimal Muscle Length' : {"Primary Source" : "Holzbaur; 2005", "Sources" : [(6.3, 'Holzbaur; 2005', 'cm')]},\
		'Optimal Tendon Length' : {	"Primary Source" : None, \
			"Sources" : \
				[\
					Input_Source(None, None, None)\
				]}, \
		'Pennation Angle' : {	"Primary Source" : None, \
		"Sources" : \
			[\
				Input_Source(None, None, None)\
			]}, \
		'PCSA' : {	"Primary Source" : None, \
		'Sources' : \
		[\
		Input_Source(None,None,None)\
		]}, \
		'Maximum Isometric Force': {	"Primary Source" : None, \
				'Sources': \
					[\
						Input_Source(None,None,None)\
					]}\
		}}

	ECRB_Settings = {\
		'Shoulder MA' : {"Primary Source" : "Est",\
		 					"Sources" : [MA_Settings(0, 'm', None, None, 'Shoulder', 'Est')]}, \
		'Elbow MA' : {"Primary Source" : "Ramsay; 2009", \
		 					"Sources" : [MA_Settings([-11.256,17.8548,1.6398,-0.5073,-2.8827,0,-0.0942,0,0,0,0,0,0,0,0,0], 'mm', 2, None, 'Elbow', 'Ramsay; 2009')]}, \
		'Spindle' : Spindle_Settings(102,132.7,0.77,"Banks; 2006"),\
		'Mass' : (32.1, "Banks; 2006", 'g'),\
		'Optimal Muscle Length' : {"Primary Source" : "Holzbaur; 2005", "Sources" : [(5.9, 'Holzbaur; 2005', 'cm')]},\
		'Optimal Tendon Length' : {	"Primary Source" : None, \
			"Sources" : \
				[\
					Input_Source(None, None, None)\
				]}, \
		'Pennation Angle' : {	"Primary Source" : None, \
		"Sources" : \
			[\
				Input_Source(None, None, None)\
			]}, \
		'PCSA' : {	"Primary Source" : None, \
		'Sources' : \
		[\
		Input_Source(None,None,None)\
		]}, \
		'Maximum Isometric Force': {	"Primary Source" : None, \
				'Sources': \
					[\
						Input_Source(None,None,None)\
					]}\
		}}

	"""
	ECRL EFE MA for Ramsay; 2009 has R² = 0.978 whereas Pigeon has R² = 0.9986. Pigeon, however, only takes elbow angle into account, whereas Ramsay; 2009 takes in variable PS angles. Additionally, Pigeon only considers elbow angles between 0 and 140 degrees and exhibits a decrease in MA as elbow angle approaches the upper bound of the ROM. This should (intiutively speaking) make the extensors MA largest, but Pigeon exhibits a drop off that may make it less favorable for movements at the boundary of the ROM. Coefficients and equation number/type are listed below to test either implementation.
	"""
	ECRL_Settings = {\
		'Shoulder MA' : {"Primary Source" : "Est",\
		 					"Sources" : [MA_Settings(0, 'm', None, None, 'Shoulder', 'Est')]}, \
		'Elbow MA' : {"Primary Source" : "Pigeon; 1996", \
		 					"Sources" : [MA_Settings(Pigeon_coeff_conversion([4.7304,1.2590,4.4347,-3.0229,0,0]), 'mm', None, None, 'Elbow', 'Pigeon; 1996'),\
							MA_Settings([-7.7034,16.3913,7.4361,-1.7566,0,-1.3336,0,0.0742,0,0,0,0,0,0,0,0], 'mm', 2, None, 'Elbow', 'Ramsay; 2009')]}, \
		'Spindle' : Spindle_Settings(74,155.2,0.48,"Banks; 2006"),\
		'Mass' : (44.3, "Banks; 2006", 'g'),\
		'Optimal Muscle Length' : {"Primary Source" : "Holzbaur; 2005", "Sources" : [(8.1, 'Holzbaur; 2005', 'cm')]},\
		'Optimal Tendon Length' : {	"Primary Source" : None, \
			"Sources" : \
				[\
					Input_Source(None, None, None)\
				]}, \
		'Pennation Angle' : {	"Primary Source" : None, \
		"Sources" : \
			[\
				Input_Source(None, None, None)\
			]}, \
		'PCSA' : {	"Primary Source" : None, \
		'Sources' : \
		[\
		Input_Source(None,None,None)\
		]}, \
		'Maximum Isometric Force': {	"Primary Source" : None, \
				'Sources': \
					[\
						Input_Source(None,None,None)\
					]}\
		}}

	FCU_Settings = {\
		'Shoulder MA' : {"Primary Source" : "Est",\
		 					"Sources" : [MA_Settings(0, 'm', None, None, 'Shoulder', 'Est')]}, \
		'Elbow MA' : {"Primary Source" : "Est", \
		 					"Sources" : [MA_Settings(1, 'cm', None, None, 'Elbow', 'Est')]}, \
		'Spindle' : Spindle_Settings(175,141.2,1.2,"Banks; 2006"),\
		'Mass' : (36.5,"Banks; 2006", 'g'),\
		'Optimal Muscle Length' : {"Primary Source" : "Holzbaur; 2005", "Sources" : [(5.1, 'Holzbaur; 2005', 'cm')]},\
		'Optimal Tendon Length' : {	"Primary Source" : None, \
			"Sources" : \
				[\
					Input_Source(None, None, None)\
				]}, \
		'Pennation Angle' : {	"Primary Source" : None, \
		"Sources" : \
			[\
				Input_Source(None, None, None)\
			]}, \
		'PCSA' : {	"Primary Source" : None, \
		'Sources' : \
		[\
		Input_Source(None,None,None)\
		]}, \
		'Maximum Isometric Force': {	"Primary Source" : None, \
				'Sources': \
					[\
						Input_Source(None,None,None)\
					]}\
		}}

	"""
	Note: only they muscle for the second digit was used for the FDS muscle. NEED TO DETERMINE IF THIS SHOULD BE A SUM OR AN AVERAGE FOR MEASURES LIKE PCSA, F_MAX, ETC.
	"""
	FDS_Settings = {\
		'Shoulder MA' : {"Primary Source" : "Est",\
		 					"Sources" : [MA_Settings(0, 'm', None, None, 'Shoulder', 'Est')]}, \
		'Elbow MA' : {"Primary Source" : "Est", \
		 					"Sources" : [MA_Settings(1, 'cm', None, None, 'Elbow', 'Est')]}, \
		'Spindle' : Spindle_Settings(356,224.9,1.6,"Banks; 2006"),\
		'Mass' : (95.2,"Banks; 2006", 'g'),\
		'Optimal Muscle Length' : {"Primary Source" : "Holzbaur; 2005", "Sources" : [(8.4, 'Holzbaur; 2005', 'cm')]},\
		'Optimal Tendon Length' : {	"Primary Source" : None, \
			"Sources" : \
				[\
					Input_Source(None, None, None)\
				]}, \
		'Pennation Angle' : {	"Primary Source" : None, \
		"Sources" : \
			[\
				Input_Source(None, None, None)\
			]}, \
		'PCSA' : {	"Primary Source" : None, \
		'Sources' : \
		[\
		Input_Source(None,None,None)\
		]}, \
		'Maximum Isometric Force': {	"Primary Source" : None, \
				'Sources': \
					[\
						Input_Source(None,None,None)\
					]}\
		}}

	PL_Settings = {\
		'Shoulder MA' : {"Primary Source" : "Est",\
		 					"Sources" : [MA_Settings(0, 'm', None, None, 'Shoulder', 'Est')]}, \
		'Elbow MA' : {"Primary Source" : "Est", \
		 					"Sources" : [MA_Settings(1, 'cm', None, None, 'Elbow', 'Est')]}, \
		'Spindle' : Spindle_Settings(None,None,None,None),\
		'Mass' : (None, "N/A","N/A"),\
		'Optimal Muscle Length' : {"Primary Source" : "Holzbaur; 2005", "Sources" : [(6.4, 'Holzbaur; 2005', 'cm')]},\
		'Optimal Tendon Length' : {	"Primary Source" : None, \
			"Sources" : \
				[\
					Input_Source(None, None, None)\
				]}, \
		'Pennation Angle' : {	"Primary Source" : None, \
		"Sources" : \
			[\
				Input_Source(None, None, None)\
			]}, \
		'PCSA' : {	"Primary Source" : None, \
		'Sources' : \
		[\
		Input_Source(None,None,None)\
		]}, \
		'Maximum Isometric Force': {	"Primary Source" : None, \
				'Sources': \
					[\
						Input_Source(None,None,None)\
					]}\
		}}

	"""
	ECU EFE MA is not listed in Ramsay; 2009 but Pigeon has a quadratic function with R² = 0.9966. Pigeon only takes elbow angle into account. Coefficients and equation number/type are listed below to test either implementation. EFE MA was estimated to be constant and -10 mm for this muscle. If you use Pigeon, make sure to only accept negative moment arm values, as this model fails outside the ROM. One option is to set the MA to the constant value (i.e., MA[theta>140°] = MA[140°])
	"""
	ECU_Settings = {\
		'Shoulder MA' : {"Primary Source" : "Est",\
		 					"Sources" : [MA_Settings(0, 'm', None, None, 'Shoulder', 'Est')]}, \
		'Elbow MA' : {"Primary Source" : "Pigeon; 1996", \
		 					"Sources" : [MA_Settings(Pigeon_coeff_conversion([-2.1826,-1.7386,1.1491,0,0,0]), 'mm', None, None, 'Elbow', 'Pigeon; 1996'),\
							MA_Settings(-1, 'cm', None, None, 'Elbow', 'Est')]}, \
		'Spindle' : Spindle_Settings(157,118,1.3,"Banks; 2006"),\
		'Mass' : (25.2,"Banks; 2006", 'g'),\
		'Optimal Muscle Length' : {"Primary Source" : "Holzbaur; 2005", "Sources" : [(6.2, 'Holzbaur; 2005', 'cm')]},\
		'Optimal Tendon Length' : {	"Primary Source" : None, \
			"Sources" : \
				[\
					Input_Source(None, None, None)\
				]}, \
		'Pennation Angle' : {	"Primary Source" : None, \
		"Sources" : \
			[\
				Input_Source(None, None, None)\
			]}, \
		'PCSA' : {	"Primary Source" : None, \
		'Sources' : \
		[\
		Input_Source(None,None,None)\
		]}, \
		'Maximum Isometric Force': {	"Primary Source" : None, \
				'Sources': \
					[\
						Input_Source(None,None,None)\
					]}\
		}}

	EDM_Settings = {\
		'Shoulder MA' : {"Primary Source" : "Est",\
		 					"Sources" : [MA_Settings(0, 'm', None, None, 'Shoulder', 'Est')]}, \
		'Elbow MA' : {"Primary Source" : "Est", \
		 					"Sources" : [MA_Settings(-1, 'cm', None, None, 'Elbow', 'Est')]}, \
		'Spindle' : Spindle_Settings(53,59.8,0.89,"Banks; 2006"),\
		'Mass' : (6.2, "Banks; 2006", 'g'),\
		'Optimal Muscle Length' : {"Primary Source" : "Holzbaur; 2005", "Sources" : [(6.8, 'Holzbaur; 2005', 'cm')]},\
		'Optimal Tendon Length' : {	"Primary Source" : None, \
			"Sources" : \
				[\
					Input_Source(None, None, None)\
				]}, \
		'Pennation Angle' : {	"Primary Source" : None, \
		"Sources" : \
			[\
				Input_Source(None, None, None)\
			]}, \
		'PCSA' : {	"Primary Source" : None, \
		'Sources' : \
		[\
		Input_Source(None,None,None)\
		]}, \
		'Maximum Isometric Force': {	"Primary Source" : None, \
				'Sources': \
					[\
						Input_Source(None,None,None)\
					]}\
		}}

	"""
	Note: only they muscle for the second digit was used for the EDC muscle. NEED TO DETERMINE IF THIS SHOULD BE A SUM OR AN AVERAGE FOR MEASURES LIKE PCSA, F_MAX, ETC.
	"""
	EDC_Settings = {\
		'Shoulder MA' : {"Primary Source" : "Est",\
		 					"Sources" : [MA_Settings(0, 'm', None, None, 'Shoulder', 'Est')]}, \
		'Elbow MA' : {"Primary Source" : "Est", \
		 					"Sources" : [MA_Settings(-1, 'cm', None, None, 'Elbow', 'Est')]}, \
		'Spindle' : Spindle_Settings(219,152.6,1.4,"Banks; 2006"),\
		'Mass' : (42.8, "Banks; 2006", 'g'),\
		'Optimal Muscle Length' : {"Primary Source" : "Holzbaur; 2005", "Sources" : [(7.0, 'Holzbaur; 2005', 'cm')]},\
		'Optimal Tendon Length' : {	"Primary Source" : None, \
			"Sources" : \
				[\
					Input_Source(None, None, None)\
				]}, \
		'Pennation Angle' : {	"Primary Source" : None, \
		"Sources" : \
			[\
				Input_Source(None, None, None)\
			]}, \
		'PCSA' : {	"Primary Source" : None, \
		'Sources' : \
		[\
		Input_Source(None,None,None)\
		]}, \
		'Maximum Isometric Force': {	"Primary Source" : None, \
				'Sources': \
					[\
						Input_Source(None,None,None)\
					]}\
		}}

	AN_Settings = {\
		'Shoulder MA' : {"Primary Source" : "Est",\
		 					"Sources" : [MA_Settings(0, 'm', None, None, 'Shoulder', 'Est')]}, \
		'Elbow MA' : {"Primary Source" : "Pigeon; 1996", \
		 					"Sources" : [MA_Settings(Pigeon_coeff_conversion([-5.3450,-2.2841,8.4297,-14.329,10.448,-2.736]), 'mm', None, None, 'Elbow', 'Pigeon; 1996')]}, \
		'Spindle' : Spindle_Settings(None,None,None,None),\
		'Mass' : (None, "N/A", "N/A"),\
		'Optimal Muscle Length' : {"Primary Source" : None, "Sources" : [(None,None,None)]},\
		'Optimal Tendon Length' : {	"Primary Source" : None, \
			"Sources" : \
				[\
					Input_Source(None, None, None)\
				]}, \
		'Pennation Angle' : {	"Primary Source" : None, \
		"Sources" : \
			[\
				Input_Source(None, None, None)\
			]}, \
		'PCSA' : {	"Primary Source" : None, \
		'Sources' : \
		[\
		Input_Source(None,None,None)\
		]}, \
		'Maximum Isometric Force': {	"Primary Source" : None, \
				'Sources': \
					[\
						Input_Source(None,None,None)\
					]}\
		}}

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
	if PreselectedMuscles==None:
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
	else:
		# assert type(PreselectedMuscles)==list and len(PreselectedMuscles)==8, "PreselectedMuscles, when used, must be a list of 8 numbers."
		assert np.array([type(MuscleNumber)==int for MuscleNumber in PreselectedMuscles]).all(),\
			"PreselectedMuscles must be a list of muscle numbers (ints)."
		assert np.array([MuscleNumber in range(1,len(AllAvailableMuscles)+1) \
			for MuscleNumber in PreselectedMuscles]).all(), \
				"PreselectedMuscles contains a muscle number outside the available muscles."
		SelectedMuscles = [AllAvailableMuscles[int(el)-1] \
								for el in PreselectedMuscles]
		MusclesToBeDeleted = [Muscle for Muscle in AllAvailableMuscles \
									if Muscle not in SelectedMuscles]
		for Muscle in MusclesToBeDeleted:
			del(AllMuscleSettings[Muscle])
	return(AllMuscleSettings)
def unit_angle_conversion(Value,Units):
	import numpy as np
	assert type(Units)==str, "Units must be a string."
	assert Units.capitalize() in ["Degrees","Deg","Degree","Radians","Radian","Rad"], "Can covert inches, cm, and mm to meters. Please use appropriate Units."

	if Units.capitalize() in ["Radians","Radian","Rad"]:
		return(Value)
	elif Units.capitalize() in ["Degrees","Deg","Degree"]:
		if type(Value)==list:
			return(list(np.array(Value)*np.pi/180))
		else:
			return(Value*np.pi/180)
def unit_length_conversion(Value,Units):
	assert type(Units)==str, "Units must be a string."
	assert Units.capitalize() in ["In","Inches","Cm","Centimeters","Centimeter","Mm","Millimeters","Millimeter","Meters","Meter","M"], "Can covert inches, cm, and mm to meters. Please use appropriate Units."

	if Units.capitalize() in ["Meter","Meters","M"]:
		return(Value)
	elif Units.capitalize() in ["In","Inches"]:
		if type(Value)==list:
			return(list(np.array(Value)*2.54/100))
		else:
			return(Value*2.54/100)
	elif Units.capitalize() in  ["Cm","Centimeters","Centimeter"]:
		if type(Value)==list:
			return(list(np.array(Value)/100))
		else:
			return(Value/100)
	elif Units.capitalize() in  ["Mm","Millimeters","Millimeter"]:
		if type(Value)==list:
			return(list(np.array(Value)/1000))
		else:
			return(Value/1000)
def unit_area_conversion(Value,Units):
	assert type(Units)==str, "Units must be a string."
	assert Units.capitalize() in ["Sq in","Squared inches","Inches squared","In sq","Cm sq","Sq cm","Centimeters squared","Squared centimeters","Mm sq","Sq mm","Millimeters squared","Squared millimeters","Meters squared","Squared meter","M sq","Sq m"], "Can covert inches², cm², and mm² to meters². Please use appropriate Units."

	if Units.capitalize() in ["Meters squared","Squared meter","M sq","Sq m"]:
		return(Value)
	elif Units.capitalize() in ["Sq in","Squared inches","Inches squared","In sq"]:
		if type(Value)==list:
			return(list(np.array(Value)*((2.54/100)**2)))
		else:
			return(Value*((2.54/100)**2))
	elif Units.capitalize() in ["Cm sq","Sq cm","Centimeters squared","Squared centimeters"]:
		if type(Value)==list:
			return(list(np.array(Value)/(100**2)))
		else:
			return(Value/(100**2))
	elif Units.capitalize() in ["Mm sq","Sq mm","Millimeters squared","Squared millimeters"]:
		if type(Value)==list:
			return(list(np.array(Value)/(1000**2)))
		else:
			return(Value/(1000**2))
def unit_mass_conversion(Value,Units):
	assert type(Units)==str, "Units must be a string."
	assert Units.capitalize() in ["Lbs","Lb","Pounds","G","Grams","Gram","Kg","Kilograms","Kilogram"], "Can covert lbs and g to kilograms. Please use appropriate Units."

	if Units.capitalize() in ["Kg","Kilograms","Kilogram"]:
		return(Value)
	elif Units.capitalize() in ["G","Grams","Gram"]:
		if type(Value)==list:
			return(list(np.array(Value)/1000))
		else:
			return(Value/1000)
	elif Units.capitalize() in  ["Lbs","Lb","Pounds"]:
		if type(Value)==list:
			return(list(np.array(Value)*0.45359237))
		else:
			return(Value*0.45359237)
def return_optimal_length(MuscleSettings):
	import numpy as np
	OptimalMuscleLengthsList = MuscleSettings["Optimal Muscle Length"]["Sources"]
	PrimarySource = MuscleSettings["Optimal Muscle Length"]["Primary Source"]
	return(OptimalMuscleLengthsList[np.where([src[1] == PrimarySource for src in OptimalMuscleLengthsList])[0][0]])
def return_MA_matrix_functions(AllMuscleSettings):
	import numpy as np
	import sympy as sp
	from sympy.utilities import lambdify
	def MA_function(Parameters):
		"""
		Note:

		Angles should be a number if Coefficients has a length of 5, or a list of length 2 when the Coefficients have lengths 16 or 18. Angles[0] will be the PRIMARY ANGLE for the DOF being considered while Angles[1] will be the secondary angle.

		Notes:

		threshold is only needed for Pigeon or Ramsay; 2009 MA functions that are invalid outside of a given value. Must be either None (default) or the radian value of the threshold.

		dof only needed for Pigeon (Ramsay; 2009 only handles EFE for this 2 DOF system). Must be either 'Shoulder' or 'Elbow'.

		eq is only needed for Ramsay; 2009 (Pigeon has one quintic polynomial). eq must be either 1, 2, or 3, with list length requirements of 5, 16, or 18, respectively.
		"""
		import sympy as sp
		import numpy as np
		ParameterList = Parameters["Sources"]
		PrimarySource = Parameters["Primary Source"]
		Parameters = ParameterList[np.where([params.Source == PrimarySource for params in ParameterList])[0][0]]
		assert str(type(Parameters))=="<class '__main__.MA_Settings'>", "Parameters are not in correct namedtuple form."
		src = Parameters.Source
		Coefficients = unit_length_conversion(Parameters.Values,Parameters.Units)
		eq = Parameters.Equation_Number
		dof = Parameters.DOF
		threshold = Parameters.Threshold

		global θ_SFE,θ_EFE,θ_PS
		assert type(src) == str, "src must be a str."
		assert src.capitalize() in ['Ramsay; 2009','Pigeon; 1996','Kuechle; 1997','Holzbaur; 2005', 'Est'], "src must be either Ramsay; 2009, Pigeon or Est (Estimate)."
		if dof != None:
			assert type(dof) == str, "dof must be a str."
			assert dof.capitalize() in ['Shoulder','Elbow'], "dof must be either Shoulder or Elbow."
		'''
		Note:
		For Kuechle and Holzbaur, where estimates or average MA were given, the format should be [MA,0,0,0,0,0] such that the function returns a constant MA function (See matrix multiplication below).
		'''
		if src.capitalize() in ['Pigeon; 1996', 'Kuechle; 1997', 'Holzbaur; 2005']:
			assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
			assert dof != None, "For Pigeon (1996), dof must be stated."
			eq = None
			if dof.capitalize() == 'Elbow MA':
				θ = θ_EFE
			else:
				θ = θ_SFE
			MomentArm = (np.matrix(Coefficients,dtype='float64')\
							*np.matrix([1,θ,θ**2,θ**3,θ**4,θ**5]).T)[0,0]
		elif src.capitalize() == 'Est' or src.capitalize() == 'Holzbaur; 2005':
			MomentArm = np.array(Coefficients,dtype='float64')
		else: #src.capitalize() == 'Ramsay; 2009'
			θ = θ_EFE
			assert type(Coefficients) == list, "Coefficients must be a list."
			assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
			assert eq in [1,2,3], "eq must be either 1, 2, or 3 when using Ramsay; 2009 (2009)."
			if eq == 1:
				assert len(Coefficients) == 5, "For Eq. 1, Coefficients must be 5 elements long."
				MomentArm = (sp.Matrix(Coefficients,dtype='float64').T\
								*sp.Matrix([1,θ,θ**2,θ**3,θ**4]))[0,0]
			elif eq == 2:
				assert len(Coefficients)==16, "For Eq. 2, Coefficients must be 16 elements long."
				MomentArm = (sp.Matrix(Coefficients,dtype='float64').T*\
								sp.Matrix([1, θ, θ_PS, θ*θ_PS, θ**2, \
											θ_PS**2, (θ**2)*θ_PS, θ*(θ_PS**2), \
											(θ**2)*(θ_PS**2), θ**3, θ_PS**3, \
											(θ**3)*θ_PS, θ*(θ_PS**3), \
											(θ**3)*(θ_PS**2), (θ**2)*(θ_PS**3), \
											(θ**3)*(θ_PS**3)]))[0, 0]
			else: # eq == 3
				assert len(Coefficients)==18, "For Eq. 3, Coefficients must be 18 elements long."
				MomentArm = (sp.Matrix(Coefficients,dtype='float64').T*\
								sp.Matrix([1, θ, θ_PS, θ*θ_PS, θ**2, \
									θ_PS**2, (θ**2)*θ_PS, θ*(θ_PS**2), (θ**2)*(θ_PS**2), \
									θ**3, (θ**3)*θ_PS, (θ**3)*(θ_PS**2), \
									θ**4, (θ**4)*θ_PS, (θ**4)*(θ_PS**2),  \
									θ**5, (θ**5)*θ_PS, (θ**5)*(θ_PS**2)]))[0, 0]
		if threshold == None:
			return(MomentArm)
		else:
			assert type(threshold) in [int,float], "threshold must be a number."
			MomentArm = sp.Piecewise((MomentArm,θ<threshold),(MomentArm.subs(θ,threshold),θ>=threshold))
			return(MomentArm)

	MuscleList = AllMuscleSettings.keys()

	RT_symbolic = sp.Matrix([MA_function(AllMuscleSettings[muscle]['Elbow MA']) for muscle in MuscleList])
	dRT_symbolic = sp.Matrix(sp.diff(RT_symbolic,θ_EFE))
	d2RT_symbolic = sp.Matrix(sp.diff(sp.diff(RT_symbolic,θ_EFE),θ_EFE))
	# RT_func = lambdify([θ_SFE,x1,θ_PS],RT_symbolic)
	# dRT_func = lambdify([θ_SFE,x1,θ_PS],dRT_symbolic)
	# import ipdb; ipdb.set_trace()
	# returns an (n,m) matrix when n is the number of muscles and m is the number of DOFS. We chose to return R.T because this is commonly utilized in muscle velocity calculations.
	return(RT_symbolic,dRT_symbolic,d2RT_symbolic)
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
		if i==1:
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

"""
################################
######## Tension Driven ########
################################

x_1 &= \theta \\
x_2 &= \dot{\theta} \\
u_1 &= T_{1} \\
u_2 &= T_{2} \\

################################
#### Muscle Velocity Driven ####
################################

x_1 &= \theta \\
x_2 &= \dot{\theta} \\
x_3 &= T_{1} \\
x_4 &= T_{2} \\
u_1 &= l_{m,1} \\
u_2 &= l_{m,2} \\

################################
###### Activation Driven #######
################################

x_1 &= \theta \\
x_2 &= \dot{\theta} \\
x_3 &= T_{1} \\
x_4 &= T_{2} \\
x_5 &= l_{m,1} \\
x_6 &= l_{m,2} \\
x_7 &= v_{m,1} \\
x_8 &= v_{m,2} \\
u_1 &= \alpha_1 \\
u_2 &= \alpha_2 \\

"""

N = 20001
Time = np.linspace(0,20,N)
dt = Time[1]-Time[0]

AllMuscleSettings = return_muscle_settings(PreselectedMuscles=[5,6])

g,L = 9.80, 0.45 #m/s², m
M = 2 # kg

α1 = 0 # 10*np.pi/180 # rads
α2 = 0 # 10*np.pi/180 # rads

m1 = 1 # kg
m2 = 1 # kg

"""
There was some debate regarding damping terms. No explicit value was given. Loeb simplifies the equation by utilizing the F_{PE_{2,i}} damping term (η) instead, as this is added to B_{m,i}*(v_{m,i}/l_{o,i}) anyways. This η is very small (0.01), so the damping is not significant. Might need to do sensitivity analysis on these two values (currently set to zero) (06/16/2018).
"""

bm1 = 0 # kg/s
bm2 = 0 # kg/s

cT = 27.8
kT = 0.0047
LrT = 0.964

β = 1.55
ω = 0.75
ρ = 2.12

V_max = -9.15
cv0 = -5.78
cv1 = 9.18
av0 = -1.53
av1 = 0
av2 = 0
bv = 0.69

c_1 = 23.0
k_1 = 0.046
Lr1 = 1.17
η = 0.01

lo1 = return_primary_source(AllMuscleSettings["BIC"]["Optimal Muscle Length"])
lo2 = return_primary_source(AllMuscleSettings["TRI"]["Optimal Muscle Length"])

lo1 = unit_length_conversion(lo1[0],lo1[2])
lo2 = unit_length_conversion(lo2[0],lo2[2])

######### NEED OPTIMAL TENDON LENGTHS FOR BIC/TRI #########
lTo1 = (1)*lo1
lTo2 = (1)*lo2
###########################################################

[[r1,r2],[dr1,dr2],[d2r1,d2r2]]= return_MA_matrix_functions(AllMuscleSettings)
PCSA1 = 25.90 # cm² (Garner & Pandy; 2003)
PCSA2 = 76.30 # cm² (Garner & Pandy; 2003)
SpecificTension = 330 #kPa (Garner & Pandy; 2003)
kPa_To_N = 0.1
F_MAX1 = 849.29 # (Garner & Pandy; 2003) ~ PCSA1*SpecificTension*kPa_To_N (in N)
F_MAX2 = 2332.92 # (Garner & Pandy; 2003) ~ PCSA2*SpecificTension*kPa_To_N (in N)

Amp = 7.5*np.pi/180
Base = 90*np.pi/180
Freq = 2*np.pi

k1,k2,k3,k4 = 100,100,10,100

MaxStep_Tension = 0.01 # percentage of positive maximum.
Tension_Bounds = [[0,F_MAX1],[0,0.10*F_MAX2]]

MaxStep_MuscleVelocity = 5 # percentage of positive maximum.
MuscleVelocity_Bounds =[[-2*lo1,2*lo1],[-0.2*lo2,0.2*lo2]]

MaxStep_Activation = 0.1 # percentage of positive maximum (1)
Activation_Bounds = [[0,1],[0,1]]

"""
c_{1} &= -\frac{3g}{2L} \\
c_{2} &= \frac{3}{ML^2} \\
c_{3} &= \cos(\rho_1) \\
c_{4} &= \cos(\rho_2)
c_{5} &= \frac{\cos(\alpha_{1})}{m_1} \\
c_{6} &= \frac{\cos^2(\alpha_{1})}{m_1} \\
c_{7} &= \frac{b_{m,1}\cos^2(\alpha_{1})}{m_1} \\
c_{8} &= \tan^2(\alpha_{1}) \\
c_{9} &= \frac{\cos(\alpha_{2})}{m_2} \\
c_{10} &= \frac{\cos^2(\alpha_{2})}{m_2} \\
c_{11} &= \frac{b_{m,2}\cos^2(\alpha_{2})}{m_2} \\
c_{12} &= \tan^2(\alpha_{2}) \\

"""

c1 = -(3*g)/(2*L)
c2 = 3/(M*L**2)
c3 = np.cos(α1)
c4 = np.cos(α2)
c5 = np.cos(α1)/m1
c6 = F_MAX1*np.cos(α1)**2/m1
c7 = F_MAX1*bm1*np.cos(α1)**2/(m1*lo1)
c8 = np.tan(α1)**2
c9 = np.cos(α2)/m2
c10 = F_MAX2*np.cos(α2)**2/m2
c11 = F_MAX2*bm2*np.cos(α2)**2/(m2*lo2)
c12 = np.tan(α2)**2

'''
R_{1} &= r_{1}\left(\theta,\theta_{PS}=\frac{\pi}{2}\right) \\
R_{2} &= r_{2}\left(\theta,\theta_{PS}=\frac{\pi}{2}\right)  \\
K_{T,1} &= \frac{F_{max,1}c^{T}}{l_{T,o,1}}\left(1 - \exp{\left(\frac{-T_1}{F_{max,1}c^{T}k^{T}}\right)}\right) \hspace{1em} \text{(Variable stiffness coefficient from Zajac (1989) ODE for tendon force)} \\
v_{MTU,1} &= \text{sgn}\left(-r_1(\theta)\right)\cdot\dot{\theta}\cdot\sqrt{\left(\frac{\partial r_1}{\partial\theta}\right)^2 + r_1^2(\theta)} \\
K_{T,2} &= \frac{F_{max,2}c^{T}}{l_{T,o,2}}\left(1 - \exp{\left(\frac{-T_2}{F_{max,2}c^{T}k^{T}}\right)}\right) \hspace{1em} \text{(Variable stiffness coefficient from Zajac (1989) ODE for tendon force)} \\
v_{MTU,2} &= \text{sgn}\left(-r_2(\theta)\right)\cdot\dot{\theta}\cdot\sqrt{\left(\frac{\partial r_2}{\partial\theta}\right)^2 + r_2^2(\theta)} \\
F_{LV,1} &= f_{L,1}(l_{m,1}) \cdot f_{V,1}(l_{m,1},v_{m,1}) \\
F_{LV,2} &= f_{L,2}(l_{m,2}) \cdot f_{V,2}(l_{m,2},v_{m,2}) \\
'''

r1 = lambdify([θ_EFE],r1.subs([(θ_PS,np.pi/2)]))
r2 = lambdify([θ_EFE],r2.subs([(θ_PS,np.pi/2)]))
dr1_dθ = lambdify([θ_EFE],dr1.subs([(θ_PS,np.pi/2)]))
dr2_dθ = lambdify([θ_EFE],dr2.subs([(θ_PS,np.pi/2)]))
d2r1_dθ2 = lambdify([θ_EFE],d2r1.subs([(θ_PS,np.pi/2)]))
d2r2_dθ2 = lambdify([θ_EFE],d2r2.subs([(θ_PS,np.pi/2)]))
# FL = lambda l,lo: 1
# FV = lambda l,v,lo: 1
FL = lambda l,lo: np.exp(-abs(((l/lo)**β-1)/ω)**ρ)
FV = lambda l,v,lo: np.piecewise(v,[v<=0, v>0],\
	[lambda v: (V_max - v/lo)/(V_max + (cv0 + cv1*(l/lo))*(v/lo)),\
	lambda v: (bv-(av0 + av1*(l/lo) + av2*(l/lo)**2)*(v/lo))/(bv + (v/lo))])

def R1(X):
	return(r1(X[0])) #
def dR1_dx1(X):
	return(dr1_dθ(X[0]))
def d2R1_dx12(X):
	return(d2r1_dθ2(X[0]))
def R2(X):
	return(r2(X[0])) #
def dR2_dx1(X):
	return(dr2_dθ(X[0]))
def d2R2_dx12(X):
	return(d2r2_dθ2(X[0]))
def KT_1(X):
	return((F_MAX1*cT/lTo1)*(1-np.exp(-X[2]/(F_MAX1*cT*kT)))) # NOT NORMALIZED (in N/m)
def dKT_1_dx3(X):
	return(1/(kT*lTo1)*np.exp(-X[2]/(F_MAX1*cT*kT))) # NOT NORMALIZED (in N/m)
def v_MTU1(X):
	return(np.sign(-R1(X))*X[1]*np.sqrt(dR1_dx1(X)**2 + R1(X)**2)) # NOT NORMALIZED (in m/s)
def KT_2(X):
	return((F_MAX2*cT/lTo2)*(1-np.exp(-X[3]/(F_MAX2*cT*kT)))) # NOT NORMALIZED (in N/m)
def v_MTU2(X):
	return(np.sign(-R2(X))*X[1]*np.sqrt(dR2_dx1(X)**2 + R2(X)**2)) # NOT NORMALIZED (in m/s)
def FLV_1(X):
	return(FL(X[4],lo1)*FV(X[4],X[6],lo1))
def FLV_2(X):
	return(FL(X[5],lo2)*FV(X[5],X[7],lo2))

"""
################################
######## Tension Driven ########
################################

\dot{x}_1 &= x_{2} \\
\dot{x}_2 &= c_{1}\sin(x_{1}) + c_{2}R_{1}u_{1} - c_{2}R_{2}u_{2} \\
u_1 &= T_{1} \\
u_2 &= T_{2} \\

################################
#### Muscle Velocity Driven ####
################################

\dot{x}_1 &= x_{2} \\
\dot{x}_2 &= c_{1}\sin(x_{1}) + c_{2}R_{1}x_{3} - c_{2}R_{2}x_{4} \\
\dot{x}_3 &= K_{T,1}(v_{MTU,1} - c_{3}u_1) \\
\dot{x}_4 &= K_{T,2}(v_{MTU,2} - c_{4}u_2) \\
u_1 &= \dot{l}_{m,1} \\
u_2 &= \dot{l}_{m,2} \\

"""

def dX1_dt(X):
	return(X[1])
def dX2_dt(X,U=None):
	if U==None:
		return(c1*np.sin(X[0]) + c2*R1(X)*X[2] + c2*R2(X)*X[3])
	else:
		return(c1*np.sin(X[0]) + c2*R1(X)*U[0] + c2*R2(X)*U[1])
def dX3_dt(X,U=None):
	if U == None:
		return(KT_1(X)*(v_MTU1(X) - c3*X[6]))
	else:
		return(KT_1(X)*(v_MTU1(X) - c3*U[0]))
def dX4_dt(X,U=None):
	if U == None:
		return(KT_2(X)*(v_MTU2(X) - c4*X[7]))
	else:
		return(KT_2(X)*(v_MTU2(X) - c4*U[1]))

r = lambda t: Amp*np.sin(Freq*t) + Base
dr = lambda t: Amp*Freq*np.cos(Freq*t)
d2r = lambda t: -Amp*Freq**2*np.sin(Freq*t)
d3r = lambda t: -Amp*Freq**3*np.cos(Freq*t)

def return_constraint_variables_tension_driven(t,X):
	def Z1(t,X):
		return(r(t) - X[0])
	def dZ1(t,X):
		return(dr(t) - dX1_dt(X))
	def A1(t,X):
		return(dr(t) + k1*Z1(t,X))
	def dA1(t,X):
		return(d2r(t) + k1*dZ1(t,X))
	def Z2(t,X):
		return(X[1] - A1(t,X))
	"""
	def dZ2(t,X,U):
		return(c1*np.sin(X[0]) + c2*R1(X)*U[0] + c2*R2(X)*U[1] - dA1(t,X))
	"""
	def A2(t,X):
		return(Z1(t,X) + dA1(t,X) - c1*np.sin(X[0]) - k2*Z2(t,X))
	Coefficient1 = c2*R1(X)
	Coefficient2 = c2*R2(X)
	Constraint = A2(t,X)
	return(Coefficient1,Coefficient2,Constraint)
def return_constraint_variables_muscle_velocity_driven(t,X):
	def Z1(t,X):
		return(r(t) - X[0])
	def dZ1(t,X):
		return(dr(t) - dX1_dt(X))
	def d2Z1(t,X):
		return(d2r(t) - dX2_dt(X))
	def A1(t,X):
		return(dr(t) + k1*Z1(t,X))
	def dA1(t,X):
		return(d2r(t) + k1*dZ1(t,X))
	def d2A1(t,X):
		return(d3r(t) + k1*d2Z1(t,X))
	def Z2(t,X):
		return(X[1] - A1(t,X))
	def dZ2(t,X):
		return(dX2_dt(X) - dA1(t,X))
	def A2(t,X):
		return(Z1(t,X) + dA1(t,X) - c1*np.sin(X[0]) - k2*Z2(t,X))
	def dA2(t,X):
		return(dZ1(t,X) + d2A1(t,X) - c1*np.cos(X[0])*dX1_dt(X) - k2*dZ2(t,X))
	def Z3(t,X):
		return(c2*R1(X)*X[2] + c2*R2(X)*X[3] - A2(t,X))
	def dZ3(t,X):
		g1 = R1(X)
		g2 = R2(X)
		g3 = KT_1(X)
		g5 = KT_2(X)
		return(c2*dR1_dx1(X)*X[1]*X[2] + c2*dR2_dx1(X)*X[1]*X[3] \
						+ c2*g1*g3*v_MTU1(X) - c2*c3*g1*g3*U[0] \
							+ c2*g2*g5*v_MTU2(X) - c2*c4*g2*g5*U[1] \
								- dA2(t,X))
	def A3(t,X):
		return(Z2(t,X) - dA2(t,X) + k3*Z3(t,X) \
		+ c2*dR1_dx1(X)*X[1]*X[2] + 	c2*dR2_dx1(X)*X[1]*X[3] \
				+ c2*R1(X)*KT_1(X)*v_MTU1(X) + c2*R2(X)*KT_2(X)*v_MTU2(X))
	Coefficient1 = c2*c3*R1(X)*KT_1(X)
	Coefficient2 = c2*c4*R2(X)*KT_2(X)
	Constraint = A3(t,X)
	return(Coefficient1,Coefficient2,Constraint)
def return_U_tension_driven(t,X,U,dt,MaxStep,Bounds,Noise):
	import random
	Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_tension_driven(t,X)
	# elif Method == "Muscle Velocity":
	# 		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t,X)
	assert np.shape(Bounds)==(2,2), "Bounds must be (2,2)."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."
	if Constraint1 != 0:
		assert Coefficient1!=0 and Coefficient2!=0, "Error with Coefficients. Shouldn't be zero with nonzero constraint."
	else:
		assert Coefficient1!=0 and Coefficient2!=0, "Error with Constraint. 0 = 0 implies all inputs valid."

	if Coefficient1 == 0:
		LowerBound = Bounds[0][0]
		UpperBound = Bounds[0][1]
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
	elif Coefficient2 == 0:
		LowerBound = Constraint1/Coefficient1
		UpperBound = Constraint1/Coefficient1
		FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
		FeasibleInput2 = (Bounds[1][1]-Bounds[1][0])*np.random.rand(1000) + Bounds[1][0]
	elif np.sign(-Coefficient1) == np.sign(Coefficient2):
		if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
			LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
		else:
			LowerBound = Bounds[0][0]

		if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
			UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
		else:
			UpperBound = Bounds[0][1]
		assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in FeasibleInput1])
	else: # np.sign(-Coefficient1) != np.sign(Coefficient2)
		if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
			LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
		else:
			LowerBound = Bounds[0][0]

		if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
			UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
		else:
			UpperBound = Bounds[0][1]
		assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in FeasibleInput1])
	#
	# if Constraint1 > 0:
	# 	assert Bounds[0] >= Constraint1/Coefficient1, "Tension cannot be made by muscle 1."
	# 	LowerBound = Constraint1/Coefficient1
	# 	if (Constraint1 - Coefficient1*Bounds[0])/Coefficient2 > Bounds[1]:
	# 		UpperBound = (Constraint1 - Coefficient2*Bounds[1])/Coefficient1
	# 	else:
	# 		UpperBound = Bounds[0]
	# elif Constraint1 < 0:
	# 	assert Bounds[1] >= Constraint1/Coefficient2, "Tension cannot be made by muscle 2."
	# 	LowerBound = 0
	# 	if (Constraint1 - Coefficient2*Bounds[1])/Coefficient1 > Bounds[0]:
	# 		UpperBound = Bounds[0]
	# 	else:
	# 		UpperBound = (Constraint1 - Coefficient2*Bounds[1])/Coefficient1
	# else: # Constraint1 == 0
	# 	LowerBound = 0
	# 	if -Coefficient1*Bounds[0]/Coefficient2 > Bounds[1]:
	# 		UpperBound = -Coefficient2*Bounds[1]/Coefficient1
	# 	else:
	# 		UpperBound = Bounds[0]
	"""
	Checking to see which inputs have the appropriate allowable step size.
	"""
	euclid_dist = np.array(list(map(lambda u1,u2: np.sqrt(((U[0]-u1)/Bounds[0][1])**2 + ((U[1]-u2)/Bounds[1][1])**2),\
							FeasibleInput1,FeasibleInput2)))

	if t<10*dt: MaxStep = 10*MaxStep
	feasible_index = np.where(euclid_dist<=MaxStep)
	# elif Method == "Muscle Velocity":
	# 	feasible_index = np.where(np.logical_and(np.logical_and(euclid_dist>=0, euclid_dist<=MaxStep),np.sign(FeasibleInput1)!=np.sign(FeasibleInput2)))
	if len(feasible_index[0]) == 0: import ipdb; ipdb.set_trace()
	next_index = random.choice(feasible_index[0])
	u1 = FeasibleInput1[next_index]
	u2 = FeasibleInput2[next_index]
	return([u1,u2])
def return_U_muscle_velocity_driven(t,X,U,dt,MaxStep,Bounds,Noise):
	"""
	Enforcing a hyperbolic domain constraint to allow for realistic lengthening/shortenting relationships.
	Input2 = (lo1*0.001)*(lo2*0.001)/Input1 = lo1*lo2/(10^6*Input1)
	"""

	import random
	Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t,X)
	assert np.shape(Bounds)==(2,2), "Bounds must be (2,2)."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."
	if Constraint1 != 0:
		assert Coefficient1!=0 or Coefficient2!=0, "Error with Coefficients. Shouldn't be zero with nonzero constraint."
	else:
		assert Coefficient1!=0 or Coefficient2!=0, "Error with Constraint. 0 = 0 implies all inputs valid."

	if abs(Coefficient1) <= 1e-7:
		LowerBound = Bounds[0][0]
		UpperBound = Bounds[0][1]
		if Constraint1/Coefficient2 > 0:
			LowerBound = Bounds[0][0]
			UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
		else:
			LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
			UpperBound = Bounds[0][1]
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
	elif abs(Coefficient2) <= 1e-7:
		LowerBound = Constraint1/Coefficient1
		UpperBound = Constraint1/Coefficient1
		FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
		if Constraint1/Coefficient1 < 0:
			LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
			UpperBound = Bounds[1][1]
		else:
			LowerBound = Bounds[1][0]
			UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
		FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
	elif np.sign(-Coefficient1) == np.sign(Coefficient2): # DIFFERENT SIGNS
		if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
			LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
		else:
			LowerBound = Bounds[0][0]

		if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
			UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
		else:
			UpperBound = Bounds[0][1]
		assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
		HyperbolicBounds = np.sort([(Constraint1 - \
										np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
											/(2*Coefficient1), \
								 	(Constraint1 + \
										np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
											/(2*Coefficient1)])
		LowerBound = max([LowerBound,HyperbolicBounds[0]])
		UpperBound = min([UpperBound,HyperbolicBounds[1]])
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in FeasibleInput1])
	else: # np.sign(-Coefficient1) != np.sign(Coefficient2) SAME SIGN
		if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
			LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
		else:
			LowerBound = Bounds[0][0]

		if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
			UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
		else:
			UpperBound = Bounds[0][1]
		assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
		if Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001) > 0:
			HyperbolicBounds = np.sort([(Constraint1 - \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1), \
									 	(Constraint1 + \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1)])

			assert LowerBound < HyperbolicBounds[0] or HyperbolicBounds[1] < UpperBound, "No feasible solutions."

			FeasibleInput1 = []
			while len(FeasibleInput1)<1000:
				Random1 = (UpperBound-LowerBound)*np.random.rand() + LowerBound
				if Random1<HyperbolicBounds[0] or Random1>HyperbolicBounds[1]: FeasibleInput1.append(Random1)
			FeasinbleInput1 = np.array(FeasibleInput1)
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		else:
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
	def plot_constraints():
		import matplotlib.pyplot as plt
		plt.figure()
		Input1 = np.linspace(LowerBound,UpperBound,1001)
		Input2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in Input1])
		plt.plot(Input1,Input2,'k--')
		plt.plot(Input1, (lo1*0.001)*(lo2*0.001)/Input1,'r')
		plt.scatter(FeasibleInput1,FeasibleInput2,c='g',marker = '.')
		plt.ylim(MuscleVelocity_Bounds[1])
		plt.show()

	"""
	Checking to see which inputs have the appropriate allowable step size. In normalized muscle velocity.
	"""
	euclid_dist = np.array(list(map(lambda u1,u2: np.sqrt(((U[0]-u1)/lo1)**2 + ((U[1]-u2)/lo2)**2),\
							FeasibleInput1,FeasibleInput2)))

	if t<30*dt: MaxStep = 10*MaxStep
	feasible_index = np.where(euclid_dist<=MaxStep)

	if len(feasible_index[0]) == 0: import ipdb; ipdb.set_trace()
	next_index = random.choice(feasible_index[0])
	u1 = FeasibleInput1[next_index]
	u2 = FeasibleInput2[next_index]
	return([u1,u2])
def plot_MA_values(Time,x1):
	import matplotlib.pyplot as plt
	import numpy as np

	fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8,6))
	plt.subplots_adjust(hspace=0.1,bottom=0.1)

	plt.suptitle("Moment arm equations")
	ax1.plot(np.linspace(0,np.pi*(160/180),1001),\
				np.array(list(map(lambda x1: R1([x1]),np.linspace(0,np.pi*(160/180),1001)))),\
				'0.70')
	ax1.plot(np.linspace(min(x1),max(x1),101),\
				np.array(list(map(lambda x1: R1([x1]),np.linspace(min(x1),max(x1),101)))),\
				'g',lw=3)
	ax1.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
	ax1.set_xticklabels([""]*len(ax1.get_xticks()))
	ax1.set_ylabel("Moment Arm for\n Muscle 1 (m)")

	ax2.plot(Time,np.array(list(map(lambda x1: R1([x1]),x1))),'g')
	ax2.set_ylim(ax1.get_ylim())
	ax2.set_yticks(ax1.get_yticks())
	ax2.set_yticklabels([""]*len(ax1.get_yticks()))
	ax2.set_xticklabels([""]*len(ax2.get_xticks()))

	ax3.plot(np.linspace(0,np.pi*(160/180),1001),\
				np.array(list(map(lambda x1: R2([x1]),np.linspace(0,np.pi*(160/180),1001)))),\
				'0.70')
	ax3.plot(np.linspace(min(x1),max(x1),101),\
				np.array(list(map(lambda x1: R2([x1]),np.linspace(min(x1),max(x1),101)))),\
				'r',lw=3)
	ax3.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
	ax3.set_xticklabels([r"$0$",r"$\frac{\pi}{4}$",r"$\frac{\pi}{2}$",r"$\frac{3\pi}{4}$",r"$\pi$"])
	ax3.set_xlabel("Joint Angle (rads)")
	ax3.set_ylabel("Moment Arm for\n Muscle 2 (m)")

	ax4.plot(Time,np.array(list(map(lambda x1: R2([x1]),x1))),'r')
	ax4.set_ylim(ax3.get_ylim())
	ax4.set_yticks(ax3.get_yticks())
	ax4.set_yticklabels([""]*len(ax3.get_yticks()))
	ax4.set_xlabel("Time (s)")
	return(fig,[ax1,ax2,ax3,ax4])
def animate_test_2(response,t,x1,x2,x3,x4,u1,u2,dt,Bounds):
	assert type(response)==bool, "Input must be either True or False."

	if response == True:
		import numpy as np
		import matplotlib.pyplot as plt
		import matplotlib.animation as animation
		import time

		fig = plt.figure(figsize=(10,8))
		ax1 = plt.gca()

		DescriptiveTitle = "Plotting Constraints vs. Time"

		ax1.set_title(DescriptiveTitle,Fontsize=20,y=0.975)

		#Hyperbolic Constraint
		Input1 = list(np.linspace(Bounds[0][0],Bounds[0][1],1000001))
		Input1.remove(0)
		Input1 = np.array(Input1)
		ax1.plot(Input1,lo1*lo2*0.001**2/Input1,'r',lw=2)
		Input2 = list(np.linspace(Bounds[0][0],Bounds[0][1],1000001))
		Input2.remove(0)
		Input2 = np.array(Input2)
		ax1.plot(lo1*lo2*0.001**2/Input2,Input2,'r',lw=2)
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[0],[x1[0],x2[0],x3[0],x4[0]])
		if abs(Coefficient1) <= 1e-7:
			LowerBound = Bounds[0][0]
			UpperBound = Bounds[0][1]
			if Constraint1/Coefficient2 > 0:
				LowerBound = Bounds[0][0]
				UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
			else:
				LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
				UpperBound = Bounds[0][1]
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
		elif abs(Coefficient2) <= 1e-7:
			LowerBound = Constraint1/Coefficient1
			UpperBound = Constraint1/Coefficient1
			FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
			if Constraint1/Coefficient1 < 0:
				LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
				UpperBound = Bounds[1][1]
			else:
				LowerBound = Bounds[1][0]
				UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
			FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		elif np.sign(-Coefficient1) == np.sign(Coefficient2): # DIFFERENT SIGNS
			if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
				LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
			else:
				LowerBound = Bounds[0][0]

			if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
				UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
			else:
				UpperBound = Bounds[0][1]
			assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
			HyperbolicBounds = np.sort([(Constraint1 - \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1), \
									 	(Constraint1 + \
											np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
												/(2*Coefficient1)])
			LowerBound = max([LowerBound,HyperbolicBounds[0]])
			UpperBound = min([UpperBound,HyperbolicBounds[1]])
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		else: # np.sign(-Coefficient1) != np.sign(Coefficient2) SAME SIGN
			if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
				LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
			else:
				LowerBound = Bounds[0][0]

			if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
				UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
			else:
				UpperBound = Bounds[0][1]
			assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
			if Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001) > 0:
				HyperbolicBounds = np.sort([(Constraint1 - \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1), \
										 	(Constraint1 + \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1)])

				assert LowerBound < HyperbolicBounds[0] or HyperbolicBounds[1] < UpperBound, "No feasible solutions."

				FeasibleInput1 = []
				while len(FeasibleInput1)<1000:
					Random1 = (UpperBound-LowerBound)*np.random.rand() + LowerBound
					if Random1<HyperbolicBounds[0] or Random1>HyperbolicBounds[1]: FeasibleInput1.append(Random1)
				FeasinbleInput1 = np.array(FeasibleInput1)
				FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
										for el in FeasibleInput1])
			else:
				FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
				FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
										for el in FeasibleInput1])
		cline, = plt.plot(FeasibleInput1,FeasibleInput2,'b',lw=2)
		TimeText = plt.text(0.1,0.1,"t = " + str(t[0]),fontsize=16)
		chosenpoint, = plt.plot(u1[0],u2[0],c='k',marker='o')
		ax1.set_xlim(Bounds[0])
		ax1.set_ylim(Bounds[1])
		ax1.spines['right'].set_visible(False)
		ax1.spines['top'].set_visible(False)
		ax1.set_aspect('equal')

		def animate(i):
			Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[i],[x1[i],x2[i],x3[i],x4[i]])
			if abs(Coefficient1) <= 1e-7:
				LowerBound = Bounds[0][0]
				UpperBound = Bounds[0][1]
				if Constraint1/Coefficient2 > 0:
					LowerBound = Bounds[0][0]
					UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
				else:
					LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
					UpperBound = Bounds[0][1]
				FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
				FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
			elif abs(Coefficient2) <= 1e-7:
				LowerBound = Constraint1/Coefficient1
				UpperBound = Constraint1/Coefficient1
				FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
				if Constraint1/Coefficient1 < 0:
					LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
					UpperBound = Bounds[1][1]
				else:
					LowerBound = Bounds[1][0]
					UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
				FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			elif np.sign(-Coefficient1) == np.sign(Coefficient2): # DIFFERENT SIGNS
				if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
					LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
				else:
					LowerBound = Bounds[0][0]

				if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
					UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
				else:
					UpperBound = Bounds[0][1]
				assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
				HyperbolicBounds = np.sort([(Constraint1 - \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1), \
										 	(Constraint1 + \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1)])
				LowerBound = max([LowerBound,HyperbolicBounds[0]])
				UpperBound = min([UpperBound,HyperbolicBounds[1]])
				FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
				FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
										for el in FeasibleInput1])
			else: # np.sign(-Coefficient1) != np.sign(Coefficient2) SAME SIGN
				if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
					LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
				else:
					LowerBound = Bounds[0][0]

				if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
					UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
				else:
					UpperBound = Bounds[0][1]
				assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
				if Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001) > 0:
					HyperbolicBounds = np.sort([(Constraint1 - \
													np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
														/(2*Coefficient1), \
											 	(Constraint1 + \
													np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
														/(2*Coefficient1)])

					assert LowerBound < HyperbolicBounds[0] or HyperbolicBounds[1] < UpperBound, "No feasible solutions."

					FeasibleInput1 = []
					while len(FeasibleInput1)<1000:
						Random1 = (UpperBound-LowerBound)*np.random.rand() + LowerBound
						if Random1<HyperbolicBounds[0] or Random1>HyperbolicBounds[1]: FeasibleInput1.append(Random1)
					FeasinbleInput1 = np.array(FeasibleInput1)
					FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
											for el in FeasibleInput1])
				else:
					FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
					FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
											for el in FeasibleInput1])
			cline.set_xdata(FeasibleInput1)
			cline.set_ydata(FeasibleInput2)
			chosenpoint.set_xdata(u1[i])
			chosenpoint.set_ydata(u2[i])
			TimeText.set_text("t = " + str(t[i]))
			return cline,chosenpoint,TimeText,


		# Init only required for blitting to give a clean slate.
		def init():
			ax1.plot(Input1,lo1*lo2*0.001**2/Input1,'r',lw=2)
			Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[0],[x1[0],x2[0],x3[0],x4[0]])
			if abs(Coefficient1) <= 1e-7:
				LowerBound = Bounds[0][0]
				UpperBound = Bounds[0][1]
				if Constraint1/Coefficient2 > 0:
					LowerBound = Bounds[0][0]
					UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
				else:
					LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
					UpperBound = Bounds[0][1]
				FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
				FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
			elif abs(Coefficient2) <= 1e-7:
				LowerBound = Constraint1/Coefficient1
				UpperBound = Constraint1/Coefficient1
				FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
				if Constraint1/Coefficient1 < 0:
					LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
					UpperBound = Bounds[1][1]
				else:
					LowerBound = Bounds[1][0]
					UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
				FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			elif np.sign(-Coefficient1) == np.sign(Coefficient2): # DIFFERENT SIGNS
				if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 < Bounds[1][0]:
					LowerBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
				else:
					LowerBound = Bounds[0][0]

				if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 > Bounds[1][1]:
					UpperBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
				else:
					UpperBound = Bounds[0][1]
				assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
				HyperbolicBounds = np.sort([(Constraint1 - \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1), \
										 	(Constraint1 + \
												np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
													/(2*Coefficient1)])
				LowerBound = max([LowerBound,HyperbolicBounds[0]])
				UpperBound = min([UpperBound,HyperbolicBounds[1]])
				FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
				FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
										for el in FeasibleInput1])
			else: # np.sign(-Coefficient1) != np.sign(Coefficient2) SAME SIGN
				if (Constraint1 - Coefficient1*Bounds[0][0])/Coefficient2 > Bounds[1][1]:
					LowerBound = (Constraint1-Coefficient2*Bounds[1][1])/Coefficient1
				else:
					LowerBound = Bounds[0][0]

				if (Constraint1 - Coefficient1*Bounds[0][1])/Coefficient2 < Bounds[1][0]:
					UpperBound = (Constraint1-Coefficient2*Bounds[1][0])/Coefficient1
				else:
					UpperBound = Bounds[0][1]
				assert UpperBound>LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
				if Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001) > 0:
					HyperbolicBounds = np.sort([(Constraint1 - \
													np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
														/(2*Coefficient1), \
											 	(Constraint1 + \
													np.sqrt(Constraint1**2 - 4*Coefficient1*Coefficient2*(lo1*0.001)*(lo2*0.001)))\
														/(2*Coefficient1)])

					assert LowerBound < HyperbolicBounds[0] or HyperbolicBounds[1] < UpperBound, "No feasible solutions."

					FeasibleInput1 = []
					while len(FeasibleInput1)<1000:
						Random1 = (UpperBound-LowerBound)*np.random.rand() + LowerBound
						if Random1<HyperbolicBounds[0] or Random1>HyperbolicBounds[1]: FeasibleInput1.append(Random1)
					FeasinbleInput1 = np.array(FeasibleInput1)
					FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
											for el in FeasibleInput1])
				else:
					FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
					FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
											for el in FeasibleInput1])
			cline, = plt.plot(FeasibleInput1,FeasibleInput2,'b',lw=2)
			cline.set_visible(False)
			chosenpoint, = plt.plot(u1[0],u2[0],c='k',marker='o')
			chosenpoint.set_visible(False)
			TimeText = plt.text(0.75,0.75,"t = " + str(t[0]),fontsize=16)
			TimeText.set_visible(False)
			return cline,chosenpoint,TimeText,

		ani = animation.FuncAnimation(fig, animate, np.arange(1, len(t),1), init_func=init,interval=1, blit=False)
		plt.show()
def plot_individual_constraint_versus_time_test_2(t,x1,x2,x3,x4,Return = False):
	import numpy as np
	import matplotlib.pyplot as plt

	DescriptiveTitle = "Plotting Coefficients/Constraints vs. Time"
	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.975)

	"""
	A⋅u₁ + B⋅u₂ = C
	"""

	A,B,C = [],[],[]
	for i in range(len(x1)):
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[i],[x1[i],x2[i],x3[i],x4[i]])
		A.append(Coefficient1)
		B.append(Coefficient2)
		C.append(Constraint1)

	ax1.plot(t[:len(x1)],A,'r',lw=2)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$1^{st}$ Coefficient")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:len(x1)],B,'b',lw=2)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$2^{nd}$ Coefficient")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	ax3.plot(t[:len(x1)],C,'k',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel("Constraint")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_individual_coefficient2_versus_time_test_2(t,x1,x2,x3,x4,Return = False):
	import numpy as np
	import matplotlib.pyplot as plt

	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
	plt.subplots_adjust(top=0.9,hspace=0.4,bottom=0.1,left=0.075,right=0.975)
	plt.suptitle(r"Plotting $2^{nd}$ Coefficient vs. Time",Fontsize=20,y=0.975)

	"""
	B = c2⋅c4⋅R2(X)⋅KT_2(X)
	"""

	r2,kt_2,B = [],[],[]
	for i in range(len(x1)):
		_,Coefficient2,_ = return_constraint_variables_muscle_velocity_driven(t[i],[x1[i],x2[i],x3[i],x4[i]])
		B.append(Coefficient2)
		r2.append(R2([x1[i],x2[i],x3[i],x4[i]]))
		kt_2.append(G5([x1[i],x2[i],x3[i],x4[i]]))

	ax1.plot(t[:len(x1)],r2,'b--',lw=2)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$g_{2}(\vec{x}(t))$")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:len(x1)],kt_2,'b:',lw=2)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$g_{5}(\vec{x}(t))$")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	ax3.plot(t[:len(x1)],B,'b',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel(r"$2^{nd}$ Coefficient")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_individual_coefficient1_versus_time_test_2(t,x1,x2,x3,x4,Return = False):
	import numpy as np
	import matplotlib.pyplot as plt

	DescriptiveTitle = "Plotting Coefficients/Constraints vs. Time"
	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
	plt.subplots_adjust(top=0.9,hspace=0.4,bottom=0.1,left=0.075,right=0.975)
	plt.suptitle(r"Plotting $1^{st}$ Coefficient vs. Time",Fontsize=20,y=0.975)

	"""
	A = c2⋅c3⋅R1(X)⋅KT_1(X)
	"""

	r1,kt_1,B = [],[],[]
	for i in range(len(x1)):
		Coefficient1,_,_ = return_constraint_variables_muscle_velocity_driven(t[i],[x1[i],x2[i],x3[i],x4[i]])
		B.append(Coefficient1)
		r1.append(R1([x1[i],x2[i],x3[i],x4[i]]))
		kt_1.append(KT_1([x1[i],x2[i],x3[i],x4[i]]))

	ax1.plot(t[:len(x1)],r1,'r--',lw=2)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$g_{2}(\vec{x}(t))$")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:len(x1)],kt_1,'r:',lw=2)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$g_{5}(\vec{x}(t))$")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	ax3.plot(t[:len(x1)],B,'r',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel(r"$2^{nd}$ Coefficient")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_states(t,X,Return=False):
	import numpy as np
	import matplotlib.pyplot as plt

	NumStates = np.shape(X)[0]
	NumRows = int(np.ceil(NumStates/5))
	if NumStates < 5:
		NumColumns = NumStates
	else:
		NumColumns = 5

	ColumnNumber = [el%5 for el in np.arange(0,NumStates,1)]
	RowNumber = [int(el/5) for el in np.arange(0,NumStates,1)]

	DescriptiveTitle = "Plotting States vs. Time"
	fig, axes = plt.subplots(NumRows,NumColumns,figsize=(3*NumColumns,2*NumRows + 2))
	plt.subplots_adjust(top=0.85,hspace=0.4,bottom=0.15,left=0.075,right=0.975)
	plt.suptitle(r"Plotting $1^{st}$ Coefficient vs. Time",Fontsize=20,y=0.975)
	if NumStates <=5:
		for j in range(NumStates):
			axes[ColumnNumber[j]].spines['right'].set_visible(False)
			axes[ColumnNumber[j]].spines['top'].set_visible(False)
			axes[ColumnNumber[j]].plot(t[:np.shape(X)[1]],X[j])
			if ColumnNumber[j]!=0:
				axes[ColumnNumber[j]].set_xticklabels(\
									[""]*len(axes[ColumnNumber[j]].get_xticks()))
			else:
				axes[ColumnNumber[j]].set_xlabel("Time (s)")
			axes[ColumnNumber[j]].set_title(r"$x_{" + str(j+1) + "}$")

	else:
		for j in range(NumStates):
			axes[RowNumber[j],ColumnNumber[j]].spines['right'].set_visible(False)
			axes[RowNumber[j],ColumnNumber[j]].spines['top'].set_visible(False)
			axes[RowNumber[j],ColumnNumber[j]].plot(t[:np.shape(X)[1]],X[j])
			if not(RowNumber[j] == RowNumber[-1] and ColumnNumber[j]==0):
				axes[RowNumber[j],ColumnNumber[j]].set_xticklabels(\
									[""]*len(axes[RowNumber[j],ColumnNumber[j]].get_xticks()))
			else:
				axes[RowNumber[j],ColumnNumber[j]].set_xlabel("Time (s)")
			axes[RowNumber[j],ColumnNumber[j]].set_title(r"$x_{" + str(j+1) + "}$")
		if NumStates%5!=0:
			[fig.delaxes(axes[RowNumber[-1],el]) for el in range(ColumnNumber[-1]+1,5)]

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_inputs(t,u1,u2,Return=False):
	import numpy as np
	import matplotlib.pyplot as plt

	DescriptiveTitle = "Plotting Inputs vs. Time"
	fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,5))
	plt.subplots_adjust(top=0.9,hspace=0.4,bottom=0.1,left=0.075,right=0.975)
	plt.suptitle("Plotting Inputs vs. Time",Fontsize=20,y=0.975)

	ax1.plot(t[:len(u1)],u1,'g--',lw=2)
	ax1.plot([-1,t[len(u1)]+1],[0,0],'k--',lw=0.5)
	ax1.set_xlim([t[0],t[len(u1)]])
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$u_1$")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:len(u2)],u2,'g',lw=2)
	ax2.plot([-1,t[len(u2)]+1],[0,0],'k--',lw=0.5)
	ax2.set_xlim([t[0],t[len(u1)]])
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$u_2$")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
x1_1,x2_1 = [Base],[Amp*Freq]
u1_1,u2_1 = [100],[10]

x1_2,x2_2,x3_2,x4_2= [Base],[Amp*Freq],[100],[70]
u1_2,u2_2 = [0.01],[0.01]


CocontractionIndex = 2

AddNoise = False
if AddNoise == True:
    np.random.seed(seed=None)
    NoiseArray = np.random.normal(loc=0.0,scale=0.2,size=(2,len(Time)))
else:
    NoiseArray = np.zeros((2,len(Time)))

def update_policy_tension_driven(t,x1_1,x2_1,dt,NoiseArray):
	import numpy as np
	Method = "Tension"
	X = [x1_1[-1],x2_1[-1]]
	U = [u1_1[-1],u2_1[-1]]
	U = return_U_tension_driven(t,X,U,dt,MaxStep_Tension,Tension_Bounds,NoiseArray[:,int(t/dt)])
	u1_1.append(U[0])
	u2_1.append(U[1])
	x2_1.append(x2_1[-1] + dX2_dt(X,U=U)*dt)
	x1_1.append(x1_1[-1] + dX1_dt(X)*dt)
def update_policy_muscle_velocity_driven(t,x1_2,x2_2,x3_2,x4_2,dt,NoiseArray):
	import numpy as np
	Method = "Muscle Velocity"
	X = [x1_2[-1],x2_2[-1],x3_2[-1],x4_2[-1]]
	U = [u1_2[-1],u2_2[-1]]
	U = return_U_muscle_velocity_driven(t,X,U,dt,MaxStep_MuscleVelocity,MuscleVelocity_Bounds,NoiseArray[:,int(t/dt)])
	u1_2.append(U[0])
	u2_2.append(U[1])
	x4_2.append(x4_2[-1] + dX4_dt(X,U=U)*dt)
	x3_2.append(x3_2[-1] + dX3_dt(X,U=U)*dt)
	x2_2.append(x2_2[-1] + dX2_dt(X)*dt)
	x1_2.append(x1_2[-1] + dX1_dt(X)*dt)

StartTime = time.time()
for t in Time[1:]:
	update_policy_tension_driven(t,x1_1,x2_1,dt,NoiseArray)
	update_policy_muscle_velocity_driven(t,x1_2,x2_2,x3_2,x4_2,dt,NoiseArray)
	statusbar(int(t/dt),len(Time),StartTime=StartTime,Title="Forced-Pendulum")

fig1,[ax1_1,ax2_1,ax3_1,ax4_1] = plot_MA_values(Time,x1_1)
fig2,[ax1_2,ax2_2,ax3_2,ax4_2] = plot_MA_values(Time,x1_2)

plt.figure()
plt.title("Underdetermined Tendon-Tension-Driven\nForced Pendulum Example",\
                fontsize=16,color='gray')
plt.plot(Time,x1_1,'b',lw=2)
plt.plot(Time,x1_2,'g',lw=2)
plt.plot(Time,r(Time),'r--')
plt.xlabel("Time (s)")
plt.ylabel("Desired Measure")
plt.legend([r"Output $y = x_{1}$ (Tension)",r"Output $y = x_{1}$ (mm Velocity)",r"Reference $r(t) = \frac{\pi}{24}\sin(2\pi t) + \frac{\pi}{2}$"],loc='best')

plt.figure()
plt.title('Error vs. Time')
plt.plot(Time, r(Time)-x1_1,color='b')
plt.plot(Time, r(Time)-x1_2,color='g')
plt.legend(["Tension Driven","Muscle Velocity Driven"],loc='best')
plt.xlabel("Time (s)")
plt.ylabel("Error")

plt.figure()
plt.plot(Time,u1_1,'g',Time,u2_1,'r')
plt.title('Tendon Tensions vs. Time')
plt.xlabel("Time (s)")
plt.ylabel("Tendon Tensions (N)")
plt.legend(["Muscle 1","Muscle 2"])

plt.figure()
plt.plot(Time,u1_2,'g',Time,u2_2,'r')
plt.title('Muscle Velocities vs. Time')
plt.xlabel("Time (s)")
plt.ylabel("Muscle Velocities (m/s)")
plt.legend(["Muscle 1","Muscle 2"])

plt.show()
