import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import sympy as sy
from sympy.utilities import lambdify
import time
from collections import namedtuple
from scipy import integrate
import matplotlib._pylab_helpers

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
	from sympy.utilities import lambdify
	import numpy as np
	from numpy import pi

	# Coefficients from observation, Ramsay; 2009, Pigeon, FVC, Holtzbaur, Garner & Pandy, or Banks.

	MA_Settings = namedtuple("MA_Settings",["Values","Source","Units","Equation_Number","Threshold","DOF"])
	Spindle_Settings = namedtuple("Spindle_Settings",["ActualNumber",'CorrectedNumber','RelativeAbundance',"Source"])
	Input_Source = namedtuple("Source_Settings",["Values","Source","Units"])

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
		}

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
		}

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
		}

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
		}

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
		}

	PRO_Settings = {\
		'Notes' : [\
					""\
					],
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
		}

	FCR_Settings = {\
		'Notes' : [\
					"FCR EFE MA is not listed in Ramsay; 2009 but Pigeon has a quadratic function with R² = 0.9975. Pigeon only takes elbow angle into account. Coefficients and equation number/type are listed below to test either implementation. EFE MA was estimated to be constant and 10 mm for this muscle. If you use Pigeon, make sure to only accept positive moment arm values, as this model fails outside the ROM. One option is to set the MA to the constant value (i.e., MA[theta>140°] = MA[140°])"\
					],\
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
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(6.3, 'Holzbaur; 2005', 'cm'),\
											Input_Source(5.10, 'Garner & Pandy; 2003', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(24.4, "Holzbaur; 2005", 'cm'),\
											Input_Source(27.08, "Garner & Pandy; 2003", 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(3, "Holzbaur; 2005", 'deg')
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(1.6,"Holzbaur; 2005",'sq cm'),\
							Input_Source(11.16,"Garner & Pandy; 2003",'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(74.0,"Holzbaur; 2005","N"),\
												Input_Source(368.63,"Garner & Pandy; 2003","N")\
											]}\
		}

	ECRB_Settings = {\
		'Notes' : [\
					"Garner and Pandy do not distinguish between ERCB and ECRL. Parameter values were included for completeness."\
					],\
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Ramsay; 2009", \
	 					"Sources" : \
							[\
								MA_Settings([-11.256,17.8548,1.6398,-0.5073,-2.8827,0,-0.0942,0,0,0,0,0,0,0,0,0], 'Ramsay; 2009', 'mm', 2, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(102,132.7,0.77,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006",\
					"Sources" : \
						[\
							Input_Source(32.1, "Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(5.9, 'Holzbaur; 2005', 'cm'),\
											Input_Source(7.28, 'Garner & Pandy; 2003', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(22.2, "Holzbaur; 2005", 'cm'),\
											Input_Source(26.80, "Garner & Pandy; 2003", 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(9, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(2.2,"Holzbaur; 2005",'sq cm'),\
							Input_Source(24.89,"Garner & Pandy; 2003",'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(100.5,"Holzbaur; 2005","N"),\
												Input_Source(755.76,"Garner & Pandy; 2003","N")\
											]}\
		}

	ECRL_Settings = {\
		'Notes' : [\
					"ECRL EFE MA for Ramsay; 2009 has R² = 0.978 whereas Pigeon has R² = 0.9986. Pigeon, however, only takes elbow angle into account, whereas Ramsay; 2009 takes in variable PS angles. Additionally, Pigeon only considers elbow angles between 0 and 140 degrees and exhibits a decrease in MA as elbow angle approaches the upper bound of the ROM. This should (intiutively speaking) make the extensors MA largest, but Pigeon exhibits a drop off that may make it less favorable for movements at the boundary of the ROM. Coefficients and equation number/type are listed below to test either implementation.",\
					"Garner and Pandy do not distinguish between ERCB and ECRL. Parameter values were included for completeness. Might need to divide by 2 to 'split' muscle."\
					],
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Pigeon; 1996", \
	 					"Sources" : \
							[\
								MA_Settings(Pigeon_coeff_conversion([4.7304,1.2590,4.4347,-3.0229,0,0]), 'Pigeon; 1996', 'mm', None, None, 'Elbow'),\
								MA_Settings([-7.7034,16.3913,7.4361,-1.7566,0,-1.3336,0,0.0742,0,0,0,0,0,0,0,0], 'Ramsay; 2009', 'mm', 2, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(74,155.2,0.48,"Banks; 2006"),\
		'Mass' : { "Primary Source" : "Banks; 2006", \
					"Sources" : \
						[\
							Input_Source(44.3, "Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(8.1, 'Holzbaur; 2005', 'cm'),\
											Input_Source(7.28, 'Garner & Pandy; 2003', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(22.4, "Holzbaur; 2005", 'cm'),\
											Input_Source(26.80, "Garner & Pandy; 2003", 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(0, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(2.2,"Holzbaur; 2005",'sq cm'),\
							Input_Source(24.89,"Garner & Pandy; 2003",'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(304.9,"Holzbaur; 2005",'N'),\
												Input_Source(755.76,"Garner & Pandy; 2003",'N')\
											]}\
		}

	FCU_Settings = {\
		'Notes' : [\
					""\
					],
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Est", \
	 					"Sources" : \
							[\
								MA_Settings(1, 'Est', 'cm', None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(175,141.2,1.2,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006", \
					"Sources" : \
						[\
							Input_Source(36.5,"Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(5.1, 'Holzbaur; 2005', 'cm'),\
											Input_Source(3.98, 'Garner & Pandy; 2003', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(26.5, "Holzbaur; 2005", 'cm'),\
											Input_Source(27.14, "Garner & Pandy; 2003", 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(12, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(2.9,"Holzbaur; 2005",'sq cm'),\
							Input_Source(16.99,"Garner & Pandy; 2003",'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(128.9,"Holzbaur; 2005",'N'),\
												Input_Source(561.00,"Garner & Pandy; 2003",'N')\
											]}\
		}

	FDS_Settings = {\
		'Notes' : [\
					"Note: only they muscle for the second digit was used for the FDS muscle. NEED TO DETERMINE IF THIS SHOULD BE A SUM OR AN AVERAGE FOR MEASURES LIKE PCSA, F_MAX, ETC.",\
					"As we are only considering one digit, do we need to sum all of the peak forces in order to get a better representation of its force producing capabilities?"\
					],
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Est", \
	 					"Sources" : \
							[\
								MA_Settings(1, 'Est', 'cm', None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(356,224.9,1.6,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006",\
					"Sources" : \
						[\
							Input_Source(95.2,"Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(8.4, 'Holzbaur; 2005', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(27.5, "Holzbaur; 2005", 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(6, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(1.4,"Holzbaur; 2005",'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(61.2,"Holzbaur; 2005",'N')\
											]}\
		}

	PL_Settings = {\
		'Notes' : [\
					""\
					],
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Est", \
	 					"Sources" : \
							[\
								MA_Settings(1, 'Est', 'cm', None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(None,None,None,None),\
		'Mass' : {	"Primary Source" : None,\
					"Sources" : \
						[\
							Input_Source(None, "N/A","N/A")\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(6.4, 'Holzbaur; 2005', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(26.9, "Holzbaur; 2005", 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(4, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(0.6,"Holzbaur; 2005",'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(26.7,"Holzbaur; 2005",'N')\
											]}\
		}

	ECU_Settings = {\
		'Notes' : [\
					"ECU EFE MA is not listed in Ramsay; 2009 but Pigeon has a quadratic function with R² = 0.9966. Pigeon only takes elbow angle into account. Coefficients and equation number/type are listed below to test either implementation. EFE MA was estimated to be constant and -10 mm for this muscle. If you use Pigeon, make sure to only accept negative moment arm values, as this model fails outside the ROM. One option is to set the MA to the constant value (i.e., MA[theta>140°] = MA[140°])"\
					],
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Pigeon; 1996", \
	 					"Sources" : \
							[\
								MA_Settings(Pigeon_coeff_conversion([-2.1826,-1.7386,1.1491,0,0,0]), 'Pigeon; 1996', 'mm', None, None, 'Elbow'),\
								MA_Settings(-1, 'Est', 'cm', None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(157,118,1.3,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006",\
					"Sources" : \
						[\
							Input_Source(25.2,"Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(6.2, 'Holzbaur; 2005', 'cm'),\
											Input_Source(3.56, 'Garner & Pandy; 2003', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(26.5, "Holzbaur; 2005", 'cm'),\
											Input_Source(28.18, "Garner & Pandy; 2003", 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(12, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(2.9,"Holzbaur; 2005",'sq cm'),\
							Input_Source(8.04,"Garner & Pandy; 2003",'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(128.9,"Holzbaur; 2005",'N'),\
												Input_Source(265.58,"Garner & Pandy; 2003",'N')\
											]}\
		}

	EDM_Settings = {\
		'Notes' : [\
					""\
					],
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Est", \
	 					"Sources" : \
							[\
								MA_Settings(-1, 'Est', 'cm', None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(53,59.8,0.89,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006", \
					"Sources" : \
						[\
							Input_Source(6.2, "Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(6.8, 'Holzbaur; 2005', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : 'Holzbaur; 2005', \
									"Sources" : \
										[\
											Input_Source(32.2, 'Holzbaur; 2005', 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : 'Holzbaur; 2005', \
								"Sources" : \
									[\
										Input_Source(3, 'Holzbaur; 2005', 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : 'Holzbaur; 2005', \
					'Sources' : \
						[\
							Input_Source(0.6,'Holzbaur; 2005','sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : 'Holzbaur; 2005', \
										'Sources': \
											[\
												Input_Source(25.3,'Holzbaur; 2005','N')\
											]}\
		}

	EDC_Settings = {\
		'Notes' : [\
					"Note: only they muscle for the second digit was used for the EDC muscle. NEED TO DETERMINE IF THIS SHOULD BE A SUM OR AN AVERAGE FOR MEASURES LIKE PCSA, F_MAX, ETC."\
					],
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Est", \
	 					"Sources" : \
							[\
								MA_Settings(-1, 'Est', 'cm', None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(219,152.6,1.4,"Banks; 2006"),\
		'Mass' : {	"Primary Source" : "Banks; 2006", \
					"Sources" : \
						[\
							Input_Source(42.8, "Banks; 2006", 'g')\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(7.0, 'Holzbaur; 2005', 'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(32.2, "Holzbaur; 2005", 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(3, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(18.3,"Holzbaur; 2005",'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(0.4,"Holzbaur; 2005",'N')\
											]}\
		}

	AN_Settings = {\
		'Notes' : [\
					""\
					],
		'Shoulder MA' : {	"Primary Source" : "Est",\
		 					"Sources" : \
								[\
									MA_Settings(0, 'Est', 'm', None, None, 'Shoulder')\
								]}, \
		'Elbow MA' : {	"Primary Source" : "Pigeon; 1996", \
	 					"Sources" : \
							[\
								MA_Settings(Pigeon_coeff_conversion([-5.3450,-2.2841,8.4297,-14.329,10.448,-2.736]), 'Pigeon; 1996', 'mm', None, None, 'Elbow')\
							]}, \
		'Spindle' : Spindle_Settings(None,None,None,None),\
		'Mass' : {	"Primary Source" : None, \
					"Sources" : \
						[\
							Input_Source(None, None, None)\
						]},\
		'Optimal Muscle Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(2.7,"Holzbaur; 2005",'cm')\
										]},\
		'Optimal Tendon Length' : {	"Primary Source" : "Holzbaur; 2005", \
									"Sources" : \
										[\
											Input_Source(1.8, "Holzbaur; 2005", 'cm')\
										]}, \
		'Pennation Angle' : {	"Primary Source" : "Holzbaur; 2005", \
								"Sources" : \
									[\
										Input_Source(0, "Holzbaur; 2005", 'deg')\
									]}, \
		'PCSA' : {	"Primary Source" : "Holzbaur; 2005", \
					'Sources' : \
						[\
							Input_Source(2.5,"Holzbaur; 2005",'sq cm')\
						]}, \
		'Maximum Isometric Force': {	"Primary Source" : "Holzbaur; 2005", \
										'Sources': \
											[\
												Input_Source(350.0,"Holzbaur; 2005",'N')\
											]}\
		}

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
	if PreselectedMuscles is None:
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
def unit_conversion(Params):
	import numpy as np
	Units = Params.Units
	Value = Params.Values
	assert Units.capitalize() in ["Degrees","Deg","Degree","Radians","Radian","Rad","Rads","In","Inches","Cm","Centimeters","Centimeter","Mm","Millimeters","Millimeter","Meters","Meter","M","Sq in","Squared inches","Inches squared","In sq","Cm sq","Sq cm","Centimeters squared","Squared centimeters","Mm sq","Sq mm","Millimeters squared","Squared millimeters","Meters squared","Squared meter","M sq","Sq m","Lbs","Lb","Pounds","G","Grams","Gram","Kg","Kilograms","Kilogram","N","Newton","Newtons"], "Improper Units Value. Please use appropriate Units."

	if Units.capitalize() in ["Degrees","Deg","Degree","Radians","Radian","Rad","Rads"]:
		if Units.capitalize() in ["Radians","Radian","Rad","Rads"]:
			return(Value)
		elif Units.capitalize() in ["Degrees","Deg","Degree"]:
			if type(Value)==list:
				return(list(np.array(Value)*np.pi/180))
			else:
				return(Value*np.pi/180)

	elif Units.capitalize() in ["In","Inches","Cm","Centimeters","Centimeter","Mm",\
					"Millimeters","Millimeter","Meters","Meter","M"]:
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

	elif Units.capitalize() in ["Sq in","Squared inches","Inches squared","In sq","Cm sq","Sq cm",\
					"Centimeters squared","Squared centimeters","Mm sq","Sq mm",\
					"Millimeters squared","Squared millimeters","Meters squared",\
					"Squared meter","M sq","Sq m"]:
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

	elif Units.capitalize() in ["Lbs","Lb","Pounds","G","Grams","Gram","Kg","Kilograms","Kilogram"]:
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

	elif Units.capitalize() in ["N","Newton","Newtons"]:
		if Units.capitalize() in ["N","Newton","Newtons"]:
			return(Value)
def return_optimal_length(MuscleSettings):
	import numpy as np
	OptimalMuscleLengthsList = MuscleSettings["Optimal Muscle Length"]["Sources"]
	PrimarySource = MuscleSettings["Optimal Muscle Length"]["Primary Source"]
	return(OptimalMuscleLengthsList[np.where([src[1] == PrimarySource for src in OptimalMuscleLengthsList])[0][0]])
def MA_function(Parameters,θ_PS=None):
	"""
	Note:

	Angles should be a number if Coefficients has a length of 5, or a list of length 2 when the Coefficients have lengths 16 or 18. Angles[0] will be the PRIMARY ANGLE for the DOF being considered while Angles[1] will be the secondary angle.

	Notes:

	threshold is only needed for Pigeon or Ramsay; 2009 MA functions that are invalid outside of a given value. Must be either None (default) or the radian value of the threshold.

	eq is only needed for Ramsay; 2009 (Pigeon has one quintic polynomial). eq must be either 1, 2, or 3, with list length requirements of 5, 16, or 18, respectively.
	"""
	import numpy as np
	Parameters = return_primary_source(Parameters)
	assert str(type(Parameters))=="<class '__main__.MA_Settings'>", "Parameters are not in correct namedtuple form."
	if θ_PS is None:
		θ_PS = np.pi
	else:
		assert type(θ_PS)==float, "θ_PS must be a float."
	src = Parameters.Source
	Coefficients = unit_conversion(Parameters)
	eq = Parameters.Equation_Number
	threshold = Parameters.Threshold

	assert type(src) == str, "src must be a str."
	assert src.capitalize() in ['Ramsay; 2009','Pigeon; 1996','Kuechle; 1997','Holzbaur; 2005', 'Est'], "src must be either Ramsay; 2009, Pigeon or Est (Estimate)."

	'''
	Note:
	For Kuechle and Holzbaur, where estimates or average MA were given, the format should be [MA,0,0,0,0,0] such that the function returns a constant MA function (See matrix multiplication below).
	'''

	if src.capitalize() in ['Pigeon; 1996', 'Kuechle; 1997', 'Holzbaur; 2005']:
		assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
		MomentArm = lambda θ: (np.matrix(Coefficients,dtype='float64')\
										*np.matrix([1,θ,θ**2,θ**3,θ**4,θ**5]).T)[0,0]
	elif src.capitalize() == 'Est':
		MomentArm = lambda θ: np.array(Coefficients,dtype='float64')
	else: #src.capitalize() == 'Ramsay; 2009'
		assert type(Coefficients) == list, "Coefficients must be a list."
		assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
		assert eq in [1,2,3], "eq must be either 1, 2, or 3 when using Ramsay; 2009 (2009)."
		if eq == 1:
			assert len(Coefficients) == 5, "For Eq. 1, Coefficients must be 5 elements long."
			MomentArm = lambda θ: \
					(np.matrix(Coefficients,dtype='float64')*\
						np.matrix([1,θ,θ**2,θ**3,θ**4]).T)[0,0]
		elif eq == 2:
			assert len(Coefficients)==16, "For Eq. 2, Coefficients must be 16 elements long."
			MomentArm = lambda θ: \
					(np.matrix(Coefficients,dtype='float64')*\
						np.matrix([1, θ, θ_PS, θ*θ_PS, θ**2, \
									θ_PS**2, (θ**2)*θ_PS, θ*(θ_PS**2), \
									(θ**2)*(θ_PS**2), θ**3, θ_PS**3, \
									(θ**3)*θ_PS, θ*(θ_PS**3), \
									(θ**3)*(θ_PS**2), (θ**2)*(θ_PS**3), \
									(θ**3)*(θ_PS**3)]).T)[0, 0]
		else: # eq == 3
			assert len(Coefficients)==18, "For Eq. 3, Coefficients must be 18 elements long."
			MomentArm = lambda θ: \
					(np.matrix(Coefficients,dtype='float64')*\
						np.matrix([1, θ, θ_PS, θ*θ_PS, θ**2, \
									θ_PS**2, (θ**2)*θ_PS, θ*(θ_PS**2), (θ**2)*(θ_PS**2), \
									θ**3, (θ**3)*θ_PS, (θ**3)*(θ_PS**2), \
									θ**4, (θ**4)*θ_PS, (θ**4)*(θ_PS**2),  \
									θ**5, (θ**5)*θ_PS, (θ**5)*(θ_PS**2)]).T)[0, 0]
	if threshold is None:
		return(MomentArm)
	else:
		assert type(threshold) in [int,float], "threshold must be a number."
		PiecewiseMomentArm = lambda θ:\
					np.piecewise(θ,[θ<threshold,θ>=threshold],\
									[MomentArm(θ),MomentArm(threshold)])
		return(PiecewiseMomentArm)
def MA_deriv(Parameters,θ_PS=None):
	"""
	Note:

	Angles should be a number if Coefficients has a length of 5, or a list of length 2 when the Coefficients have lengths 16 or 18. Angles[0] will be the PRIMARY ANGLE for the DOF being considered while Angles[1] will be the secondary angle.

	Notes:

	threshold is only needed for Pigeon or Ramsay; 2009 MA functions that are invalid outside of a given value. Must be either None (default) or the radian value of the threshold.

	eq is only needed for Ramsay; 2009 (Pigeon has one quintic polynomial). eq must be either 1, 2, or 3, with list length requirements of 5, 16, or 18, respectively.
	"""
	import numpy as np
	Parameters = return_primary_source(Parameters)
	assert str(type(Parameters))=="<class '__main__.MA_Settings'>", "Parameters are not in correct namedtuple form."
	if θ_PS is None:
		θ_PS = np.pi
	else:
		assert type(θ_PS)==float, "θ_PS must be a float."
	src = Parameters.Source
	Coefficients = unit_conversion(Parameters)
	eq = Parameters.Equation_Number
	threshold = Parameters.Threshold

	assert type(src) == str, "src must be a str."
	assert src.capitalize() in ['Ramsay; 2009','Pigeon; 1996','Kuechle; 1997','Holzbaur; 2005', 'Est'], "src must be either Ramsay; 2009, Pigeon or Est (Estimate)."

	'''
	Note:
	For Kuechle and Holzbaur, where estimates or average MA were given, the format should be [MA,0,0,0,0,0] such that the function returns a constant MA function (See matrix multiplication below).
	'''

	if src.capitalize() in ['Pigeon; 1996', 'Kuechle; 1997', 'Holzbaur; 2005']:
		assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
		"""
		(d/dθ) [MomentArm] = (np.matrix(Coefficients,dtype='float64')\
						*np.matrix([(d/dθ)[1],(d/dθ)[θ],(d/dθ)[θ**2],\
										(d/dθ)[θ**3],(d/dθ)[θ**4],(d/dθ)[θ**5]]).T)[0,0]
		"""
		Derivative = lambda θ: (np.matrix(Coefficients,dtype='float64')\
						*np.matrix([0,1,(2*θ),(3*θ**2),(4*θ**3),(5*θ**4)]).T)[0,0]
	elif src.capitalize() == 'Est':
		"""
		(d/dθ)[MomentArm] = np.array((d/dθ)[Coefficients],dtype='float64')
		"""
		Derivative = lambda θ:  0
	else: #src.capitalize() == 'Ramsay; 2009'
		assert type(Coefficients) == list, "Coefficients must be a list."
		assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
		assert eq in [1,2,3], "eq must be either 1, 2, or 3 when using Ramsay; 2009 (2009)."
		if eq == 1:
			assert len(Coefficients) == 5, "For Eq. 1, Coefficients must be 5 elements long."
			"""
			(d/dθ)[MomentArm] = (np.matrix(Coefficients,dtype='float64')\
									*np.matrix([	(d/dθ)[1],		(d/dθ)[θ],\
													(d/dθ)[θ**2],	(d/dθ)[θ**3],\
													(d/dθ)[θ**4]]).T)[0,0]
			"""
			Derivative = lambda θ: (np.matrix(Coefficients,dtype='float64')\
							*np.matrix([	0,				1,\
											(2*θ),			(3*θ**2),\
											(4*θ**3)]).T)[0,0]
		elif eq == 2:
			"""
			This is only good for this ReferenceTracking Ex where PS is fixed. Derivative is only wrt one DOF.
			"""
			assert len(Coefficients)==16, "For Eq. 2, Coefficients must be 16 elements long."
			"""
			(d/dθ)[MomentArm] = \
					(np.matrix(Coefficients,dtype='float64')*\
							np.matrix([(d/dθ)[1], 					(d/dθ)[θ], \
										(d/dθ)[θ_PS], 				(d/dθ)[θ*θ_PS],\
										(d/dθ)[θ**2],				(d/dθ)[θ_PS**2],\
										(d/dθ)[(θ**2)*θ_PS],		(d/dθ)[θ*(θ_PS**2)],\
										(d/dθ)[(θ**2)*(θ_PS**2)],	(d/dθ)[θ**3], \
										(d/dθ)[θ_PS**3], 			(d/dθ)[(θ**3)*θ_PS],\
										(d/dθ)[θ*(θ_PS**3)],		(d/dθ)[(θ**3)*(θ_PS**2)],\
										(d/dθ)[(θ**2)*(θ_PS**3)],	(d/dθ)[(θ**3)*(θ_PS**3)]\
										]).T)[0, 0]
			"""
			Derivative = lambda θ: (np.matrix(Coefficients,dtype='float64')*\
							np.matrix([ 0, 					1,\
										0,  				θ_PS, \
										(2*θ), 				0, \
										(2*θ)*θ_PS, 		(θ_PS**2), 	\
										(2*θ)*(θ_PS**2), 	(3*θ**2), \
										0, 					(3*θ**2)*θ_PS, \
										(θ_PS**3), 			(3*θ**2)*(θ_PS**2), \
										(2*θ)*(θ_PS**3),	(3*θ**2)*(θ_PS**3)]).T)[0, 0]
		else: # eq == 3
			assert len(Coefficients)==18, "For Eq. 3, Coefficients must be 18 elements long."
			"""
			(d/dθ)[MomentArm] = \
					(np.matrix(Coefficients,dtype='float64')*\
							np.matrix([(d/dθ)[1], 					(d/dθ)[θ], \
										(d/dθ)[θ_PS], 				(d/dθ)[θ*θ_PS],\
										(d/dθ)[θ**2], 				(d/dθ)[θ_PS**2],\
										(d/dθ)[(θ**2)*θ_PS], 		(d/dθ)[θ*(θ_PS**2)],\
										(d/dθ)[(θ**2)*(θ_PS**2)], 	(d/dθ)[θ**3],\
										(d/dθ)[(θ**3)*θ_PS], 		(d/dθ)[(θ**3)*(θ_PS**2)],\
										(d/dθ)[θ**4], 				(d/dθ)[(θ**4)*θ_PS],\
										(d/dθ)[(θ**4)*(θ_PS**2)],  	(d/dθ)[θ**5],\
										(d/dθ)[(θ**5)*θ_PS], 		(d/dθ)[(θ**5)*(θ_PS**2)\
										]]).T)[0, 0]
			"""
			Derivative = lambda θ: (np.matrix(Coefficients,dtype='float64')*\
								np.matrix([	0, 					1,\
											0, 					θ_PS,\
										  	(2*θ),				0,\
									 		(2*θ)*θ_PS, 		(θ_PS**2),\
									  		(2*θ)*(θ_PS**2),	(3*θ**2),\
									 		(3*θ**2)*θ_PS, 		(3*θ**2)*(θ_PS**2),\
											(4*θ**3), 			(4*θ**3)*θ_PS,\
									 		(4*θ**3)*(θ_PS**2),	(5*θ**4),\
									 		(5*θ**4)*θ_PS, 		(5*θ**4)*(θ_PS**2)]).T)[0, 0]
	if threshold is None:
		return(Derivative)
	else:
		assert type(threshold) in [int,float], "threshold must be a number."
		PiecewiseDerivative = lambda θ:\
									np.piecewise(θ,[θ<threshold,θ>=threshold],\
													[Derivative(θ),0])
		return(PiecewiseDerivative)
def MA_2nd_deriv(Parameters,θ_PS=None):
	"""
	Note:

	Angles should be a number if Coefficients has a length of 5, or a list of length 2 when the Coefficients have lengths 16 or 18. Angles[0] will be the PRIMARY ANGLE for the DOF being considered while Angles[1] will be the secondary angle.

	Notes:

	threshold is only needed for Pigeon or Ramsay; 2009 MA functions that are invalid outside of a given value. Must be either None (default) or the radian value of the threshold.

	eq is only needed for Ramsay; 2009 (Pigeon has one quintic polynomial). eq must be either 1, 2, or 3, with list length requirements of 5, 16, or 18, respectively.
	"""

	import numpy as np
	Parameters = return_primary_source(Parameters)
	assert str(type(Parameters))=="<class '__main__.MA_Settings'>", "Parameters are not in correct namedtuple form."
	if θ_PS is None:
		θ_PS = np.pi
	else:
		assert type(θ_PS)==float, "θ_PS must be a float."
	src = Parameters.Source
	Coefficients = unit_conversion(Parameters)
	eq = Parameters.Equation_Number
	threshold = Parameters.Threshold

	assert type(src) == str, "src must be a str."
	assert src.capitalize() in ['Ramsay; 2009','Pigeon; 1996','Kuechle; 1997','Holzbaur; 2005', 'Est'], "src must be either Ramsay; 2009, Pigeon or Est (Estimate)."

	'''
	Note:
	For Kuechle and Holzbaur, where estimates or average MA were given, the format should be [MA,0,0,0,0,0] such that the function returns a constant MA function (See matrix multiplication below).
	'''

	if src.capitalize() in ['Pigeon; 1996', 'Kuechle; 1997', 'Holzbaur; 2005']:
		assert len(Coefficients)==6, 'For Pigeon (1996) the list of Coefficients must be 6 elements long. Insert zeros (0) for any additional empty coefficients.'
		"""
		(d²/dθ²) [MomentArm] \
			= (np.matrix(Coefficients,dtype='float64')\
						*np.matrix([(d²/dθ²)[1],		(d²/dθ²)[θ],\
									(d²/dθ²)[θ**2],		(d²/dθ²)[θ**3],\
									(d²/dθ²)[θ**4],		(d²/dθ²)[θ**5]]).T)[0,0]

			= (d/dθ)[Derivative]
			= (np.matrix(Coefficients,dtype='float64')\
						*np.matrix([(d/dθ)[0],			(d/dθ)[1],\
									(d/dθ)[(2*θ)],		(d/dθ)[(3*θ**2)],\
									(d/dθ)[(4*θ**3)],	(d/dθ)[(5*θ**4)]]).T)[0,0]
		"""

		SecondDerivative = lambda θ: (np.matrix(Coefficients,dtype='float64')\
											*np.matrix([0,			0,\
														2,			6*θ,\
														(12*θ**2),	(20*θ**3)]).T)[0,0]
	elif src.capitalize() == 'Est':
		"""
		(d²/dθ²)[MomentArm] = np.array((d²/dθ²)[Coefficients],dtype='float64')
		"""
		# (d/dθ)[Derivative] = lambda θ:  (d/dθ)[0]
		SecondDerivative = lambda θ:  0
	else: #src.capitalize() == 'Ramsay; 2009'
		assert type(Coefficients) == list, "Coefficients must be a list."
		assert len(Coefficients) in [5,16,18], "Coefficients as a list must be of length 5, 16, or 18."
		assert eq in [1,2,3], "eq must be either 1, 2, or 3 when using Ramsay; 2009 (2009)."
		if eq == 1:
			assert len(Coefficients) == 5, "For Eq. 1, Coefficients must be 5 elements long."
			"""
			(d²/dθ²)[MomentArm] \
				= (np.matrix(Coefficients,dtype='float64')\
							*np.matrix([(d²/dθ²)[1],		(d²/dθ²)[θ],\
										(d²/dθ²)[θ**2],		(d²/dθ²)[θ**3],\
										(d²/dθ²)[θ**4]]).T)[0,0]

				= (d/dθ)[Derivative] \
				= (np.matrix(Coefficients,dtype='float64')\
							*np.matrix([(d/dθ)[0],			(d/dθ)[1],\
										(d/dθ)[(2*θ)],  	(d/dθ)[(3*θ**2)],\
										(d/dθ)[(4*θ**3)]]).T)[0,0]
			"""
			SecondDerivative = lambda θ: \
						(np.matrix(Coefficients,dtype='float64')\
								*np.matrix([0,			0,\
											2, 			6*θ,\
											(12*θ**2)				]).T)[0,0]
		elif eq == 2:
			"""
			This is only good for this ReferenceTracking Ex where PS is fixed. Derivative is only wrt one DOF.
			"""
			assert len(Coefficients)==16, "For Eq. 2, Coefficients must be 16 elements long."
			"""
			(d²/dθ²)[MomentArm] \
				= (np.matrix(Coefficients,dtype='float64')*\
						np.matrix([	(d²/dθ²)[1], 				(d²/dθ²)[θ], \
									(d²/dθ²)[θ_PS], 			(d²/dθ²)[θ*θ_PS],\
									(d²/dθ²)[θ**2],				(d²/dθ²)[θ_PS**2],\
									(d²/dθ²)[(θ**2)*θ_PS],		(d²/dθ²)[θ*(θ_PS**2)],\
									(d²/dθ²)[(θ**2)*(θ_PS**2)],	(d²/dθ²)[θ**3], \
									(d²/dθ²)[θ_PS**3], 			(d²/dθ²)[(θ**3)*θ_PS],\
									(d²/dθ²)[θ*(θ_PS**3)],		(d²/dθ²)[(θ**3)*(θ_PS**2)],\
									(d²/dθ²)[(θ**2)*(θ_PS**3)],	(d²/dθ²)[(θ**3)*(θ_PS**3)]\
									]).T)[0, 0]

				= (d/dθ)[Derivative] \
				= (np.matrix(Coefficients,dtype='float64')*\
						np.matrix([ (d/dθ)[0], 					(d/dθ)[1], \
									(d/dθ)[0], 					(d/dθ)[θ_PS], \
									(d/dθ)[(2*θ)], 				(d/dθ)[0], \
									(d/dθ)[(2*θ)*θ_PS], 		(d/dθ)[(θ_PS**2)], \
									(d/dθ)[(2*θ)*(θ_PS**2)], 	(d/dθ)[(3*θ**2)], \
									(d/dθ)[0], 					(d/dθ)[(3*θ**2)*θ_PS], \
									(d/dθ)[(θ_PS**3)], 			(d/dθ)[(3*θ**2)*(θ_PS**2)], \
									(d/dθ)[(2*θ)*(θ_PS**3)],	(d/dθ)[(3*θ**2)*(θ_PS**3)]\
									]).T)[0, 0]
			"""
			SecondDerivative = lambda θ: \
					(np.matrix(Coefficients,dtype='float64')*\
						np.matrix([ 0, 					0,\
						 			0,					0,\
									2,					0,\
				 					2*θ_PS,				0,\
						 			2*(θ_PS**2),		(6*θ),\
									0,					(6*θ)*θ_PS,\
									0,					(6*θ)*(θ_PS**2),\
									2*(θ_PS**3),		(6*θ)*(θ_PS**3)\
									]).T)[0, 0]
		else: # eq == 3
			assert len(Coefficients)==18, "For Eq. 3, Coefficients must be 18 elements long."
			"""
			(d²/dθ²)[MomentArm] \
				= (np.matrix(Coefficients,dtype='float64')*\
						np.matrix([(d²/dθ²)[1], 					(d²/dθ²)[θ], \
									(d²/dθ²)[θ_PS], 				(d²/dθ²)[θ*θ_PS],\
									(d²/dθ²)[θ**2], 				(d²/dθ²)[θ_PS**2],\
									(d²/dθ²)[(θ**2)*θ_PS], 			(d²/dθ²)[θ*(θ_PS**2)],\
									(d²/dθ²)[(θ**2)*(θ_PS**2)], 	(d²/dθ²)[θ**3],\
									(d²/dθ²)[(θ**3)*θ_PS], 			(d²/dθ²)[(θ**3)*(θ_PS**2)],\
									(d²/dθ²)[θ**4], 				(d²/dθ²)[(θ**4)*θ_PS],\
									(d²/dθ²)[(θ**4)*(θ_PS**2)], 	(d²/dθ²)[θ**5],\
									(d²/dθ²)[(θ**5)*θ_PS], 			(d²/dθ²)[(θ**5)*(θ_PS**2)\
									]]).T)[0, 0]

				= (d/dθ)[Derivative] \
				= (np.matrix(Coefficients,dtype='float64')*\
						np.matrix([	(d/dθ)[0], 					(d/dθ)[1],\
									(d/dθ)[0], 					(d/dθ)[θ_PS],\
								  	(d/dθ)[(2*θ)],				(d/dθ)[0],\
							 		(d/dθ)[(2*θ)*θ_PS], 		(d/dθ)[(θ_PS**2)],\
							  		(d/dθ)[(2*θ)*(θ_PS**2)],	(d/dθ)[(3*θ**2)],\
							 		(d/dθ)[(3*θ**2)*θ_PS], 		(d/dθ)[(3*θ**2)*(θ_PS**2)],\
									(d/dθ)[(4*θ**3)], 			(d/dθ)[(4*θ**3)*θ_PS],\
							 		(d/dθ)[(4*θ**3)*(θ_PS**2)],	(d/dθ)[(5*θ**4)],\
							 		(d/dθ)[(5*θ**4)*θ_PS], 		(d/dθ)[(5*θ**4)*(θ_PS**2)]\
									]).T)[0, 0]
			"""
			SecondDerivative = lambda θ: \
					(np.matrix(Coefficients,dtype='float64')*\
							np.matrix([	0, 						0,\
										0, 						0,\
										2,						0,\
										2*θ_PS, 				0,\
										2*(θ_PS**2),			(6*θ),\
										(6*θ)*θ_PS, 			(6*θ)*(θ_PS**2),\
										(12*θ**2), 				(12*θ**2)*θ_PS,\
										(12*θ**2)*(θ_PS**2),	(20*θ**3),\
										(20*θ**3)*θ_PS, 		(20*θ**3)*(θ_PS**2)\
										]).T)[0, 0]
	if threshold is None:
		return(SecondDerivative)
	else:
		assert type(threshold) in [int,float], "threshold must be a number."
		PiecewiseSecondDerivative = lambda θ:\
									np.piecewise(θ,[θ<threshold,θ>=threshold],\
													[SecondDerivative(θ),0])
		return(PiecewiseSecondDerivative)
def return_MA_matrix_functions(AllMuscleSettings,ReturnMatrixFunction=False,θ_PS=None):
	import numpy as np
	import sympy as sp
	from sympy.utilities import lambdify

	MuscleList = AllMuscleSettings.keys()
	if ReturnMatrixFunction == False:
		R_Transpose = np.matrix([\
			[MA_function(AllMuscleSettings[muscle]["Shoulder MA"],θ_PS=θ_PS), \
			MA_function(AllMuscleSettings[muscle]["Elbow MA"],θ_PS=θ_PS)]\
				for muscle in MuscleList])
		dR_Transpose = np.matrix([\
			[MA_deriv(AllMuscleSettings[muscle]["Shoulder MA"],θ_PS=θ_PS), \
			MA_deriv(AllMuscleSettings[muscle]["Elbow MA"],θ_PS=θ_PS)]\
				for muscle in MuscleList])
		d2R_Transpose = np.matrix([\
			[MA_2nd_deriv(AllMuscleSettings[muscle]["Shoulder MA"],θ_PS=θ_PS), \
			MA_2nd_deriv(AllMuscleSettings[muscle]["Elbow MA"],θ_PS=θ_PS)]\
					for muscle in MuscleList])
	else:
		R_Transpose = lambda θ_SFE, θ_EFE: \
			np.matrix([[MA_function(AllMuscleSettings[muscle]["Shoulder MA"],θ_PS=θ_PS)(θ_SFE),\
						MA_function(AllMuscleSettings[muscle]["Elbow MA"],θ_PS=θ_PS)(θ_EFE)] \
							for muscle in MuscleList])
		dR_Transpose = lambda θ_SFE, θ_EFE, θ_PS: \
			np.matrix([[MA_deriv(AllMuscleSettings[muscle]["Shoulder MA"],θ_PS=θ_PS)(θ_SFE),\
						MA_deriv(AllMuscleSettings[muscle]["Elbow MA"],θ_PS=θ_PS)(θ_EFE)] \
							for muscle in MuscleList])
		d2R_Transpose = lambda θ_SFE, θ_EFE, θ_PS: \
			np.matrix([[MA_2nd_deriv(AllMuscleSettings[muscle]["Shoulder MA"],θ_PS=θ_PS)(θ_SFE),\
						MA_2nd_deriv(AllMuscleSettings[muscle]["Elbow MA"],θ_PS=θ_PS)(θ_EFE)] \
							for muscle in MuscleList])
	# returns an (n,m) matrix when n is the number of muscles and m is the number of DOFS. We chose to return R.T because this is commonly utilized in muscle velocity calculations.
	return(R_Transpose,dR_Transpose,d2R_Transpose)
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

"""

################################
#### Muscle Velocity Driven ####
################################

x_1 &= \theta \\
x_2 &= \dot{\theta} \\
x_3 &= T_{1} \\
x_4 &= T_{2} \\
u_1 &= l_{m,1} \\
u_2 &= l_{m,2} \\

"""
N_seconds = 2

AllMuscleSettings = return_muscle_settings(PreselectedMuscles=[5,6])

# g,L = 9.80, 0.45 #m/s², m
g,L = 0, 0.45 #m/s², m REMOVING GRAVITY
M = 1.6 # kg

α1 = unit_conversion(return_primary_source(AllMuscleSettings["BIC"]["Pennation Angle"])) # rads
α2 = unit_conversion(return_primary_source(AllMuscleSettings["TRI"]["Pennation Angle"])) # rads

"""
There was some debate regarding damping terms. No explicit value was given. Loeb simplifies the equation by utilizing the F_{PE_{2,i}} damping term (η) instead, as this is added to B_{m,i}*(v_{m,i}/l_{o,i}) anyways. This η is very small (0.01), so the damping is not significant. Might need to do sensitivity analysis on these two values (currently set to zero) (06/16/2018).
"""

cT = 27.8
kT = 0.0047
# LrT = 0.964 # Not used in this formulation because of the created ODE.

lo1 = unit_conversion(return_primary_source(AllMuscleSettings["BIC"]["Optimal Muscle Length"]))
lo2 = unit_conversion(return_primary_source(AllMuscleSettings["TRI"]["Optimal Muscle Length"]))

lTo1 = unit_conversion(return_primary_source(AllMuscleSettings["BIC"]["Optimal Tendon Length"]))
lTo2 = unit_conversion(return_primary_source(AllMuscleSettings["TRI"]["Optimal Tendon Length"]))

R_Transpose, dR_Transpose, _ = \
			return_MA_matrix_functions(AllMuscleSettings,ReturnMatrixFunction=False,θ_PS=np.pi/2)
"""
R_Transpose, dR_Transpose, and d2R_Transpose are of the form (n,m), where n is the number of muscles and m in the number of joints. In order to unpack the two muscles used in this model, we first must get the elbow MA functions R_Transpose[:,1], then change to a 1xn matrix (by the transpose), and then change to an array to reduce the ndmin from 2 to 1.
"""

r1,r2 = np.array(R_Transpose[:,1].T)[0]
dr1_dθ, dr2_dθ = np.array(dR_Transpose[:,1].T)[0]

PCSA1 = unit_conversion(return_primary_source(AllMuscleSettings["BIC"]["PCSA"]))
PCSA2 = unit_conversion(return_primary_source(AllMuscleSettings["TRI"]["PCSA"]))

F_MAX1 = unit_conversion(return_primary_source(AllMuscleSettings["BIC"]["Maximum Isometric Force"]))
F_MAX2 = unit_conversion(return_primary_source(AllMuscleSettings["TRI"]["Maximum Isometric Force"]))

Amp = 7.5*np.pi/180
Base = 90*np.pi/180
Freq = 2*np.pi

k1,k2,k3 = 100,100,100

if g == 0:
	MaxStep_Tension = 0.03 # percentage of positive maximum.
	Tension_Bounds = [[0,F_MAX1],[0,F_MAX2]]

	MaxStep_MuscleVelocity = 0.05 # percentage of positive maximum.
	MuscleVelocity_Bounds =[[-5*lo1,5*lo1],[-5*lo2,5*lo2]]
else:
	MaxStep_Tension = 0.01 # percentage of positive maximum.
	Tension_Bounds = [[0,F_MAX1],[0,0.10*F_MAX2]]

	MaxStep_MuscleVelocity = 1 # percentage of positive maximum.
	MuscleVelocity_Bounds =[[-5*lo1,5*lo1],[-1*lo2,1*lo2]]

"""
c_{1} &= -\frac{3g}{2L} \\
c_{2} &= \frac{3}{ML^2} \\
c_{3} &= \cos(\rho_1) \\
c_{4} &= \cos(\rho_2)

"""

c1 = -(3*g)/(2*L)
c2 = 3/(M*L**2)
c3 = np.cos(α1)
c4 = np.cos(α2)

'''
R_{1} &= r_{1}\left(\theta,\theta_{PS}=\frac{\pi}{2}\right) \\
R_{2} &= r_{2}\left(\theta,\theta_{PS}=\frac{\pi}{2}\right)  \\
K_{T,1} &= \frac{F_{max,1}c^{T}}{l_{T,o,1}}\left(1 - \exp{\left(\frac{-T_1}{F_{max,1}c^{T}k^{T}}\right)}\right) \hspace{1em} \text{(Variable stiffness coefficient from Zajac (1989) ODE for tendon force)} \\
v_{MTU,1} &= \text{sgn}\left(-r_1(\theta)\right)\cdot\dot{\theta}\cdot\sqrt{\left(\frac{\partial r_1}{\partial\theta}\right)^2 + r_1^2(\theta)} \\
K_{T,2} &= \frac{F_{max,2}c^{T}}{l_{T,o,2}}\left(1 - \exp{\left(\frac{-T_2}{F_{max,2}c^{T}k^{T}}\right)}\right) \hspace{1em} \text{(Variable stiffness coefficient from Zajac (1989) ODE for tendon force)} \\
v_{MTU,2} &= \text{sgn}\left(-r_2(\theta)\right)\cdot\dot{\theta}\cdot\sqrt{\left(\frac{\partial r_2}{\partial\theta}\right)^2 + r_2^2(\theta)} \\

'''

def R1(X):
	return(r1(X[0])) #
def dR1_dx1(X):
	return(dr1_dθ(X[0]))
def R2(X):
	return(r2(X[0])) #
def dR2_dx1(X):
	return(dr2_dθ(X[0]))
def KT_1(X):
	return((F_MAX1*cT/lTo1)*(1-np.exp(-X[2]/(F_MAX1*cT*kT)))) # NOT NORMALIZED (in N/m)
def v_MTU1(X,r1,dr1_dx1):
	return(np.sign(-r1)*X[1]*np.sqrt(dr1_dx1**2 + r1**2)) # NOT NORMALIZED (in m/s)
def KT_2(X):
	return((F_MAX2*cT/lTo2)*(1-np.exp(-X[3]/(F_MAX2*cT*kT)))) # NOT NORMALIZED (in N/m)
def v_MTU2(X,r2,dr2_dx1):
	return(np.sign(-r2)*X[1]*np.sqrt(dr2_dx1**2 + r2**2)) # NOT NORMALIZED (in m/s)

"""

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
def d2X1_dt2(X):
	return(dX2_dt(X))
def dX2_dt(X,U=None):
	if U is None:
		return(c1*np.sin(X[0]) + c2*R1(X)*X[2] + c2*R2(X)*X[3])
	else:
		return(c1*np.sin(X[0]) + c2*R1(X)*U[0] + c2*R2(X)*U[1])
def d2X2_dt2(X):
	return(c1*np.cos(X[0])*dX1_dt(X) + c2*dR1_dx1(X)*dX1_dt(X)*X[2] + c2*R1(X)*dX3_dt(X)\
			+ c2*dR2_dx1(X)*dX1_dt(X)*X[3] + c2*R2(X)*dX4_dt(X))
def dX3_dt(X,U=None):
	if U is None:
		return(KT_1(X)*(v_MTU1(X,R1(X),dR1_dx1(X)) - c3*X[6]))
	else:
		return(KT_1(X)*(v_MTU1(X,R1(X),dR1_dx1(X)) - c3*U[0]))
def dX4_dt(X,U=None):
	if U is None:
		return(KT_2(X)*(v_MTU2(X,R2(X),dR2_dx1(X)) - c4*X[7]))
	else:
		return(KT_2(X)*(v_MTU2(X,R2(X),dR2_dx1(X)) - c4*U[1]))

### Reference Trajectory ###

r = lambda t: Amp*np.sin(Freq*t) + Base
dr = lambda t: Amp*Freq*np.cos(Freq*t)
d2r = lambda t: -Amp*Freq**2*np.sin(Freq*t)
d3r = lambda t: -Amp*Freq**3*np.cos(Freq*t)

############################

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
	"""
	dZ2(t,X,U) = c1*np.sin(X[0]) + c2*R1(X)*U[0] + c2*R2(X)*U[1] - dA1(t,X)
	"""
	return(dX2_dt(X) - dA1(t,X))
def A2(t,X):
	return(Z1(t,X) + dA1(t,X) - c1*np.sin(X[0]) - k2*Z2(t,X))
def dA2(t,X):
	return(dZ1(t,X) + d2A1(t,X) - c1*np.cos(X[0])*dX1_dt(X) - k2*dZ2(t,X))
def Z3(t,X):
	return(c2*R1(X)*X[2] + c2*R2(X)*X[3] - A2(t,X))
def A3(t,X,r1,r2,kt_1,kt_2):
	dr1_dx1 = dR1_dx1(X)
	dr2_dx1 = dR2_dx1(X)
	return(Z2(t,X) - dA2(t,X) + k3*Z3(t,X) \
		+ c2*dr1_dx1*dX1_dt(X)*X[2] + 	c2*dr2_dx1*dX1_dt(X)*X[3] \
			+ c2*r1*kt_1*v_MTU1(X,r1,dr1_dx1) + c2*r2*kt_2*v_MTU2(X,r2,dr2_dx1))

def return_constraint_variables_tension_driven(t,X):
	Coefficient1 = c2*R1(X)
	Coefficient2 = c2*R2(X)
	Constraint = A2(t,X)
	return(Coefficient1,Coefficient2,Constraint)
def return_constraint_variables_muscle_velocity_driven(t,X):
	r1 = R1(X)
	r2 = R2(X)
	kt_1 = KT_1(X)
	kt_2 = KT_2(X)
	Coefficient1 = c2*c3*r1*kt_1
	Coefficient2 = c2*c4*r2*kt_2
	Constraint = A3(t,X,r1,r2,kt_1,kt_2)
	return(Coefficient1,Coefficient2,Constraint)

def return_U_muscle_velocity_driven(t:float,X,U,**kwargs):
	"""
	Takes in time scalar (float) (t), state numpy.ndarray (X) of shape (4,), and previous input numpy.ndarray (U) of shape (2,) and returns the input for this time step.

	Enforcing a hyperbolic domain constraint to allow for realistic lengthening/shortenting relationships.
	Input2 = (lo1*0.001)*(lo2*0.001)/Input1 = lo1*lo2/(10^6*Input1)

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Noise - must be an numpy.ndarray of shape (2,). Default is np.zeros((1,2)).

	2) Seed - must be a scalar value. Default is None.

	3) Bounds - must be a (2,2) list with each row in ascending order. Default is given by MuscleVelocity_Bounds.

	4) MaxStep - must be a scalar (int or float). Default is MaxStep_MuscleVelocity.

	"""
	import random
	import numpy as np

	assert np.shape(X) == (4,) and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (4,) numpy.ndarray"
	assert np.shape(U) == (2,) and str(type(U)) == "<class 'numpy.ndarray'>", "U must be a (2,) numpy.ndarray"

	dt = Time2[1]-Time2[0]

	Noise = kwargs.get("Noise",np.zeros((2,)))
	assert np.shape(Noise) == (2,) and str(type(Noise)) == "<class 'numpy.ndarray'>", "Noise must be a (2,) numpy.ndarray"

	Seed = kwargs.get("Seed",None)
	assert type(Seed) in [float,int] or Seed is None, "Seed must be a float or an int or None."
	np.random.seed(Seed)

	Bounds = kwargs.get("Bounds",MuscleVelocity_Bounds)
	assert type(Bounds) == list and np.shape(Bounds) == (2,2), "Bounds for Muscle Velocity Control must be a (2,2) list."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."

	MaxStep = kwargs.get("MaxStep",MaxStep_MuscleVelocity)
	assert type(MaxStep) in [int,float], "MaxStep for Muscle Velocity Controller should be an int or float."

	Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t,X)

	if Constraint1 != 0:
		assert Coefficient1!=0 or Coefficient2!=0, "Error with Coefficients. Shouldn't be zero with nonzero constraint."
	else:
		assert Coefficient1!=0 or Coefficient2!=0, "Error with Constraint. 0 = 0 implies all inputs valid."

	Roots = np.sort(\
				np.array(\
	     			list(\
	     				set(\
		     				np.roots(\
					     				[1,\
						     				-Constraint1/Coefficient1,\
						     					Coefficient2*lo1*lo2*(10**-6)/Coefficient1]\
																							)))))
	Roots = Roots[np.isreal(Roots)]

	if Coefficient1 == 0:
		if Constraint1/Coefficient2 > 0:
			LowerBound = Bounds[0][0]
			UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
		else:
			LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient2)
			UpperBound = Bounds[0][1]
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2]*1000)
	elif Coefficient2 == 0:
		FeasibleInput1 = np.array([Constraint1/Coefficient1]*1000)
		if Constraint1/Coefficient1 < 0:
			LowerBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
			UpperBound = Bounds[1][1]
		else:
			LowerBound = Bounds[1][0]
			UpperBound = (lo1*(0.001)*lo2*(0.001))/(Constraint1/Coefficient1)
		FeasibleInput2 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
	else:
		assert 0 not in Roots, "Zero should not be a root. (Implies Coefficient2 == 0)"
		if len(Roots) in [0,1]:
			SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
										(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
			LowerBound = max(Bounds[0][0], SortedBounds[0])
			UpperBound = min(Bounds[0][1], SortedBounds[1])
			assert UpperBound >= LowerBound, "Error generating bounds. Not feasible!"
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		elif (Roots<0).all() or (Roots>0).all():
			SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
										(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
			LowerBound = max(Bounds[0][0], SortedBounds[0])
			UpperBound = min(Bounds[0][1], SortedBounds[1])
			ConstraintLength1 = Coefficient1/(2*Coefficient2)*(LowerBound**2-Roots[0]**2) \
									- Constraint1/Coefficient2*(LowerBound-Roots[0])
			ConstraintLength1 = ConstraintLength1*(ConstraintLength1>0)
			ConstraintLength2 = Coefficient1/(2*Coefficient2)*(Roots[1]**2-UpperBound**2) \
									- Constraint1/Coefficient2*(Roots[1]-UpperBound)
			ConstraintLength2 = ConstraintLength2*(ConstraintLength2>0)
			assert ConstraintLength1!=0 or ConstraintLength2!=0, \
								"Error generating bounds. Not feasible!"
			N1 = int(np.round(1000*ConstraintLength1/(ConstraintLength1+ConstraintLength2)))
			N2 = 1000-N1
			FeasibleInput1_1 = (Roots[0]-LowerBound)*np.random.rand(N1) + LowerBound
			FeasibleInput1_2 = (UpperBound-Roots[1])*np.random.rand(N2) + Roots[1]
			FeasibleInput1 = np.concatenate([FeasibleInput1_1,FeasibleInput1_2])
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		else: # not((Roots<0).all()) and not((Roots>0).all()):
			SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
					(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
			LowerBound = max(Bounds[0][0], SortedBounds[0],Roots[0])
			UpperBound = min(Bounds[0][1], SortedBounds[1],Roots[1])
			assert UpperBound >= LowerBound, "Error with Bounds. Infeasible!"
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
		plt.ylim(Bounds[1])
		plt.show()

	"""
	Checking to see which inputs have the appropriate allowable step size. In normalized muscle velocity.
	"""
	euclid_dist = np.array(list(map(lambda u1,u2: np.sqrt(((U[0]-u1)/lo1)**2 + ((U[1]-u2)/lo2)**2),\
							FeasibleInput1,FeasibleInput2)))

	feasible_index = np.where(euclid_dist <= \
									(MaxStep*(t>=100*dt) + \
									 	10.0*MaxStep*(50*dt<=t<100*dt) + \
											50.0*MaxStep*(t<50*dt)))
	next_index = np.random.choice(feasible_index[0])
	u1 = FeasibleInput1[next_index]
	u2 = FeasibleInput2[next_index]
	return(np.array([u1,u2]))

def return_initial_U_tension_driven(t:float,X_o,**kwargs):
	"""
	Takes in time scalar (float) (t), initial state numpy.ndarray (X_o) of shape (2,) and returns an initial input (2,) numpy.ndarray.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Seed - must be a scalar value. Default is None.

	2) Bounds - must be a (2,2) list with each row in ascending order. Default is given by Tension_Bounds.

	"""
	import random
	import numpy as np

	Seed = kwargs.get("Seed",None)
	assert type(Seed) in [float,int] or Seed is None, "Seed must be a float or an int or None."
	np.random.seed(Seed)

	Bounds = kwargs.get("Bounds",Tension_Bounds)
	assert type(Bounds) == list and np.shape(Bounds) == (2,2), "Bounds for Tension Control must be a (2,2) list."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."

	assert np.shape(X_o) == (2,) and str(type(X_o)) == "<class 'numpy.ndarray'>", "X_o must be a (2,) numpy.ndarray"

	Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_tension_driven(t,X_o)
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
	else:
		SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
									(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
		LowerBound = max(Bounds[0][0], SortedBounds[0])
		UpperBound = min(Bounds[0][1], SortedBounds[1])
		assert UpperBound >= LowerBound, "Error generating bounds. Not feasible!"
		FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
		FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
								for el in FeasibleInput1])
	index = np.random.choice(range(1000))
	u1 = FeasibleInput1[index]
	u2 = FeasibleInput2[index]
	return(np.array([u1,u2]))
def return_initial_U_muscle_velocity_driven(t:float,X_o,**kwargs):
	"""
	Takes in time scalar (float) (t), initial state numpy.ndarray (X_o) of shape (4,) and returns an initial input (2,) numpy.ndarray.

	Enforcing a hyperbolic domain constraint to allow for realistic lengthening/shortenting relationships.
	Input2 = (lo1*0.001)*(lo2*0.001)/Input1 = lo1*lo2/(10^6*Input1)

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Seed - must be a scalar value. Default is None.

	2) Bounds - must be a (2,2) list with each row in ascending order. Default is given by MuscleVelocity_Bounds.

	"""
	import random
	import numpy as np

	Seed = kwargs.get("Seed",None)
	assert type(Seed) in [float,int] or Seed is None, "Seed must be a float or an int or None."
	np.random.seed(Seed)

	Bounds = kwargs.get("Bounds",MuscleVelocity_Bounds)
	assert type(Bounds) == list and np.shape(Bounds) == (2,2), "Bounds for Muscle Velocity Control must be a (2,2) list."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."

	assert np.shape(X_o) == (4,) and str(type(X_o)) == "<class 'numpy.ndarray'>", "X_o must be a (2,) numpy.ndarray"

	Coefficient1,Coefficient2,Constraint1 = \
	 									return_constraint_variables_muscle_velocity_driven(t,X_o)
	assert np.shape(Bounds)==(2,2), "Bounds must be (2,2)."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."
	if Constraint1 != 0:
		assert Coefficient1!=0 or Coefficient2!=0, "Error with Coefficients. Shouldn't be zero with nonzero constraint."
	else:
		assert Coefficient1!=0 or Coefficient2!=0, "Error with Constraint. 0 = 0 implies all inputs valid."

	Roots = np.sort(\
				np.array(\
	     			list(\
	     				set(\
		     				np.roots(\
					     				[1,\
						     				-Constraint1/Coefficient1,\
						     					Coefficient2*lo1*lo2*(10**-6)/Coefficient1]\
																							)))))
	Roots = Roots[np.isreal(Roots)]

	if Coefficient1 == 0:
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
	elif Coefficient2 == 0:
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
	else:
		assert 0 not in Roots, "Zero should not be a root. (Implies Coefficient2 == 0)"
		if len(Roots) in [0,1]:
			SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
										(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
			LowerBound = max(Bounds[0][0], SortedBounds[0])
			UpperBound = min(Bounds[0][1], SortedBounds[1])
			assert UpperBound >= LowerBound, "Error generating bounds. Not feasible!"
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		elif (Roots<0).all() or (Roots>0).all():
			SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
										(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
			LowerBound = max(Bounds[0][0], SortedBounds[0])
			UpperBound = min(Bounds[0][1], SortedBounds[1])
			ConstraintLength1 = Coefficient1/(2*Coefficient2)*(LowerBound**2-Roots[0]**2) \
									- Constraint1/Coefficient2*(LowerBound-Roots[0])
			ConstraintLength1 = ConstraintLength1*(ConstraintLength1>0)
			ConstraintLength2 = Coefficient1/(2*Coefficient2)*(Roots[1]**2-UpperBound**2) \
									- Constraint1/Coefficient2*(Roots[1]-UpperBound)
			ConstraintLength2 = ConstraintLength2*(ConstraintLength2>0)
			assert ConstraintLength1!=0 or ConstraintLength2!=0, \
								"Error generating bounds. Not feasible!"
			N1 = int(np.round(1000*ConstraintLength1/(ConstraintLength1+ConstraintLength2)))
			N2 = 1000-N1
			FeasibleInput1_1 = (Roots[0]-LowerBound)*np.random.rand(N1) + LowerBound
			FeasibleInput1_2 = (UpperBound-Roots[1])*np.random.rand(N2) + Roots[1]
			FeasibleInput1 = np.concatenate([FeasibleInput1_1,FeasibleInput1_2])
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])
		else: # not((Roots<0).all()) and not((Roots>0).all()):
			SortedBounds = np.sort([(Constraint1-Coefficient2*Bounds[1][0])/Coefficient1,\
										(Constraint1-Coefficient2*Bounds[1][1])/Coefficient1])
			LowerBound = max(Bounds[0][0], SortedBounds[0],Roots[0])
			UpperBound = min(Bounds[0][1], SortedBounds[1],Roots[1])
			assert UpperBound >= LowerBound, "Error with Bounds. Infeasible!"
			FeasibleInput1 = (UpperBound-LowerBound)*np.random.rand(1000) + LowerBound
			FeasibleInput2 = np.array([Constraint1/Coefficient2 - (Coefficient1/Coefficient2)*el \
									for el in FeasibleInput1])

	index = np.random.choice(range(1000))
	u1 = FeasibleInput1[index]
	u2 = FeasibleInput2[index]
	return(np.array([u1,u2]))

def plot_MA_values(t,X,**kwargs):
	"""
	Take the numpy.ndarray time array (t) of size (N,) and the state space numpy.ndarray (X) of size (2,N), (4,N), or (8,N), and plots the moment are values of the two muscles versus time and along the moment arm function.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) InputString - must be a string. Used to alter the figure Title. Default is None.
	"""
	import matplotlib.pyplot as plt
	import numpy as np

	assert (np.shape(X)[0] in [2,4,8]) \
				and (np.shape(X)[1] == len(t)) \
					and (str(type(X)) == "<class 'numpy.ndarray'>"), \
			"X must be a (2,N), (4,N), or (8,N) numpy.ndarray, where N is the length of t."

	assert np.shape(t) == (len(t),) and str(type(t)) == "<class 'numpy.ndarray'>", "t must be a (N,) numpy.ndarray."

	InputString = kwargs.get("InputString",None)
	assert InputString is None or type(InputString)==str, "InputString must either be a string or None."
	if InputString is None:
		DescriptiveTitle = "Moment arm equations"
	else:
		assert type(InputString)==str, "InputString must be a string"
		DescriptiveTitle = "Moment arm equations\n(" + InputString + " Driven)"

	fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8,6))
	plt.subplots_adjust(left = 0.15,hspace=0.1,bottom=0.1)
	plt.suptitle(DescriptiveTitle)

	ax1.plot(np.linspace(0,np.pi*(160/180),1001),\
				np.array(list(map(lambda x1: R1([x1]),np.linspace(0,np.pi*(160/180),1001)))),\
				'0.70')
	ax1.plot(np.linspace(min(X[0,:]),max(X[0,:]),101),\
				np.array(list(map(lambda x1: R1([x1]),np.linspace(min(X[0,:]),max(X[0,:]),101)))),\
				'g',lw=3)
	ax1.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
	ax1.set_xticklabels([""]*len(ax1.get_xticks()))
	ax1.set_ylabel("Moment Arm for\n Muscle 1 (m)")

	"""
	Note: Need to Transpose X in order for Map to work.
	"""

	ax2.plot(t,np.array(list(map(lambda X: R1(X),X.T))),'g')
	ax2.set_ylim(ax1.get_ylim())
	ax2.set_yticks(ax1.get_yticks())
	ax2.set_yticklabels([""]*len(ax1.get_yticks()))
	ax2.set_xticklabels([""]*len(ax2.get_xticks()))

	ax3.plot(np.linspace(0,np.pi*(160/180),1001),\
				np.array(list(map(lambda x1: R2([x1]),np.linspace(0,np.pi*(160/180),1001)))),\
				'0.70')
	ax3.plot(np.linspace(min(X[0,:]),max(X[0,:]),101),\
				np.array(list(map(lambda x1: R2([x1]),np.linspace(min(X[0,:]),max(X[0,:]),101)))),\
				'r',lw=3)
	ax3.set_xticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi])
	ax3.set_xticklabels([r"$0$",r"$\frac{\pi}{4}$",r"$\frac{\pi}{2}$",r"$\frac{3\pi}{4}$",r"$\pi$"])
	ax3.set_xlabel("Joint Angle (rads)")
	ax3.set_ylabel("Moment Arm for\n Muscle 2 (m)")

	ax4.plot(t,np.array(list(map(lambda X: R2(X),X.T))),'r')
	ax4.set_ylim(ax3.get_ylim())
	ax4.set_yticks(ax3.get_yticks())
	ax4.set_yticklabels([""]*len(ax3.get_yticks()))
	ax4.set_xlabel("Time (s)")
	return(fig,[ax1,ax2,ax3,ax4])

def animate_muscle_velocity_driven(t,X,U,**kwargs):
	"""
	Takes in Time (t - numpy.ndarray of shape (N,)), the state array (X - numpy.ndarray of shape (4,N)), and the input array (U - numpy.ndarray of shape (2,N)) and animates constraint equation over time.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Bounds - must be a (2,2) list with each row in ascending order. Default is given by MuscleVelocity_Bounds.

	"""

	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib.animation as animation
	import matplotlib.patches as patches
	import time

	assert np.shape(X) == (4,len(t)) and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (4,N) numpy.ndarray"

	Bounds = kwargs.get("Bounds",MuscleVelocity_Bounds)
	assert type(Bounds) == list and np.shape(Bounds) == (2,2), "Bounds for Muscle Velocity Control must be a (2,2) list."
	assert Bounds[0][0]<Bounds[0][1],"Each set of bounds must be in ascending order."
	assert Bounds[1][0]<Bounds[1][1],"Each set of bounds must be in ascending order."

	dt = t[1]-t[0]
	fig = plt.figure(figsize=(10,8))
	ax1 = plt.gca()

	DescriptiveTitle = "Plotting Constraints vs. Time\nMuscle Velocity Driven"

	ax1.set_title(DescriptiveTitle,Fontsize=20,y=1)

	#Hyperbolic Constraint/Bounding Constraints
	Input1 = list(np.linspace(Bounds[0][0],Bounds[0][1],1000001))
	Input1.remove(0)
	Input1 = np.array(Input1)
	ax1.plot(Input1,lo1*lo2*0.001**2/Input1,'r',lw=2)
	Input2 = list(np.linspace(Bounds[1][0],Bounds[1][1],1000001))
	Input2.remove(0)
	Input2 = np.array(Input2)
	ax1.plot(lo1*lo2*0.001**2/Input2,Input2,'r',lw=2)

	ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][0],Bounds[1][0]],'k--')
	ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][1],Bounds[1][1]],'k--')
	ax1.plot([Bounds[0][0],Bounds[0][0]],[Bounds[1][0],Bounds[1][1]],'k--')
	ax1.plot([Bounds[0][1],Bounds[0][1]],[Bounds[1][0],Bounds[1][1]],'k--')

	Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[0],X[:,0])
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
		assert UpperBound>=LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
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
		assert UpperBound>=LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
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
	feasible = plt.Circle((U[:,0]),radius=MaxStep_MuscleVelocity,Color='b',alpha=0.5)
	ax1.add_patch(feasible)
	cline, = plt.plot(FeasibleInput1,FeasibleInput2,'b',lw=2)
	TimeText = plt.text(0.1,0.1,"t = " + str(t[0]),fontsize=16)
	chosenpoint, = plt.plot(U[:,0],c='k',marker='o')
	ax1.set_xlabel(r'$v_{m,1}$',fontsize=14)
	ax1.set_ylabel(r'$v_{m,2}$',fontsize=14)
	ax1.set_xlim([Bounds[0][0]-0.10*(np.diff(Bounds[0])[0]/2),\
					Bounds[0][1]+0.10*(np.diff(Bounds[0])[0]/2)])
	ax1.set_ylim([Bounds[1][0]-0.10*(np.diff(Bounds[1])[0]/2),\
					Bounds[1][1]+0.10*(np.diff(Bounds[1])[0]/2)])
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_aspect('equal')

	def animate(i):
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[i],X[:,i])
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
			assert UpperBound>=LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
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
			assert UpperBound>=LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
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
		feasible.center = (U[:,i])
		if i<10:
			feasible.radius = 10*MaxStep_MuscleVelocity
		else:
			feasible.radius = MaxStep_MuscleVelocity
		cline.set_xdata(FeasibleInput1)
		cline.set_ydata(FeasibleInput2)
		chosenpoint.set_xdata(U[0,i])
		chosenpoint.set_ydata(U[1,i])
		TimeText.set_text("t = " + str(t[i]))
		return feasible,cline,chosenpoint,TimeText,


	# Init only required for blitting to give a clean slate.
	def init():
		ax1.plot(Input1,lo1*lo2*0.001**2/Input1,'r',lw=2)
		ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][0],Bounds[1][0]],'k--')
		ax1.plot([Bounds[0][0],Bounds[0][1]],[Bounds[1][1],Bounds[1][1]],'k--')
		ax1.plot([Bounds[0][0],Bounds[0][0]],[Bounds[1][0],Bounds[1][1]],'k--')
		ax1.plot([Bounds[0][1],Bounds[0][1]],[Bounds[1][0],Bounds[1][1]],'k--')
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[0],X[:,0])
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
			assert UpperBound>=LowerBound, "Error creating bounds - UpperBound should always be greater than LowerBound."
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
		feasible = plt.Circle((U[:,0]),radius=MaxStep_MuscleVelocity,Color='b',alpha=0.5)
		feasible.set_visible(False)
		cline, = plt.plot(FeasibleInput1,FeasibleInput2,'b',lw=2)
		cline.set_visible(False)
		chosenpoint, = plt.plot(U[:,0],c='k',marker='o')
		chosenpoint.set_visible(False)
		TimeText = plt.text(0.75,0.75,"t = " + str(t[0]),fontsize=16)
		TimeText.set_visible(False)
		return feasible,cline,chosenpoint,TimeText,

	ani = animation.FuncAnimation(fig, animate, np.arange(1, np.shape(X)[1],1), init_func=init,interval=1, blit=False)
	plt.show()
def plot_individual_constraint_versus_time_muscle_velocity_driven(t,X,**kwargs):
	"""
	A⋅u₁ + B⋅u₂ = C

	Takes in Time (t - numpy.ndarray of shape (N,)), the state array (X - numpy.ndarray of shape (4,N)) and plots the constraint equation and its coefficients over time.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines whether or not a function handle is returned. Default is False.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert np.shape(X)[0] == 4 and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (4,N) numpy.ndarray"

	Return = kwargs.get("Return",False)
	assert type(Return) == bool, "Return must be either True or False."

	DescriptiveTitle = "Plotting Coefficients/Constraints vs. Time\nMuscle Velocity Driven\n$A\cdot u_{1} + B\cdot u_{2} = C$"
	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
	plt.subplots_adjust(wspace=0.4,top=0.8)
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.95)

	A,B,C = [],[],[]
	for i in range(np.shape(X)[1]):
		Coefficient1,Coefficient2,Constraint1 = return_constraint_variables_muscle_velocity_driven(t[i],X[:,i])
		A.append(Coefficient1)
		B.append(Coefficient2)
		C.append(Constraint1)

	ax1.plot(t[:np.shape(X)[1]],A,'r',lw=2)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$1^{st}$ Coefficient")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:np.shape(X)[1]],B,'b',lw=2)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$2^{nd}$ Coefficient")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	ax3.plot(t[:np.shape(X)[1]],C,'k',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel("Constraint")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_individual_coefficient2_versus_time_muscle_velocity_driven(t,X,**kwargs):
	"""
	B = c2⋅c4⋅R2(X)⋅KT_2(X)

	Takes in Time (t - numpy.ndarray of shape (N,)), the state array (X - numpy.ndarray of shape (4,N)) and plots the 2nd Coefficient of the Constraint Equation over time as well as its components.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines whether or not a function handle is returned. Default is False.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert np.shape(X)[0] == 4 and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (4,N) numpy.ndarray"

	Return = kwargs.get("Return",False)
	assert type(Return) == bool, "Return must be either True or False."

	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
	plt.subplots_adjust(top=0.8,hspace=0.4,bottom=0.1,left=0.1,right=0.975,wspace=0.4)
	DescriptiveTitle = "Plotting $2^{nd}$ Coefficient vs. Time\n$B=c_{2}c_{4}R_{2}(\\vec{x}(t))K_{T,2}(\\vec{x}(t))$"
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.95)

	r2,kt_2,B = [],[],[]
	for i in range(np.shape(X)[1]):
		_,Coefficient2,_ = return_constraint_variables_muscle_velocity_driven(t[i],X[:,i])
		B.append(Coefficient2)
		r2.append(R2(X[:,i]))
		kt_2.append(KT_2(X[:,i]))

	ax1.plot(t[:np.shape(X)[1]],r2,'b',lw=2)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$R_{2}(\vec{x}(t))$")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:np.shape(X)[1]],kt_2,'b',lw=2)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$K_{T,2}(\vec{x}(t))$")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	ax3.plot(t[:np.shape(X)[1]],B,'b',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel(r"$2^{nd}$ Coefficient")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_individual_coefficient1_versus_time_muscle_velocity_driven(t,X,**kwargs):
	"""
	A = c2⋅c3⋅R1(X)⋅KT_1(X)

	Takes in Time (t - numpy.ndarray of shape (N,)), the state array (X - numpy.ndarray of shape (4,N)) and plots the 1st Coefficient of the Constraint Equation over time as well as its components.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines whether or not a function handle is returned. Default is False.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert np.shape(X)[0] == 4 and str(type(X)) == "<class 'numpy.ndarray'>", "X must be a (4,N) numpy.ndarray"

	Return = kwargs.get("Return",False)
	assert type(Return) == bool, "Return must be either True or False."

	fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
	plt.subplots_adjust(top=0.8,hspace=0.4,bottom=0.1,left=0.075,right=0.975,wspace=0.4)
	DescriptiveTitle = "Plotting $1^{st}$ Coefficient vs. Time\n$A=c_{2}c_{3}R_{1}(\\vec{x}(t))K_{T,1}(\\vec{x}(t))$"
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.95)

	r1,kt_1,B = [],[],[]
	for i in range(np.shape(X)[1]):
		Coefficient1,_,_ = return_constraint_variables_muscle_velocity_driven(t[i],X[:,i])
		B.append(Coefficient1)
		r1.append(R1(X[:,i]))
		kt_1.append(KT_1(X[:,i]))

	ax1.plot(t[:np.shape(X)[1]],r1,'r',lw=2)
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$R_{1}(\vec{x}(t))$")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t[:np.shape(X)[1]],kt_1,'r',lw=2)
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$K_{T,1}(\vec{x}(t))$")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	ax3.plot(t[:np.shape(X)[1]],B,'r',lw=2)
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.set_ylabel(r"$1^{st}$ Coefficient")
	ax3.set_xticks(ax1.get_xticks())
	ax3.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()

def plot_states(t,X,**kwargs):
	"""
	Takes in a numpy.ndarray for time (t) of shape (N,) and the numpy.ndarray for the state space (X) of shape (M,N), where M is the number of states and N is the same length as time t. Returns a plot of the states.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines if the function returns a function handle. Default is False.

	2) InputString - must be a string. Input to the DescriptiveTitle that can be used to personalize the title. Default is None.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert (np.shape(X)[0] in [2,4,8]) \
				and (np.shape(X)[1] == len(t)) \
					and (str(type(X)) == "<class 'numpy.ndarray'>"), \
			"X must be a (2,N), (4,N), or (8,N) numpy.ndarray, where N is the length of t."


	Return = kwargs.get("Return",False)
	assert type(Return)==bool, "Return must be either True or False."

	InputString = kwargs.get("InputString",None)
	assert InputString is None or type(InputString)==str, "InputString must either be None or a str."

	NumStates = np.shape(X)[0]
	NumRows = int(np.ceil(NumStates/5))
	if NumStates < 5:
		NumColumns = NumStates
	else:
		NumColumns = 5

	ColumnNumber = [el%5 for el in np.arange(0,NumStates,1)]
	RowNumber = [int(el/5) for el in np.arange(0,NumStates,1)]
	Units = ["(Rads)","(Rads/s)","(N)","(N)","(m)","(m)","(m/s)","(m/s)"]
	if InputString is None:
		DescriptiveTitle = "Plotting States vs. Time"
	else:
		assert type(InputString)==str, "InputString must be a string"
		DescriptiveTitle = InputString + " Driven"
	fig, axes = plt.subplots(NumRows,NumColumns,figsize=(3*NumColumns,2*NumRows + 2))
	plt.subplots_adjust(top=0.85,bottom=0.15,left=0.075,right=0.975)
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.975)
	if NumStates <=5:
		for j in range(NumStates):
			axes[ColumnNumber[j]].spines['right'].set_visible(False)
			axes[ColumnNumber[j]].spines['top'].set_visible(False)
			axes[ColumnNumber[j]].plot(t,X[j,:])
			if ColumnNumber[j]!=0:
				axes[ColumnNumber[j]].set_xticklabels(\
									[""]*len(axes[ColumnNumber[j]].get_xticks()))
			else:
				axes[ColumnNumber[j]].set_xlabel("Time (s)")
			axes[ColumnNumber[j]].set_title(r"$x_{" + str(j+1) + "}$ " + Units[j])

	else:
		for j in range(NumStates):
			axes[RowNumber[j],ColumnNumber[j]].spines['right'].set_visible(False)
			axes[RowNumber[j],ColumnNumber[j]].spines['top'].set_visible(False)
			axes[RowNumber[j],ColumnNumber[j]].plot(t,X[j,:])
			if not(RowNumber[j] == RowNumber[-1] and ColumnNumber[j]==0):
				axes[RowNumber[j],ColumnNumber[j]].set_xticklabels(\
									[""]*len(axes[RowNumber[j],ColumnNumber[j]].get_xticks()))
			else:
				axes[RowNumber[j],ColumnNumber[j]].set_xlabel("Time (s)")
			axes[RowNumber[j],ColumnNumber[j]].set_title(r"$x_{" + str(j+1) + "}$ "+ Units[j])
		if NumStates%5!=0:
			[fig.delaxes(axes[RowNumber[-1],el]) for el in range(ColumnNumber[-1]+1,5)]

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_inputs(t,U,**kwargs):
	"""
	Takes in a numpy.ndarray for time (t) of shape (N,) and the numpy.ndarray for the input (U) (NOT NECESSARILY THE SAME LENGTH AS t). Returns a plot of the states.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines if the function returns a function handle. Default is False.

	2) InputString - must be a string. Input to the DescriptiveTitle that can be used to personalize the title. Default is None.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert (np.shape(U)[0] == 2) \
				and (np.shape(U)[1] == len(t)) \
					and (str(type(U)) == "<class 'numpy.ndarray'>"), \
			"X must be a (2,N) numpy.ndarray, where N is the length of t."

	Return = kwargs.get("Return",False)
	assert type(Return)==bool, "Return must be either True or False."

	InputString = kwargs.get("InputString",None)
	assert InputString is None or type(InputString)==str, "InputString must either be None or a str."

	if InputString is None:
		DescriptiveTitle = "Plotting Inputs vs. Time"
	else:
		assert type(InputString)==str, "InputString must be a string"
		DescriptiveTitle = InputString + " vs. Time"
	fig, (ax1,ax2) = plt.subplots(1,2,figsize=(13,5))
	plt.subplots_adjust(top=0.9,hspace=0.4,bottom=0.1,left=0.075,right=0.975)
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.975)

	ax1.plot(t,U[0,:],'r',lw=2)
	ax1.plot([-1,t[-1]+1],[0,0],'k--',lw=0.5)
	ax1.set_xlim([t[0],t[-1]])
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.set_ylabel(r"$u_1$")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t,U[1,:],'g',lw=2)
	ax2.plot([-1,t[-1]+1],[0,0],'k--',lw=0.5)
	ax2.set_xlim([t[0],t[-1]])
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.set_ylabel(r"$u_2$")
	ax2.set_xticks(ax1.get_xticks())
	ax2.set_xticklabels([""]*len(ax1.get_xticks()))

	if Return == True:
		return(fig)
	else:
		plt.show()
def plot_l_m_comparison(t,X,**kwargs):

	"""
	Takes in a numpy.ndarray for time (t) of shape (N,) and the numpy.ndarray for the input (U) (NOT NECESSARILY THE SAME LENGTH AS t). Returns a plot of the states.

	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	1) Return - must be a bool. Determines if the function returns a function handle. Default is False.

	2) InputString - must be a string. Input to the DescriptiveTitle that can be used to personalize the title. Default is None.

	"""
	import numpy as np
	import matplotlib.pyplot as plt

	assert (np.shape(X)[0] >= 2) \
				and (np.shape(X)[1] == len(t)) \
					and (str(type(X)) == "<class 'numpy.ndarray'>"), \
			"X must be a (M,N) numpy.ndarray, where M is greater than or equal to 2 and N is the length of t."

	Return = kwargs.get("Return",False)
	assert type(Return)==bool, "Return must be either True or False."

	InputString = kwargs.get("InputString",None)
	assert InputString is None or type(InputString)==str, "InputString must either be None or a str."
	if InputString is None:
		DescriptiveTitle = "Muscle vs. Musculotendon Lengths"
	else:
		DescriptiveTitle = "Muscle vs. Musculotendon Lengths\n" + InputString + " Driven"

	L_m = kwargs.get("MuscleLengths",None)
	assert (L_m is None) or (str(type(L_m))=="<class 'numpy.ndarray'>" and np.shape(L_m)==(2,len(t))), "L_m must either be a numpy.ndarray of size (2,N) or left as None (Default)."

	V_m = kwargs.get("MuscleVelocities",None)
	assert (V_m is None) or (str(type(V_m))=="<class 'numpy.ndarray'>" and np.shape(V_m)==(2,len(t))), "V_m must either be a numpy.ndarray of size (2,N) or left as None (Default)."

	assert L_m is not None or V_m is not None, "Error! Need to input some length/velocity measurement for the muscles."

	fig,[[ax1,ax2],[ax3,ax4]] = plt.subplots(2,2,figsize = (14,7))
	plt.suptitle(DescriptiveTitle,Fontsize=20,y=0.975)

	if L_m is None:
		"""
		This is for the muscle velocity driven controller. These values of initial muscle length are estimates taken to be the optimal muscle lengths. We will need to run some sensitivity analysis to ensure that this does not drastically effect the deviations from the MTU estimate.
		"""
		l_m1 = integrate.cumtrapz(V_m[0,:],t,initial = 0) + np.ones(len(t))*lo1
		l_m2 = integrate.cumtrapz(V_m[1,:],t,initial = 0) + np.ones(len(t))*lo2
	else:
		l_m1 = L_m[0,:]
		l_m2 = L_m[1,:]

	"""
	Note: X must be transposed in order to run through map()
	"""
	ax1.plot(t,l_m1,'r',t,integrate.cumtrapz(np.array(list(map(lambda X: v_MTU1(X,R1(X),dR1_dx1(X)),X.T))),\
						t,initial=0) + np.ones(len(t))*l_m1[0], 'b')
	ax1.set_ylabel(r"$l_{m,1}/l_{MTU,1}$ (m)")
	ax1.set_xlabel("Time (s)")

	ax2.plot(t,l_m1-integrate.cumtrapz(np.array(list(map(lambda X: v_MTU1(X,R1(X),dR1_dx1(X)),X.T))),\
						t,initial=0) - np.ones(len(t))*l_m1[0], 'k')
	ax2.set_ylabel("Error (m)")
	ax2.set_xlabel("Time (s)")

	ax3.plot(t,l_m2,'r',t,integrate.cumtrapz(np.array(list(map(lambda X: v_MTU2(X,R2(X),dR2_dx1(X)),X.T))),\
						t,initial=0) + np.ones(len(t))*l_m2[0], 'b')
	ax3.set_ylabel(r"$l_{m,2}/l_{MTU,2}$ (m)")
	ax3.set_xlabel("Time (s)")

	ax4.plot(t,l_m2-integrate.cumtrapz(np.array(list(map(lambda X: v_MTU2(X,R2(X),dR2_dx1(X)),X.T))),\
						t,initial=0) - np.ones(len(t))*l_m2[0], 'k')
	ax4.set_ylabel("Error (m)")
	ax4.set_xlabel("Time (s)")

	if Return == True:
		return(fig)
	else:
		plt.show()

def save_figures(BaseFileName,figs):
	import os.path
	from matplotlib.backends.backend_pdf import PdfPages
	i = 1
	FileName = BaseFileName + "_" + "{:0>2d}".format(i) + ".pdf"
	if os.path.exists(FileName) == True:
		while os.path.exists(FileName) == True:
			i += 1
			FileName = BaseFileName + "_" + "{:0>2d}".format(i) + ".pdf"
	PDFFile = PdfPages(FileName)
	if len(figs)==1:
		PDFFile.savefig(figs[0])
	else:
		[PDFFile.savefig(fig) for fig in figs]
	PDFFile.close()

def return_length_of_nonzero_array(X):
	"""
	Takes in a numpy.ndarray X of shape (m,n) and returns the length of the array that removes any trailing zeros.
	"""
	assert str(type(X))=="<class 'numpy.ndarray'>", "X should be a numpy array"
	assert np.shape(X)[1]!=1, "X should be a wide rectangular array. (m,1) is a column, therefore a nonzero X of this shape will return 1 (trivial solution). Transpose X to properly identify nonzero array length."
	assert np.shape(X)!=(1,1), "Check input. Should not be of shape (1,1) (trivial solution)."
	if (X!=np.zeros(np.shape(X))).all():
		return(np.shape(X)[1])
	else:
		return(np.argmax((X == np.zeros(np.shape(X))).sum(axis=0) == np.shape(X)[0]))

AnotherIteration2 = True
AttemptNumber2 = 0
while AnotherIteration2 == True:
	N2 = N_seconds*100 + 1
	Time2 = np.linspace(0,N_seconds,N2)
	dt = Time2[1]-Time2[0]
	"""
	l_m1[0] = lo1 and l_m2[0] = lo2. This is a floating parameter that will need sensitivity analysis!
	"""
	X2 = np.zeros((4,len(Time2)))
	U2 = np.zeros((2,len(Time2)))
	temp_Tension = return_initial_U_tension_driven(0,np.array([Base, Amp*Freq]))
	X2[:,0] = [Base, Amp*Freq, temp_Tension[0], temp_Tension[1]]
	U2[:,0] = return_initial_U_muscle_velocity_driven(Time2[0],X2[:,0])

	AddNoise2 = False
	if AddNoise2 == True:
	    np.random.seed(seed=None)
	    NoiseArray2 = np.random.normal(loc=0.0,scale=0.2,size=(2,len(Time2)))
	else:
	    NoiseArray2 = np.zeros((2,len(Time2)))

	try:
		StartTime = time.time()
		for i in range(len(Time2)-1):
			U2[:,i+1] = return_U_muscle_velocity_driven(Time2[i],X2[:,i],U2[:,i],Noise = NoiseArray2[:,i])
			X2[:,i+1] = X2[:,i] + dt*np.array([dX1_dt(X2[:,i]),\
												dX2_dt(X2[:,i]),\
												dX3_dt(X2[:,i],U=U2[:,i+1]),\
												dX4_dt(X2[:,i],U=U2[:,i+1])])
			statusbar(i,len(Time2)-1,StartTime=StartTime,Title="Vm Controlled")
		AnotherIteration2 = False
	except:
		AttemptNumber2 += 1
		print("\n")
		print("Attempt #" + str(int(AttemptNumber2)) + " Failed.\n")
		if AttemptNumber2 == 20: AnotherIteration2 = False
ArrayLength2 = return_length_of_nonzero_array(X2)
X2 = X2[:,:ArrayLength2]
U2 = U2[:,:ArrayLength2]
Time2 = Time2[:ArrayLength2]
print('\n')

if ArrayLength2>50:
	plt.figure(figsize = (9,7))
	plt.title("Underdetermined Forced-Pendulum Example",\
	                fontsize=16,color='gray')
	plt.plot(Time2,X2[0,:],'g',lw=2)
	plt.plot(np.linspace(0,Time2[-1],1001),\
			r(np.linspace(0,Time2[-1],1001)),\
				'r--')
	plt.xlabel("Time (s)")
	plt.ylabel("Desired Measure")
	plt.legend([r"Output $y = x_{1}$ (mm Velocity)",r"Reference $r(t) = \frac{\pi}{24}\sin(2\pi t) + \frac{\pi}{2}$"],loc='best')

	plt.figure(figsize = (9,7))
	plt.title('Error vs. Time')
	plt.plot(Time2, r(Time2)-X2[0,:],color='g')
	plt.xlabel("Time (s)")
	plt.ylabel("Error")
	plt.legend(["Muscle Velocity Driven"],loc='best')

	fig1_2,[ax1_2,ax2_2,ax3_2,ax4_2] = plot_MA_values(Time2,X2,InputString = "Normalized Muscle Velocities")
	fig2_2 = plot_states(Time2,X2,Return=True,InputString = "Normalized Muscle Velocities")
	fig3_2 = plot_inputs(Time2,U2*np.concatenate([np.ones((1,len(Time2)))/lo1,np.ones((1,len(Time2)))/lo2],axis=0),\
							Return=True,InputString = "Normalized Muscle Velocities")
	fig4_2 = plot_l_m_comparison(Time2,X2,MuscleVelocities = U2, Return=True, InputString = "Muscle Velocity")

BaseFileName = "ReferenceTracking_MuscleVelocityDrivenPendulumExample"
figs=[manager.canvas.figure
         for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
if len(figs)>=1:
	save_figures(BaseFileName,figs)
plt.close('all')
# plt.show()
