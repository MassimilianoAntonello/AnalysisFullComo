def choosecalibration(material, answermincalib):
	"""Function to return calibration constants"""
	if material == "G4_Fe":
		scincalibconstant = 22.459 #MeV/MeV
		chercalibconstant = 21.654 #MeV/Cpe
		Chi = 0.40
		HEcher = 0.41
		HEscin = 0.76
		hadconteinment = 0.9367
		emconteinment = 0.9982
		if answermincalib == "20":
			a = 28.44
			b = -3.73
		if answermincalib == "40":
			a = 29.76
			b = -5.4
		if answermincalib == "60":
			a = 30.54
			b = -6.41
		if answermincalib == "80":
			a = 31.04
			b = -7.06
	if material == "G4_Pb":
		scincalibconstant = 28.0705 #MeV/MeV
		chercalibconstant = 28.135 #MeV/Cpe
		Chi = 0.2
		HEcher = 0.3093
		HEscin = 0.8635
		hadconteinment = 0.964
		emconteinment = 0.997
		if answermincalib == "20":
			a = 33.15
			b = -4.07
		if answermincalib == "40":
			a = 33.83
			b = -4.9
		if answermincalib == "60":
			a = 34.17
			b = -5.38
		if answermincalib == "80":
			a = 34.53
			b = -5.8
	if material == "Brass": 
		scincalibconstant = 23.832 #MeV/MeV
		chercalibconstant = 22.838 #MeV/Cpe
		Chi = 0.37
		HEcher = 0.37
		HEscin = 0.77
		hadconteinment = 0.9403
		emconteinment = 0.99817
		if answermincalib == "20":
			a = 30.01
			b = -4.04
		if answermincalib == "40":
			a = 31.22
			b = -5.65
		if answermincalib == "60":
			a = 31.96
			b = -6.56
		if answermincalib == "80":
			a = 32.48
			b = -7.2
	if material == "G4_W": 
		scincalibconstant = 46.19 #MeV/MeV
		chercalibconstant = 46.65 #MeV/Cpe
		Chi = 0.18
		HEcher = 0.22
		HEscin = 0.86
		hadconteinment = 0.969
		emconteinment = 0.997
		if answermincalib == "20":
			a =  53.05
			b =  -4.97		
		if answermincalib == "40":
			a = 54.26 #at 40 GeV pi
			b = -6.79 #at 40 GeV pi
		if answermincalib == "60":
			a = 55.02
			b = -7.94
		if answermincalib == "80":
			a = 55.21
			b = -8.21
	if material == "G4_Pt": 
		scincalibconstant = 51.104 #MeV/MeV
		chercalibconstant = 50.561 #MeV/Cpe
		Chi = 0.12
		HEcher = 0.238
		HEscin = 0.910
		hadconteinment = 0.97139
		emconteinment = 0.996529
		if answermincalib == "20":
			a =  55.49
			b =  -3.22		
		if answermincalib == "40":
			a =  56.64 #at 40 GeV pi
			b = -4.87 #at 40 GeV pi
		if answermincalib == "60":
			a = 56.96
			b = -5.40
		if answermincalib == "80":
			a = 57.15
			b = -5.61
	if material == "G4_Cu": 
		scincalibconstant = 24.584 #MeV/MeV
		chercalibconstant = 23.526 #MeV/Cpe
		Chi = 0.37
		HEcher = 0.37
		HEscin = 0.767
		hadconteinment = 0.944
		emconteinment = 0.998
		if answermincalib == "20":
			a =  0
			b =  0		
		if answermincalib == "40":
			a = 0 #at 40 GeV pi
			b = 0 #at 40 GeV pi
		if answermincalib == "60":
			a = 0 
			b = 0
		if answermincalib == "80":
			a = 0 
			b = 0
	print "Using "+ str(material) + " calibration constants."
	return scincalibconstant, chercalibconstant, Chi, HEcher, HEscin, hadconteinment, emconteinment, a, b