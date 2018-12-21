import map

def group(vectorsignals, Treshold = 0):
	"""Function that returns the grouped vectors of signals,
	grouping is done module by module"""
	GroupedVectorSignals = [0.]*(111*111)
	for indexfiber in range(len(vectorsignals)):
		indexmodule,indexfiber2 = map.mapmodulefiber(indexfiber)
		if vectorsignals[indexfiber] >= Treshold:
			GroupedVectorSignals[indexmodule] += vectorsignals[indexfiber]
	return GroupedVectorSignals

def mapgroupedXY(indexmodule):
	"""Function that given the index module return X andd Y of the module"""
	ymodule = (1.2*(indexmodule//111)) + 0.6
	xmodule = (1.2*(indexmodule-((indexmodule//111)*111))) + 0.6
	return xmodule,ymodule

def distancegrouped(index1,index2):
	"""Function that given two indeces returns the distance between
	corresponding modules (both scintillating and Cherenkvo)"""
	X1,Y1 = mapgroupedXY(index1)
	X2,Y2 = mapgroupedXY(index2)
	distance = ((X1-X2)**2+(Y1-Y2)**2)**0.5
	return distance

def isinradiusgrouped(radius,fixedindex,newindex):
	"""Function that finds if a fiber is in a given circle respect
	to a fixed fiber"""
	if distancegrouped(fixedindex,newindex)<radius:
		return True
	else:
		return False

def energyatdistancegrouped(radius,fixedindex,vectorsignals):
	"""Function that computes the energy deposited in a circle 
	delimited by two fibers, vectorsignals must be a python list"""
	energy = 0.0
	for indexfibers in range(len(vectorsignals)):
		if isinradiusgrouped(radius,fixedindex,indexfibers):
			energy += vectorsignals[indexfibers] 
	return energy
