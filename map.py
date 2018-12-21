import math

def mapmodulefiber(index):
	"""Function that given the VectorSignals returns the
	fiber number (from 0 to 63) and the module number (from 0 to 71x71)"""
	indexmodule = (index)//64
	indexfiber = (index)-indexmodule*64
	return indexmodule,indexfiber;

def mapXY(index):
	"""Function that given the index module and index fiber
	returns X Y positions of the fiber (both Scintillating and Cherenkov)"""
	indexmodule,indexfiber = mapmodulefiber(index)
	yfiber = 0.075+(0.15*(indexfiber//8))
	xfiber = 0.075+(0.15*(indexfiber-((indexfiber//8)*8)))
	ymodule = (1.2*(indexmodule//111))
	xmodule = (1.2*(indexmodule-((indexmodule//111)*111)))
	X = xmodule+xfiber
	Y = ymodule+yfiber
	return X,Y

def distance(index1,index2):
	"""Function that given two indeces returns the distance between 
	corresponding fibers (both scintillating and Cherenkov)"""
	X1,Y1 = mapXY(index1)
	X2,Y2 = mapXY(index2)
	distance = ((X1-X2)**2+(Y1-Y2)**2)**0.5
	return distance

def areaatdistance(distance,radiusstep):
	"""Function that computes the area of the ring delimited by two
	fibers"""
	previouseradius = (distance//radiusstep)*radiusstep 
	area1 = math.pi*(previouseradius**2)
	area2 = math.pi*(distance**2)
	areaatdistance = area2-area1
	return areaatdistance 

def isinradius(radius,fixedindex,newindex,radiusstep = 1.0):
	"""Function that finds if a fiber is in a given ring respect
	to a fixed fiber"""
	previouseradius = (radius//radiusstep)*radiusstep
	nextradius = ((radius//radiusstep)*radiusstep)+radiusstep
	if distance(fixedindex,newindex)>previouseradius and distance(fixedindex,newindex)<nextradius:
		return True
	else:
		return False

def energyatdistance(distance,radiusstep,fixedindex,vectorsignals):
	"""Function that computes the energy deposited in a ring 
	delimited by two fibers, vectorsignals must be a python list"""
	energy = 0.0
	for indexfibers in range(len(vectorsignals)):
		if isinradius(distance,fixedindex,indexfibers):
			energy += vectorsignals[indexfibers] 
	return energy

def energydensityatdistance(distance,radiusstep,vectorsignals):
	"""Function that computes the energy deposited in a ring
	delimited by two fibers divided by the area of the ring,
	vectorsignals must be a python list"""
	return energydensityatdistance(distance,radiusstep,fixedindex,vectorsignals)/areaatdistance(distance,radiusstep)
