import map
import mapgroup
import numpy as np
import machinelearning

def buildcluster(fixedindex,vectorsignals,treshold = 0.001):
	"""Function that given a fixed fiber builds the 
	cluster around it, a cluster is given by a fixed 
	central fiber (or module), a radius and a signal"""
	signalcluster = 0
	maxdistance = 0
	radius = 0.6
	gain = 0
	faildattempt = 0
	while True: #find radius of cluster
		energyatdistance = mapgroup.energyatdistancegrouped(radius,fixedindex,vectorsignals)
		gain = energyatdistance - signalcluster
		if gain >= treshold: 
			signalcluster += gain
			radius += 0.6
			#print str(radius) + " " + str(signalcluster) + " " + str(gain) #if you want to check step by step uncomment
			faildattempt = 0
		else:
			faildattempt = faildattempt + 1
			print "One failed attempt to build cluster"
			radius += 0.6
			if faildattempt == 2:
				maxdistance = radius - 1.2
				print "----------Building cluster complete----------"
				print(fixedindex, maxdistance, signalcluster)
				break
			else:
				continue
	cluster = Cluster(fixedindex, maxdistance, signalcluster)
	return cluster

class Cluster:
	"""Cluster class, a cluster is made of a fixed central index,
	a radius around it, a scintillation signal, a Cherenkov signal,
	a list of fibers (or modules), an ID (electron, hadron),
	an energy calibrated for Scin and Cher and a total reconstructed
	energy"""
	def __init__(self, fixedindex, radius, scinsignal):
		self.fixedindex = fixedindex
		self.radius = radius
		self.scinsignal = scinsignal
		self.chersignal = None
		self.modules = []
		self.ID = None
		self.EnergyScin = None
		self.EnergyCher = None
		self.Energy = None
		self.MLEnergy = None
		self.MLScinEnergy = None
		self.MLCherEnergy = None
		self.secondMLEnergy = None
		self.secondMLScinEnergy = None
		self.secondMLCherEnergy = None

	def computeEnergy(self, Chi, sigmaScin = 0, sigmaCher = 0, weighted = False):
		"""Function to compute the energy of the cluster, to properly perform it
		the cluster ID is needed"""
		if self.ID == "electron":
			if weighted == True:
				self.Energy = (self.EnergyScin*sigmaScin+self.EnergyCher*sigmaCher)/(sigmaCher+sigmaScin)
			else:
				self.Energy = (self.EnergyScin+self.EnergyCher)/2
		elif self.ID == "hadron":
			self.Energy = (self.EnergyScin-Chi*self.EnergyCher)/(1-Chi)

	def compute_ID(self, namefile):
		self.ID = machinelearning.find_ID(self, namefile)

	def compute_MLEnergy(self, namefile):
		self.MLEnergy, self.MLScinEnergy, self.MLCherEnergy, self.secondMLEnergy, self.secondMLScinEnergy, self.secondMLCherEnergy = machinelearning.find_energy(self, namefile)

	def computeEnergyScin(self, scincalibconstant):
		self.EnergyScin = self.scinsignal*scincalibconstant

	def computeEnergyCher(self, chercalibconstant):
		self.EnergyCher = self.chersignal*chercalibconstant

	def compute_chersignal(self, chervectorsignal):
		chersignal = 0
		for indexmodule in range(len(chervectorsignal)):
			if indexmodule in self.modules:
				chersignal += chervectorsignal[indexmodule]
		self.chersignal = chersignal

	def find_modules(self, vectorsignal):
		"""Function to find modules belonging to a cluster,
		result is a list of modules where only the indeces of
		modules in the cluster appear"""
		for indexmodule in range(len(vectorsignal)):
			if mapgroup.isinradiusgrouped(self.radius,self.fixedindex,indexmodule):
				self.modules.append(indexmodule)
			else:
				self.modules.append(0.0) 

	def printposition(self):
		"""Function to print cluster position"""
		X,Y = mapgroup.mapgroupedXY(self.fixedindex)
		print "Cluster at x=" + str(X) + " cm y= " + str(Y) + " cm " + str(self.radius) + " cm large\n"

'''
Event class was thought to work with many clusters, not implemented yet
class Event:
	"""Event class, an event is an Event number, a list of clusters,
	a remaining vector signal scintillating and Cherenkov and
	a total energy reconstructed"""
	def __init__(self, eventnumber):
		self.eventnumber = eventnumber
		self.eventclusters = []
		self.remainingVectorSignal = []
		self.remainingVectorSignalCher = []
		self.Energy = None

	def addcluster(self, cluster):
		"""Function to add one cluster in the event's eventclusters"""
		self.eventclusters.append(cluster)

	def createevent(self, vectorsignals, treshold = 20):
		"""Function to add all the clusters in event's eventclusters"""
		if max(vectorsignals) > treshold:
			startindex = vectorsignals.index(max(vectorsignals)) 
			startcluster = buildcluster(startindex, vectorsignals)
			startcluster.find_modules(vectorsignals)
			self.addcluster(startcluster)
			vectorsignals = list(np.subtract(vectorsignals,startcluster.modules))
		for moduleindex in range(len(vectorsignals)):
			if vectorsignals[moduleindex] > treshold:
				cluster = buildcluster(moduleindex, vectorsignals)
				cluster.find_modules(vectorsignals)
				self.addcluster(cluster)
				vectorsignals = list(np.subtract(vectorsignals,cluster.modules))
		print "----------Event " + str(self.eventnumber) +" complete: " + str(len(self.eventclusters)) + " clusters found"

	def computeremainingvectorsignal(self, vectorsignals):
		"""Function to compute remaining scin vector signal after cluster subtraction"""
		self.remainingVectorSignal = vectorsignals
		for cluster in self.eventclusters:
			self.remainingVectorSignal = list(np.subtract(self.remainingVectorSignal,cluster.modules))

	def computeremainingvectorsignalCher(self, vectorsignalscher):
		"""Function to compute remaining Cher vector signal after cluster subtraction"""
		self.remainingVectorSignalCher = vectorsignalscher
		for cluster in self.eventclusters:
			self.remainingVectorSignalCher = list(np.subtract(self.remainingVectorSignalCher,cluster.modules))

	def compute_EventEnergy(self, vectorsignalcher, Chi, scincalibconstant, chercalibconstant, namefile):
		"""Function to compute Event energy"""
		Energy = 0
		for cluster in self.eventclusters:
			cluster.compute_chersignal(vectorsignalcher)
			cluster.computestep1(chervectorsignal, scincalibconstant, chercalibconstant)
			cluster.computeEnergyScin(scincalibconstant)
			cluster.computeEnergyCher(chercalibconstant)
			cluster.compute_ID(namefile)
			cluster.computeEnergy(Chi)
			Energy = Energy + cluster.computeEnergy
		#self.Energy = Energy + sum(remainingVectorSignal)*scincalibconstant + sum(remainingVectorSignalCher)*chercalibconstant
		self.Energy = Energy
'''
