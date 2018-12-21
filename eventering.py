import map
import numpy as np
import machinelearning

def buildevent(vectorsignal, vectorsignalcher):
	"""Function to create an event"""

	scinsignal = sum(vectorsignal)
	chersignal = sum(vectorsignalcher)

	event = Event(scinsignal, chersignal) 

	return event

class Event:
	"""Event class, an event is made of a scintillation signal, a Cherenkov signal,
	an ID (electron, hadron), an energy calibrated for Scin and Cher and a total reconstructed
	energy"""
	def __init__(self, scinsignal, chersignal):
		self.scinsignal = scinsignal
		self.chersignal = chersignal
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
		self.TMLEnergy = None
		self.TMLScinEnergy = None
		self.TMLCherEnergy = None
		self.TsecondMLEnergy = None
		self.TsecondMLScinEnergy = None
		self.TsecondMLCherEnergy = None
		self.minEnergy = None

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

	def computeminEnergy(self, a, b):
		"""Function to compute energy from parameters found with minimization"""
		self.minEnergy = a*self.scinsignal+b*self.chersignal

	def compute_ID(self, namefile, noftrained):
		self.ID = machinelearning.find_ID(self, namefile, noftrained)

	def compute_MLEnergy(self, namefile, noftrained):
		self.MLEnergy, self.MLScinEnergy, self.MLCherEnergy, self.secondMLEnergy, self.secondMLScinEnergy, self.secondMLCherEnergy, MLFEM = machinelearning.find_energy(self, namefile, noftrained)
		return MLFEM
    #not using it -> going to cancel it
	#def compute_TMLEnergy(self, Tnamefile):
	#	self.TMLEnergy, self.TMLScinEnergy, self.TMLCherEnergy, self.TsecondMLEnergy, self.TsecondMLScinEnergy, self.TsecondMLCherEnergy = machinelearning.find_energy(self, Tnamefile)

	def computeEnergyScin(self, scincalibconstant):
		self.EnergyScin = self.scinsignal*scincalibconstant

	def computeEnergyCher(self, chercalibconstant):
		self.EnergyCher = self.chersignal*chercalibconstant

class Trainedevent:
	"""Trained Event class: event + truth FEM"""
	def __init__(self, scinsignal, chersignal, FEM):
		self.FEM = FEM
		self.scinsignal = scinsignal
		self.chersignal = chersignal
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
		self.TMLEnergy = None
		self.TMLScinEnergy = None
		self.TMLCherEnergy = None
		self.TsecondMLEnergy = None
		self.TsecondMLScinEnergy = None
		self.TsecondMLCherEnergy = None
		self.minEnergy = None

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

	def computeminEnergy(self, a, b):
		"""Function to compute energy from parameters found with minimization"""
		self.minEnergy = a*self.scinsignal+b*self.chersignal

	def compute_ID(self, namefile, noftrained):
		self.ID = machinelearning.find_ID(self, namefile, noftrained)

	def compute_MLEnergy(self, namefile, noftrained):
		self.MLEnergy, self.MLScinEnergy, self.MLCherEnergy, self.secondMLEnergy, self.secondMLScinEnergy, self.secondMLCherEnergy = machinelearning.find_energy(self, namefile, noftrained)
    #not using it -> going to cancel it
	#def compute_TMLEnergy(self, Tnamefile):
	#	self.TMLEnergy, self.TMLScinEnergy, self.TMLCherEnergy, self.TsecondMLEnergy, self.TsecondMLScinEnergy, self.TsecondMLCherEnergy = machinelearning.find_energy(self, Tnamefile)

	def computeEnergyScin(self, scincalibconstant):
		self.EnergyScin = self.scinsignal*scincalibconstant

	def computeEnergyCher(self, chercalibconstant):
		self.EnergyCher = self.chersignal*chercalibconstant

def buildTrainedevent(vectorsignal, vectorsignalcher, FEM):
	"""Function to create a TRAINED event: enent + truth fem"""

	scinsignal = sum(vectorsignal)
	chersignal = sum(vectorsignalcher)
                   
	trainedevent = Trainedevent(scinsignal, chersignal, FEM) 

	return trainedevent
