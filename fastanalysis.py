"""
Fast Analysis of a Dual-Readout (Como Version) - Lorenzo Pezzotti
Version 1.0

Usage:
  fastanalysis.py -h | --help
  fastanalysis.py --version
  fastanalysis.py [-d <WhatToDo>]
  fastanalysis.py [-m <Material>]

Options:
  -h --help            Display this help and exit
  --version            Display version and exit
  -d <WhatToDo>        calibration, training, analysis, minimization, resolution
"""

import numpy as np 
import sys
from scipy.optimize import minimize
import random
import os
import shutil
import glob
import docopt
import datetime
from array import array
from ROOT import gROOT, TFile, TTree, TF1, TFitter, TF2, gPad, gStyle, gSystem
import map
import mapgroup
import clustering
import machinelearning
import ROOTHistograms
import service
import eventering
import questions
import calibrationconstants
import calibrationconstantsQGSP

argv = docopt.docopt(__doc__, version="fastanalysis.py 1.0")
ToDo = str(argv['-d'])
print "Doing "+str(ToDo)+". \n"
physicslist = raw_input("Insert physics list (FTFP, QGSP): ")
machine = raw_input("On what machine running? (mac, linux, office) ")
if machine == "mac":
	path = str("/Users/lorenzo/cernbox/work/Git-to-Mac/AnalysisFullCalorimeter")
if machine == "linux":
	path = str("/home/lorenzo/cernbox/work/Git-to-Mac/AnalysisFullCalorimeter")
if machine == "office":
	path = str("/media/geant4-mc-infn/DataStorage/lorenzo/cernbox/work/Git-to-Mac/AnalysisFullCalorimeter")
if str(physicslist) == "QGSP":
	if machine == "mac":
		path = str("/Users/lorenzo/cernbox/work/Git-to-Mac/AnalysisFullCalorimeter/QGSP")
	if machine == "linux":
		path = str("/home/lorenzo/cernbox/work/Git-to-Mac/AnalysisFullCalorimeter/QGSP")
	if machine == "office":
		path = str("/media/geant4-mc-infn/DataStorage/lorenzo/cernbox/work/Git-to-Mac/AnalysisFullCalorimeter/QGSP")

endtime = datetime.datetime.now()

#---------------------------------------------------------------------------------------------------
def calibration(answer):
	"""Function to perform calibration with single answer file"""
	#Electromagnetic Energy containments

	#Take electron file
	gROOT.Reset()
	inputFile = TFile(str(answer.electronnamefile)) #root input file
	tree = TTree() #Tree name B4
	inputFile.GetObject("B4", tree)

	EnergiesEscaped = []
	for Event in range(int(answer.electronNofEvents)):

		tree.GetEntry(Event)
		EnergyEscaped = tree.EscapedEnergy #Total energy escaped
		EnergiesEscaped.append(EnergyEscaped)
		if Event % 500 == 0:
			print "Estimating EM containment: Event " + str(Event) + " of " + str(NofEventsProcessed) + " \n"
		if int(Event) == 0:
			PrimaryEnergy = tree.PrimaryParticleEnergy #MC truth primary particle energy

	MeanEnergyLeak = np.mean(EnergiesEscaped) #energy leak average
	EMEnergyContainment = (PrimaryEnergy-MeanEnergyLeak)/PrimaryEnergy
	print "EM energy containment = " + str(EMEnergyContainment) + " percent \n"

	#Set parameters
	fastchercalibconstant = 0 #Done without building clusters
	fastscincalibconstant = 0 #Done without building clusters
	firstcounter = 0

	#Estimate Cherenkov and Scintillation calibration constants
	for Event in range(int(answer.electronNofEvents)):

		tree.GetEntry(Event)

		#Set values of the tree
		PrimaryParticleName = tree.PrimaryParticleName # MC truth: primary particle Geant4 name
		EnergyEscaped = tree.EscapedEnergy #Total energy escaped
		PrimaryParticleEnergy = (tree.PrimaryParticleEnergy) - EnergyEscaped # MC truth: primary particle energy, energy containment FOR ELECTRONS
		EnergyTot = tree.EnergyTot # Total energy deposited in calorimeter
		Energyem = tree.Energyem # Energy deposited by the em component
		EnergyScin = tree.EnergyScin # Energy deposited in Scin fibers (not Birk corrected)
		EnergyCher = tree.EnergyCher # Energy deposited in Cher fibers (not Birk corrected)
		NofCherenkovDetected = tree.NofCherenkovDetected # Total Cher p.e. detected
		VectorSignals = tree.VectorSignals # Vector of energy deposited in Scin fibers (Birk corrected)
		VectorSignalsCher = tree.VectorSignalsCher # Vector of Cher p.e. detected in Cher fibers

		#For output file
		FirstPrimaryParticleName = PrimaryParticleName
		FirstPrimaryParticleEnergy = tree.PrimaryParticleEnergy

		if Event % 500 == 0:
			print "--------------------------------------------------------\n"
			print "Processing event " + str(Event) + " of " + str(NofEventsProcessed) + ": " + PrimaryParticleName + " energy " + str(PrimaryParticleEnergy) + " MeV" +"\n"

		#Estimating calibration constants using all the calorimeter signal
		if sum(VectorSignalsCher) == 0.0:
			continue
			print "One event with no Cherenkov signal"

		fastchercalibconstant += PrimaryParticleEnergy/sum(VectorSignalsCher)
		fastscincalibconstant += PrimaryParticleEnergy/sum(VectorSignals)

		firstcounter = firstcounter+1

	#Finalize calibration constants
	fastchercalibconstant = fastchercalibconstant/firstcounter
	fastscincalibconstant = fastscincalibconstant/firstcounter

	print "End of calibration constants -> now estimating h/e values."

	#Estimating h/e Cher and scin values and Chi factor
	#Set root file and tree
	gROOT.Reset()
	inputFile = TFile(answer.pionnamefile) #root input file
	tree = TTree()
	inputFile.GetObject("B4", tree)

	#Hadronic Energy containments
	EnergiesEscaped[:] = []
	for Event in range(int(answer.pionNofEvents)):

		tree.GetEntry(Event)
		EnergyEscaped = tree.EscapedEnergy #Total energy escaped
		EnergiesEscaped.append(EnergyEscaped)
		if Event % 500 == 0:
			print "Estimating HAD containment: Event " + str(Event) + " of " + str(NofEventsProcessed) + " \n"
		if Event == 0:
			PrimaryEnergy = tree.PrimaryParticleEnergy #MC truth primary particle energy

	MeanEnergyLeak = np.mean(EnergiesEscaped) #energy leak average
	HADEnergyContainment = (PrimaryEnergy-MeanEnergyLeak)/PrimaryEnergy
	print "Hadronic energy containment = " + str(HADEnergyContainment) + " percent \n"


	#For parameters estimation
	SecondPrimaryParticleName = PrimaryParticleName
	SecondPrimaryParticleEnergy = tree.PrimaryParticleEnergy

	#Set parameters and counters
	NofEvents = tree.GetEntries()
	counter1 = 0
	counter2 = 0
	counter3 = 0
	counter4 = 0
	counter5 = 0
	counter6 = 0
	counter7 = 0
	counter8 = 0
	HistFastScinenergy1 = []
	HistFastScinenergy2 = []
	HistFastScinenergy3 = []
	HistFastScinenergy4 = []
	HistFastCherenergy1 = []
	HistFastCherenergy2 = []
	HistFastCherenergy3 = []
	HistFastCherenergy4 = []
	HistFastScinenergy5 = []
	HistFastScinenergy6 = []
	HistFastScinenergy7 = []
	HistFastScinenergy8 = []
	HistFastCherenergy5 = []
	HistFastCherenergy6 = []
	HistFastCherenergy7 = []
	HistFastCherenergy8 = []
	rmsScin1 = 0
	rmsCher1 = 0
	rmsScin2 = 0
	rmsCher2 = 0
	rmsScin3 = 0
	rmsCher3 = 0
	rmsScin4 = 0
	rmsCher4 = 0
	rmsScin5 = 0
	rmsCher5 = 0
	rmsScin6 = 0
	rmsCher6 = 0
	rmsScin7 = 0
	rmsCher7 = 0
	rmsScin8 = 0
	rmsCher8 = 0
	entriesScin1 = 0
	entriesCher1 = 0
	entriesScin2 = 0
	entriesCher2 = 0
	entriesScin3 = 0
	entriesCher3 = 0
	entriesScin4 = 0
	entriesCher4 = 0
	entriesScin5 = 0
	entriesCher5 = 0
	entriesScin6 = 0
	entriesCher6 = 0
	entriesScin7 = 0
	entriesCher7 = 0
	entriesScin8 = 0
	entriesCher8 = 0
	fem = []
	fem1 = []
	fem2 = []
	fem3 = []
	fem4 = []
	fem5 = []
	fem6 = []
	fem7 = []
	fem8 = []
	Fem1 = 0
	Fem2 = 0
	Fem3 = 0
	Fem4 = 0
	Fem5 = 0
	Fem6 = 0
	Fem7 = 0
	Fem8 = 0
	rmsFem1 = 0
	rmsFem2 = 0
	rmsFem3 = 0
	rmsFem4 = 0
	entriesFem1 = 0
	entriesFem2 = 0
	entriesFem3 = 0
	entriesFem4 = 0
	entriesFem5 = 0
	entriesFem6 = 0
	entriesFem7 = 0
	entriesFem8 = 0
	Histhescin = []
	Histhecher = []
	HistChi = []
	hescin = 0
	hecher = 0
	hescinerror = 0
	hechererror = 0
	chi = 0
	chierror = 0
	entrieshescin = 0
	entrieshecher = 0
	entrieschi = 0

	for Event in range(int(answer.pionNofEvents)):

		tree.GetEntry(Event)

		#Set values of tree
		PrimaryParticleName = tree.PrimaryParticleName # MC truth: primary particle Geant4 name
		EnergyEscaped = tree.EscapedEnergy #Escaped energy
		PrimaryParticleEnergy = tree.PrimaryParticleEnergy - EnergyEscaped # MC truth: primary particle energy, 99% energy containment
		EnergyTot = tree.EnergyTot # Total energy deposited in calorimeter
		Energyem = tree.Energyem2 # Energy deposited by the em component
		EnergyScin = tree.EnergyScin # Energy deposited in Scin fibers (not Birk corrected)
		EnergyCher = tree.EnergyCher # Energy deposited in Cher fibers (not Birk corrected)
		NofCherenkovDetected = tree.NofCherenkovDetected # Total Cher p.e. detected
		VectorSignals = tree.VectorSignals # Vector of energy deposited in Scin fibers (Birk corrected)
		VectorSignalsCher = tree.VectorSignalsCher # Vector of Cher p.e. detected in Cher fibers

		if Event % 500 == 0:
			print "--------------------------------------------------------\n"
			print "Processing event " + str(Event) + " of " + str(NofEventsProcessed) + ": " + PrimaryParticleName + " energy " + str(PrimaryParticleEnergy) + " MeV" +"\n"

		#hes hec and chi histogram part
		if PrimaryParticleEnergy == 0.0:
			continue
		Histhescin.append((sum(VectorSignals)*fastscincalibconstant/PrimaryParticleEnergy - Energyem/PrimaryParticleEnergy)/(1-(Energyem/PrimaryParticleEnergy)))
		Histhecher.append(((sum(VectorSignalsCher)*fastchercalibconstant)/PrimaryParticleEnergy - Energyem/PrimaryParticleEnergy)/(1-(Energyem/PrimaryParticleEnergy)))
		#don't nedd to know fem to estimate Chi
		HistChi.append((sum(VectorSignals)*fastscincalibconstant-PrimaryParticleEnergy)/(sum(VectorSignalsCher)*fastchercalibconstant-PrimaryParticleEnergy))

		#Create histograms
		if PrimaryParticleEnergy == 0.0:
			continue
		Fem = Energyem/PrimaryParticleEnergy #Primary particle energy corrected for event containment
		fem.append(Fem)
		if Fem < 0.25 and Fem > 0.225:

			HistFastScinenergy1.append(fastscincalibconstant*sum(VectorSignals))
			HistFastCherenergy1.append(fastchercalibconstant*sum(VectorSignalsCher))

			counter1 += 1

			fem1.append(Energyem/PrimaryParticleEnergy)

		elif Fem < 0.35 and Fem > 0.325:

			HistFastScinenergy2.append(fastscincalibconstant*sum(VectorSignals))
			HistFastCherenergy2.append(fastchercalibconstant*sum(VectorSignalsCher))

			counter2 += 1

			fem2.append(Energyem/PrimaryParticleEnergy)

		elif Fem < 0.45 and Fem > 0.425:

			HistFastScinenergy3.append(fastscincalibconstant*sum(VectorSignals))
			HistFastCherenergy3.append(fastchercalibconstant*sum(VectorSignalsCher))

			counter3 += 1

			fem3.append(Energyem/PrimaryParticleEnergy)

		elif Fem < 0.55 and Fem > 0.525:

			HistFastScinenergy4.append(fastscincalibconstant*sum(VectorSignals))
			HistFastCherenergy4.append(fastchercalibconstant*sum(VectorSignalsCher))

			counter4 += 1

			fem4.append(Energyem/PrimaryParticleEnergy)

		elif Fem < 0.65 and Fem > 0.625:

			HistFastScinenergy5.append(fastscincalibconstant*sum(VectorSignals))
			HistFastCherenergy5.append(fastchercalibconstant*sum(VectorSignalsCher))

			counter5 += 1

			fem5.append(Energyem/PrimaryParticleEnergy)

		elif Fem < 0.75 and Fem > 0.725:

			HistFastScinenergy6.append(fastscincalibconstant*sum(VectorSignals))
			HistFastCherenergy6.append(fastchercalibconstant*sum(VectorSignalsCher))

			counter6 += 1

			fem6.append(Energyem/PrimaryParticleEnergy)

		elif Fem < 0.85 and Fem > 0.825:

			HistFastScinenergy7.append(fastscincalibconstant*sum(VectorSignals))
			HistFastCherenergy7.append(fastchercalibconstant*sum(VectorSignalsCher))

			counter7 += 1

			fem7.append(Energyem/PrimaryParticleEnergy)

	'''
	#Print statistics on fem events
	print "Found " + str(counter1) + " events with 0.2< fem <0.4 \n"
	print "      " + str(counter2) + " events with 0.4< fem <0.6 \n"
	print "      " + str(counter3) + " events with 0.6< fem <0.8 \n"
	print "      " + str(counter4) + " events with 0.8< fem <1.0 \n"
	'''

	#Print histograms 
	ROOTHistograms.create_fastroothistograms(fem, "fem", "fem", "Entries", "Femall.pdf")

	FastScinenergy1, rmsScin1, entriesScin1 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastScinenergy1, "0.225< fem < 0.25", "Energy Scin", "Entries", "0.225-fem-0.25Scin.pdf")
	FastCherenergy1, rmsCher1, entriesCher1 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastCherenergy1, "0.225< fem < 0.25", "Energy Cher", "Entries", "0.225-fem-0.25Cher.pdf")
	
	FastScinenergy2, rmsScin2, entriesScin2 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastScinenergy2, "0.325< fem < 0.35", "Energy Scin", "Entries", "0.325-fem-0.35Scin.pdf")
	FastCherenergy2, rmsCher2, entriesCher2 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastCherenergy2, "0.325< fem < 0.35", "Energy Cher", "Entries", "0.325-fem-0.35Cher.pdf")
	
	FastScinenergy3, rmsScin3, entriesScin3 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastScinenergy3, "0.425< fem < 0.45", "Energy Scin", "Entries", "0.425-fem-0.45Scin.pdf")
	FastCherenergy3, rmsCher3, entriesCher3 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastCherenergy3, "0.425< fem < 0.45", "Energy Cher", "Entries", "0.425-fem-0.45Cher.pdf")
	
	FastScinenergy4, rmsScin4, entriesScin4 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastScinenergy4, "0.525< fem < 0.55", "Energy Scin", "Entries", "0.525-fem-0.55Scin.pdf")
	FastCherenergy4, rmsCher4, entriesCher4 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastCherenergy4, "0.525< fem < 0.55", "Energy Cher", "Entries", "0.525-fem-0.55Cher.pdf")
	
	FastScinenergy5, rmsScin5, entriesScin5 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastScinenergy5, "0.625< fem < 0.65", "Energy Scin", "Entries", "0.6-fem-0.65Scin.pdf")
	FastCherenergy5, rmsCher5, entriesCher5 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastCherenergy5, "0.625< fem < 0.65", "Energy Cher", "Entries", "0.6-fem-0.65Cher.pdf")
	
	FastScinenergy6, rmsScin6, entriesScin6 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastScinenergy6, "0.725< fem < 0.75", "Energy Scin", "Entries", "0.725-fem-0.75Scin.pdf")
	FastCherenergy6, rmsCher6, entriesCher6 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastCherenergy6, "0.725< fem < 0.75", "Energy Cher", "Entries", "0.725-fem-0.75Cher.pdf")
	
	FastScinenergy7, rmsScin7, entriesScin7 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastScinenergy7, "0.825< fem < 0.85", "Energy Scin", "Entries", "0.825-fem-0.85Scin.pdf")
	FastCherenergy7, rmsCher7, entriesCher7 =  ROOTHistograms.create_fastroothistogramsgaus(HistFastCherenergy7, "0.825< fem < 0.85", "Energy Cher", "Entries", "0.825-fem-0.85Cher.pdf")
	
	Fem1, rmsFem1, entriesFem1 = ROOTHistograms.create_fastroothistograms(fem1, "0.225< fem <0.25", "fem", "Entries", "fem1.pdf")
	Fem2, rmsFem2, entriesFem2 = ROOTHistograms.create_fastroothistograms(fem2, "0.325< fem <0.35", "fem", "Entries", "fem2.pdf")
	Fem3, rmsFem3, entriesFem3 = ROOTHistograms.create_fastroothistograms(fem3, "0.425< fem <0.45", "fem", "Entries", "fem3.pdf")
	Fem4, rmsFem4, entriesFem4 = ROOTHistograms.create_fastroothistograms(fem4, "0.525< fem <0.55", "fem", "Entries", "fem4.pdf")
	Fem5, rmsFem5, entriesFem5 = ROOTHistograms.create_fastroothistograms(fem5, "0.625< fem <0.65", "fem", "Entries", "fem5.pdf")
	Fem6, rmsFem6, entriesFem6 = ROOTHistograms.create_fastroothistograms(fem6, "0.725< fem <0.75", "fem", "Entries", "fem6.pdf")
	Fem7, rmsFem7, entriesFem7 = ROOTHistograms.create_fastroothistograms(fem7, "0.825< fem <0.85", "fem", "Entries", "fem7.pdf")

	#Print histograms hesci hecher and chi
	hescin, hescinerror, entrieshescin = ROOTHistograms.create_fastroothistograms(Histhescin, "hescin", "he", "Entries", "hescin.pdf")
	hecher, hechererror, entrieshecher = ROOTHistograms.create_fastroothistograms(Histhecher, "hecher", "he", "Entries", "hecher.pdf")
	chi, chierror, entrieschi = ROOTHistograms.create_fastroothistograms(HistChi, "chi", "chi", "Entries", "chi.pdf")
	hescinerror = hescinerror/(entrieshescin**0.5)
	hechererror = hechererror/(entrieshecher**0.5)
	chierror = chierror/(entrieschi**0.5)
	'''
	#Values grouped in list for ROOT graphs
	FastScinenergy = array('d',[ FastScinenergy2, FastScinenergy3, FastScinenergy4, FastScinenergy5, FastScinenergy6 ])
	FastScinenergyError = array('d',[ rmsScin2/(entriesScin2**0.5), rmsScin3/(entriesScin3**0.5), rmsScin4/(entriesScin4**0.5), rmsScin5/(entriesScin5**0.5), rmsScin6/(entriesScin6**0.5) ])
	FastCherenergy = array('d',[ FastCherenergy2, FastCherenergy3, FastCherenergy4, FastCherenergy5, FastCherenergy6 ])
	FastCherenergyError = array('d',[ rmsCher2/(entriesCher2**0.5), rmsCher3/(entriesCher3**0.5), rmsCher4/(entriesCher4**0.5), rmsCher5/(entriesCher5**0.5), rmsCher6/(entriesCher6**0.5) ])
	fem = array('d',[ Fem2, Fem3, Fem4, Fem5, Fem6 ])
	femError = array('d',[ rmsFem2/(entriesFem2**0.5), rmsFem3/(entriesFem3**0.5), rmsFem4/(entriesFem4**0.5), rmsFem5/(entriesFem5**0.5), rmsFem6/(entriesFem6**0.5) ])
	'''
	#Values grouped in list for ROOT graphs
	FastScinenergy = array('d',[ FastScinenergy2, FastScinenergy3, FastScinenergy4, FastScinenergy5, FastScinenergy6 ])
	FastScinenergyError = array('d',[ rmsScin2/2, rmsScin3/2, rmsScin4/2, rmsScin5/2, rmsScin6/2 ])
	FastCherenergy = array('d',[ FastCherenergy2, FastCherenergy3, FastCherenergy4, FastCherenergy5, FastCherenergy6 ])
	FastCherenergyError = array('d',[ rmsCher2/2, rmsCher3/2, rmsCher4/2, rmsCher5/2, rmsCher6/2 ])
	fem = array('d',[ Fem2, Fem3, Fem4, Fem5, Fem6 ])
	femError = array('d',[ rmsFem2/2, rmsFem3/2, rmsFem4/2, rmsFem5/2, rmsFem6/2 ])

	#Create graphs, compute he Cher and scin values and Chi factor
	Fasthescin, Fasthecher, Fasthescinerror, Fasthechererror, Fasthescin2, Fasthecher2, hesnew, hecnew, hesnewerror, hecnewerror = ROOTHistograms.create_fastgrapherror(FastScinenergy, FastCherenergy, FastScinenergyError, FastCherenergyError, fem, femError, HADEnergyContainment*PrimaryEnergy)
	fastChi = (1 - Fasthescin)/(1 - Fasthecher)
	fastChierror = (Fasthescinerror/Fasthescin+Fasthechererror/Fasthecher)*fastChi
	fastChinew = (1-hesnew)/(1-hecnew)
	fastChinewerror = (hesnewerror/hesnew+hecnewerror/hecnew)*fastChinew

	#Print calibration results on file
	outputfile = open("fastcalibration.txt","w+")
	outputfile.write("Calibration constants for electron + hadron of energies: \n")
	outputfile.write(str(FirstPrimaryParticleEnergy) + " MeV " + "EM containment " + str(EMEnergyContainment) + " \n")
	outputfile.write(str(SecondPrimaryParticleEnergy) + " MeV " + "HAD containment " + str(HADEnergyContainment) + " \n")
	outputfile.write("-------------------------FAST ANALYSIS---------------------------------------------------------\n")
	outputfile.write("fastchercalibconstant = " + str(fastchercalibconstant) + " MeV/Cpe" + " fastscincalibconstant = " + str(fastscincalibconstant) + " MeV/MeV\n")
	outputfile.write("With graph (pol1 fit): fasthescin = " + str(Fasthescin)+"+-"+str(Fasthescinerror) + " fasthecher = " + str(Fasthecher)+"+-"+str(Fasthechererror) + " fastChi = " + str(fastChi)+"+-"+str(fastChierror) + " \n")
	outputfile.write("         (custom fit): fasthescin = " + str(hesnew)+"+-"+str(hesnewerror) + " fasthecher = " + str(hecnew)+"+-"+str(hecnewerror) + " fastchi = " + str(fastChinew)+"+-"+str(fastChinewerror)+ " \n")
	outputfile.write("-----------------------------------------------------------------------------------------------\n")
	outputfile.write("With histogram: fasthescin = " + str(hescin)+"+-"+str(hescinerror) + " fasthecher = " + str(hecher)+"+-"+str(hechererror) + " fastchi = " + str((1-hescin)/(1-hecher))+"+-"+str((hescinerror/(1-hescin)+hechererror/(1-hecher))*(1-hescin)/(1-hecher)) + " \n")
	outputfile.write("-----------------------------------------------------------------------------------------------\n")
	outputfile.write("With direct Chi estimation: "	+ " Chi = " + str(chi)+"+-"+str(chierror) + " \n")
	outputfile.close()

	#Create folder calibration and move output file (eps and txt)
	if not os.path.exists(str(path)+"/"+"fastcalibration"):
		os.makedirs(str(path)+"/"+"fastcalibration")
	if not os.path.exists(str(path)+"/"+"fastcalibration/"+str(answer.absorbermaterial)):
		os.makedirs(str(path)+"/"+"fastcalibration/"+str(answer.absorbermaterial))
	if not os.path.exists(str(path)+"/"+"fastcalibration/"+str(answer.absorbermaterial)+"/"+str(endtime.strftime("%H-%M_%d_%m_%Y"))):
		os.makedirs(str(path)+"/"+"fastcalibration/"+str(answer.absorbermaterial)+"/"+str(endtime.strftime("%H-%M_%d_%m_%Y")))
	foldername = str(path)+"/"+"fastcalibration/"+str(answer.absorbermaterial)+"/"+str(endtime.strftime("%H-%M_%d_%m_%Y"))+"/"+str(answer.FolderName) # foldername asked before
	os.makedirs(str(foldername))
	if not os.path.exists(str(foldername)+"/result"):
		os.makedirs(str(foldername)+"/result")
	for pdf_file in glob.iglob('*.pdf'):
		shutil.move(pdf_file, str(foldername)+"/result")
	for txt_file in glob.iglob('*.txt'):
		shutil.move(txt_file, str(foldername)+"/result")
	'''
	if not os.path.exists(str(foldername)+"/fem"):
		os.makedirs(str(foldername)+"/fem")
	os.chdir(str(foldername)+"/result")
	for fem_file in glob.iglob('*fem*'):
		shutil.move(fem_file, "../fem")
	'''
	inputFile.Close()
	print "End of h/e values estimation. \n"
	return None
#---------------------------------------------------------------------------------------------------
#Section used only if you want to estimate the calibration constants and Chi factor of a module.
#DREAM is calibrated with electrons. The input file must be from electron events with the same energy.
#To estimate h/e Cher and scin values and the Chi factor a second root file with pions (+ or -) is needed.
if ToDo == "calibration":

	answers = [] #list with answers for calibration
	Files = [] #list of all hadron files to study
	Files = raw_input("List all the hadron files to study: ").split()
	NofFiles = int(len(Files))
	electronfile = raw_input("Input electron name file: ")
	counter = 0

	for File in range(NofFiles):

		if counter == 0:
			absorbermaterial = raw_input("Choose absorber material: ")
		
		#Choose ROOT electron file and pick file
		gROOT.Reset()
		electronfilename = str(path)+"/"+"data-rotated/"+str(absorbermaterial)+"/"+str(electronfile)+"/B4.root"
		inputFile = TFile(str(electronfilename)) #root input file
		tree = TTree() #Tree name B4
		inputFile.GetObject("B4", tree)

		#Check if MATERIAL is correct, very important to not mix things!!!
		tree.GetEntry(0)
		geant4material = str(tree.AbsorberMaterial)
		print "Found Geant4 absorber material: "+str(geant4material)+" \n"
		geant4material = geant4material[0:len(absorbermaterial)]
		if geant4material != str(absorbermaterial):
			sys.exit()
		
		#Ask for how many events to be processed
		if counter == 0:
			NofEvents = tree.GetEntries()
			NofEventsProcessed = raw_input(str(NofEvents) + " Events, how many to process: ")

		#To later perform h/e Cher and scin computation and Chi factor	
		namepionfile = str(path)+"/"+"data-rotated/"+str(absorbermaterial)+"/"+str(Files[counter])+"/B4.root"
		gROOT.Reset()
		inputFile = TFile(str(namepionfile)) #root input file
		tree = TTree() #Tree name B4
		inputFile.GetObject("B4", tree)

		#Check if MATERIAL is correct, very important to not mix things!!!
		tree.GetEntry(0)
		geant4material = str(tree.AbsorberMaterial)
		print "Found Geant4 absorber material: "+str(geant4material)+" \n"
		geant4material = geant4material[0:len(absorbermaterial)]
		if geant4material != str(absorbermaterial):
			sys.exit()
		
		#Ask for how many events to be processed
		if counter == 0:
			pionNofEvents = tree.GetEntries()
			pionNofEventsProcessed = raw_input(str(pionNofEvents) + " Events, how many to process: ")

		inputFile.Close()

		#Ask for folder name where to store results
		foldername = str(electronfile)+"-"+str(Files[counter])

		counter = counter + 1

		#Create answer class with parameters for file
		answer = questions.Calibrationanswer(absorbermaterial, electronfilename, NofEventsProcessed, namepionfile, pionNofEventsProcessed, foldername)
		answers.append(answer)
		#End of answers filling
	
	#Now perform calibration on all files
	for answer in answers:
		print "calibration with: "+str(answer.electronnamefile)+" "+str(answer.pionnamefile)
		calibration(answer)	
#---------------------------------------------------------------------------------------------------	
def training(answer):

	if str(physicslist) == "QGSP":
		#charge calibration constants for had containment
		scincalibconstant, chercalibconstant, Chi, HEcher, HEscin, hadconteinment, emconteinment, a, b = calibrationconstantsQGSP.choosecalibration(answer.absorbermaterial, "40")
	else:
		scincalibconstant, chercalibconstant, Chi, HEcher, HEscin, hadconteinment, emconteinment, a, b = calibrationconstants.choosecalibration(answer.absorbermaterial, "40")		

	#list of trained events
	listoftrainedevents = [] 

	#Choose ROOT file and pick file
	gROOT.Reset()
	filename = str(path)+"/"+"data-rotated/"+str(answer.absorbermaterial)+"/"+str(answer.namefile)+"/B4.root"
	inputFile = TFile(str(filename)) #root input file
	tree = TTree() #Tree name B4
	inputFile.GetObject("B4", tree)

	for Event in range(int(answer.NofEvents)):
		
		tree.GetEntry(Event)

		#Set values of tree
		PrimaryParticleName = tree.PrimaryParticleName # MC truth: primary particle Geant4 name
		EnergyEscaped = tree.EscapedEnergy #Escaped energy
		PrimaryParticleEnergy = tree.PrimaryParticleEnergy - EnergyEscaped # MC truth: primary particle energy - escaped energy
		EnergyTot = tree.EnergyTot # Total energy deposited in calorimeter
		Energyem = tree.Energyem # Energy deposited by the em component
		EnergyScin = tree.EnergyScin # Energy deposited in Scin fibers (not Birk corrected)
		EnergyCher = tree.EnergyCher # Energy deposited in Cher fibers (not Birk corrected)
		NofCherenkovDetected = tree.NofCherenkovDetected # Total Cher p.e. detected
		VectorSignals = tree.VectorSignals # Vector of energy deposited in Scin fibers (Birk corrected)
		VectorSignalsCher = tree.VectorSignalsCher # Vector of Cher p.e. detected in Cher fibers

		#if doing training with jets primary particle energy is energy of a single track
		#I have to fix it
		if "jet" in answer.namefile: #if it is a jet file
			b = [s for s in answer.namefile[0:3] if s.isdigit()]
			c = ""
			for s in range(len(b)):
				c = c+b[s]
			d = int(c)
			PrimaryParticleEnergy = float(d*1000 - EnergyEscaped)
		#end of finxing jet energy problem

		if Event % 500 == 0:
			print "--------------------------------------------------------\n"
			print "Processing event " + str(Event) + " of " + str(answer.NofEvents) + ": " + PrimaryParticleName + " energy " + str(PrimaryParticleEnergy) + " MeV" +"\n"

		#Create Trainedevent
		if PrimaryParticleEnergy == 0.0:
			print "One event with wrong primary energy"
			continue
		FEM = Energyem/PrimaryParticleEnergy
		Trainedevent = eventering.buildTrainedevent(VectorSignals, VectorSignalsCher, FEM)
		Trainedevent.computeEnergyScin(scincalibconstant)
		Trainedevent.computeEnergyCher(chercalibconstant)
		if Trainedevent.chersignal == 0.0:
			print "One cluster with no Cherenkov signal"
			continue
		Trainedevent.Energy = PrimaryParticleEnergy # MC truth to create trained events
		if particletype == "electron":
			Trainedevent.ID = "electron"
		else:
			Trainedevent.ID = "hadron"

		#Creat list of trained events, random shuffle it and save in pickle file
		listoftrainedevents.append(Trainedevent)

	random.shuffle(listoftrainedevents)
	machinelearning.picklelistofclusters(listoftrainedevents, answer.trainednamefile)
	print "Pickled " + str(len(listoftrainedevents)) + " trained events in " + answer.trainednamefile +". \n"

	#Create folder trainedevents and move output file into
	if not os.path.exists(str(path)+"/"+"fasttrained-rotated"):
		os.makedirs(str(path)+"/"+"fasttrained-rotated")
	if not os.path.exists(str(path)+"/"+"fasttrained-rotated/"+str(answer.absorbermaterial)):
		os.makedirs(str(path)+"/"+"fasttrained-rotated/"+str(answer.absorbermaterial))
	foldername = str(path)+"/"+"fasttrained-rotated/"+str(answer.absorbermaterial)
	for p_file in glob.iglob('*.p'):
		shutil.move(p_file, str(foldername))
	print "End of training."
	return listoftrainedevents
#---------------------------------------------------------------------------------------------------
#Section used only if you want to ceate a trained list of clusters to later perform machine learning
#ID and energy reconstruction. ROOT events must be single particle events to be classified in as "electron"
#or "hadron" with a given MC truth energy. Enable this part with Docopt.
if ToDo == "training":
	answers = [] #list of answers
	Files = [] #list of all hadron files to study
	Files = raw_input("List all the hadron files to study: ").split()
	NofFiles = int(len(Files))
	Types = [] #list of particle types
	Types = raw_input("List all the particle types: ").split()
	counter = 0

	for File in range(NofFiles):

		if counter == 0:
			absorbermaterial = raw_input("Choose absorber material: ")


		#Set name of pickle file containing trained clusters
		namefile = str(Files[counter])
		particletype = str(Types[counter])

		#Choose ROOT electron file and pick file
		gROOT.Reset()
		file = str(Files[counter])
		filename = str(path)+"/"+"data-rotated/"+str(absorbermaterial)+"/"+str(file)+"/B4.root"
		inputFile = TFile(str(filename)) #root input file
		tree = TTree() #Tree name B4
		inputFile.GetObject("B4", tree)

		#Check if MATERIAL is correct, very important to not mix things!!!
		tree.GetEntry(0)
		geant4material = str(tree.AbsorberMaterial)
		print "Found Geant4 absorber material: "+str(geant4material)+" \n"
		geant4material = geant4material[0:len(absorbermaterial)]
		if geant4material != str(absorbermaterial):
			sys.exit()
		
		#Ask for how many events to be processed
		if counter == 0:
			NofEvents = tree.GetEntries()
			NofEventsProcessed = raw_input(str(NofEvents) + " Events, how many to process: ")

		inputFile.Close()

		if counter == 0:
			picklename = raw_input("Insert pickle name (.p) ")

		counter = counter + 1

		#save answers
		answer = questions.Traininganswer(absorbermaterial, file, namefile, particletype, NofEventsProcessed)
		answers.append(answer)

	#Now perform training on all files
	listofalltrainedevents = [] #list with all trained events over more files
	for answer in answers:
		print "Training with: "+str(answer.namefile)+" \n"
		listofalltrainedevents = listofalltrainedevents + training(answer)
	if len(answers)>1:
		random.shuffle(listofalltrainedevents)
		machinelearning.picklelistofclusters(listofalltrainedevents, str(picklename))
		foldername = str(path)+"/"+"fasttrained-rotated/"+str(answers[0].absorbermaterial)
		for p_file in glob.iglob('*.p'):
			shutil.move(p_file, str(foldername))
		print "Pickled " + str(len(listofalltrainedevents)) + " trained events. \n"
#----------------------------------------------------------------------------------------------------
def minimization(answer):
	
	if str(physicslist) == "QGSP":
		#charge calibration constants for had containment
		scincalibconstant, chercalibconstant, Chi, HEcher, HEscin, hadconteinment, emconteinment, a, b = calibrationconstantsQGSP.choosecalibration(answer.absorbermaterial, "40")
	else:
		scincalibconstant, chercalibconstant, Chi, HEcher, HEscin, hadconteinment, emconteinment, a, b = calibrationconstants.choosecalibration(answer.absorbermaterial, "40")		

	a = 0
	b = 0

	#Parameters
	SumE_real = 0.0
	SumE_real2 = 0.0
	SumC = 0.0
	SumS = 0.0
	SumS2 = 0.0
	SumC2 = 0.0
	SumSC = 0.0
	SumE_realS = 0.0
	SumE_realC = 0.0
	counter = 0

	#Choose ROOT file and pick file
	gROOT.Reset()
	filename = str(path)+"/"+"data-rotated/"+str(answer.absorbermaterial)+"/"+str(answer.namefile)+"/B4.root"
	inputFile = TFile(str(filename)) #root input file
	tree = TTree() #Tree name B4
	inputFile.GetObject("B4", tree)

	#I have to minimize Sum_i(E_real-(aS_i)**2-(bC_i)**2-2abS_iC_i)
	for Event in range(int(answer.NofEvents)):
		tree.GetEntry(Event)

		#Set values of tree -> pass from MeV to GeV scale
		PrimaryParticleName = tree.PrimaryParticleName # MC truth: primary particle Geant4 name
		EnergyEscaped = tree.EscapedEnergy #Energy escaped
		PrimaryParticleEnergy = (tree.PrimaryParticleEnergy-EnergyEscaped)/1000 # MC truth: primary particle energy
		EnergyTot = tree.EnergyTot # Total energy deposited in calorimeter
		Energyem = tree.Energyem # Energy deposited by the em component
		EnergyScin = tree.EnergyScin # Energy deposited in Scin fibers (not Birk corrected)
		EnergyCher = tree.EnergyCher # Energy deposited in Cher fibers (not Birk corrected)
		NofCherenkovDetected = tree.NofCherenkovDetected # Total Cher p.e. detected
		VectorSignals = (tree.VectorSignals) # Vector of energy deposited in Scin fibers (Birk corrected)
		VectorSignalsCher = (tree.VectorSignalsCher) # Vector of Cher p.e. detected in Cher fibers

		if Event % 500 == 0:
			print "--------------------------------------------------------\n"
			print "Processing event " + str(Event) + " of " + str(answer.NofEvents) + ": " + PrimaryParticleName + " energy " + str(PrimaryParticleEnergy) + " GeV" +"\n"

		#This version of the code deals with single particle events
		event = eventering.buildevent(VectorSignals, VectorSignalsCher)

		if event.chersignal == 0.0:
			print "Event with no Cherenkov signal"
			continue
		if PrimaryParticleEnergy == 0.0:
			print "Event with wrong PrimaryParticleEnergy"
			continue
		SumE_real += PrimaryParticleEnergy
		SumE_real2 += PrimaryParticleEnergy**2
		SumS2 += (sum(VectorSignals)/1000)**2
		SumS += sum(VectorSignals)/1000
		SumC += sum(VectorSignalsCher)/1000
		SumC2 += (sum(VectorSignalsCher)/1000)**2
		SumSC += (sum(VectorSignals)/1000)*(sum(VectorSignalsCher)/1000)
		SumE_realS += PrimaryParticleEnergy*(sum(VectorSignals)/1000)
		SumE_realC += PrimaryParticleEnergy*(sum(VectorSignalsCher)/1000)
		counter += 1

	print "Found "+str(counter)+" good events over "+str(answer.NofEvents)+" given events."

	#Print parameters for minimization
	print "Doing minimization with true energy: " + str(PrimaryParticleEnergy)
	print "SumE_real2 = "+str(SumE_real2)
	print "Using: SumS2 = "+str(SumS2)+ " SumC2 = "+str(SumC2)
	print "using: SumSC = "+str(SumSC)+ " SumE_realS = "+str(SumE_realS)+ " SumE_realC = "+str(SumE_realC)

	#Python minimization
	#Create function to minimize
	def functomin(x):
		#E=x[0]S+x[1]C
		#minimizing Sum_i(E_real-E_i)**2/N
		return (SumE_real2+SumS2*x[0]**2+SumC2*x[1]**2+2*SumSC*x[0]*x[1]-2*SumE_realS*x[0]-2*SumE_realC*x[1])/counter
		
	#Starting point
	x0 = [None,None]
	x0[0] = random.uniform(-10,10)
	x0[1] = random.uniform(-10,10)
	print "Starting minimization with SolCG method "+str(x0)+" : "+str(functomin(x0))

	#Minimization
	solCG = minimize(functomin,x0,method='CG')
	print str(solCG.fun)+" minimum found with CG method" #minimum found
	print str(solCG.x)+" parameters found with CG method" #parameters found

	x0[0] = random.uniform(-10,10)
	x0[1] = random.uniform(-10,10)
	print "Starting minimization with solSLSQP method "+str(x0)+" : "+str(functomin(x0))

	solSLSQP = minimize(functomin,x0,method='SLSQP')
	print str(solSLSQP.fun)+" minimum found with SLSQP method" #minimum found
	print str(solSLSQP.x)+ " parameters found with SLSQP method" #parameters found

	func = TF2("",'([0]+[1]*x**2+[2]*y**2+2*[3]*x*y-2*[4]*x-2*[5]*y)/[6]',-10,60,-25,50)
	func.SetParameter(0, SumE_real2)
	func.SetParameter(1, SumS2)
	func.SetParameter(2, SumC2)
	func.SetParameter(3, SumSC)
	func.SetParameter(4, SumE_realS)
	func.SetParameter(5, SumE_realC)
	func.SetParameter(6, counter)

	XAxis = func.GetXaxis()
	XAxis.SetTitle("a")
	YAxis = func.GetYaxis()
	YAxis.SetTitle("b")
	func.Draw("cont4 z")
	gPad.SaveAs("functomin.pdf")

	#Print parameters results on file
	outputfile = open("fastparameters.txt","w+")
	outputfile.write("Calibration paramters with hadron of energies: \n")
	outputfile.write(str(PrimaryParticleEnergy) + " MeV\n")
	outputfile.write("-------------------------FAST ANALYSIS---------------------------------------------------------\n")
	outputfile.write("fasta-CG = " + str(solCG.x[0]) + "GeV/MeV" + " fastb-CG = " + str(solCG.x[1]) + "GeV/Cpe" + " \n")
	outputfile.write("fasta-SLSQP = " + str(solSLSQP.x[0]) + "GeV/MeV" + " fastb-SLSQP = " + str(solSLSQP.x[1]) + "GeV/Cpe" + " \n")
	outputfile.close()

	#Create folder calibration and move output file (eps and txt)
	if not os.path.exists(str(path)+"/"+"fastparameters"):
		os.makedirs(str(path)+"/"+"fastparameters")
	if not os.path.exists(str(path)+"/"+"fastparameters/"+str(answer.absorbermaterial)):
		os.makedirs(str(path)+"/"+"fastparameters/"+str(answer.absorbermaterial))
	if not os.path.exists(str(path)+"/"+"fastcalibration/"+str(answer.absorbermaterial)+"/"+str(endtime.strftime("%H-%M_%d_%m_%Y"))):
		os.makedirs(str(path)+"/"+"fastcalibration/"+str(answer.absorbermaterial)+"/"+str(endtime.strftime("%H-%M_%d_%m_%Y")))
	foldername = str(path)+"/"+"fastparameters/"+str(answer.absorbermaterial)+"/"+str(endtime.strftime("%H-%M_%d_%m_%Y"))+"/"+answer.FolderName # foldername asked before
	os.makedirs(str(foldername))
	for pdf_file in glob.iglob('*.pdf'):
		shutil.move(pdf_file, str(foldername))
	for txt_file in glob.iglob('*.txt'):
		shutil.move(txt_file, str(foldername))
		
	print "End of minimization"
#----------------------------------------------------------------------------------------------------
#Section used for the miniminaztion of parameters to find the best way to estimate hadron energy.
if ToDo == "minimization":
	answers = [] #list of answers
	Files = [] #list of all hadron files to study
	Files = raw_input("List all the hadron files to study: ").split()
	NofFiles = int(len(Files))
	counter = 0

	for File in range(NofFiles):

		if counter == 0:
			absorbermaterial = raw_input("Choose absorber material: ")

		#Choose ROOT electron file and pick file
		gROOT.Reset()
		file = str(Files[counter])
		filename = str(path)+"/"+"data-rotated/"+str(absorbermaterial)+"/"+str(file)+"/B4.root"
		inputFile = TFile(str(filename)) #root input file
		tree = TTree() #Tree name B4
		inputFile.GetObject("B4", tree)

		#Check if MATERIAL is correct, very important to not mix things!!!
		tree.GetEntry(0)
		geant4material = str(tree.AbsorberMaterial)
		print "Found Geant4 absorber material: "+str(geant4material)+" \n"
		geant4material = geant4material[0:len(absorbermaterial)]
		if geant4material != str(absorbermaterial):
			sys.exit()
		
		#Ask for how many events to be processed
		if counter == 0:
			NofEvents = tree.GetEntries()
			NofEventsProcessed = raw_input(str(NofEvents) + " Events, how many to process: ")

		inputFile.Close()

		#Ask for folder name where to store results
		foldername = str(Files[counter])

		counter = counter + 1

		#save answers
		answer = questions.Minimizationanswer(absorbermaterial, file, NofEventsProcessed, foldername)
		answers.append(answer)

	#Now perform training on all files
	for answer in answers:
		print "Minimization with: "+str(answer.namefile)+" \n"
		minimization(answer)		
#----------------------------------------------------------------------------------------------------
def analysis(answer):
	
	if str(physicslist) == "QGSP":
		#charge calibration constants for had containment
		scincalibconstant, chercalibconstant, Chi, HEcher, HEscin, hadconteinment, emconteinment, a, b = calibrationconstantsQGSP.choosecalibration(answer.absorbermaterial, "40")
	else:
		scincalibconstant, chercalibconstant, Chi, HEcher, HEscin, hadconteinment, emconteinment, a, b = calibrationconstants.choosecalibration(answer.absorbermaterial, "40")		

	#Choose ROOT file and pick file
	gROOT.Reset()
	filename = str(path)+"/"+"data-rotated/"+str(answer.absorbermaterial)+"/"+str(answer.namefile)+"/B4.root"
	inputFile = TFile(str(filename)) #root input file
	tree = TTree() #Tree name B4
	inputFile.GetObject("B4", tree)

	#Set Parameters, containers and counters
	traditionalreconstructedenergy = []
	machinelearningreconstructedenergy = []
	secondmlreconstructedenergy = []
	traditionalscinreconstructedenergy = []
	machinelearningscinreconstructedenergy = []
	secondmlscinreconstructedenergy = []
	traditionalcherreconstructedenergy = []
	machinelearningcherreconstructedenergy = []
	secondmlcherreconstructedenergy = []
	minEnergy = []
	correctML_ID = 0
	wrongML_ID = 0
	trueemfraction = []
	dremfraction = []
	differenceemfraction = []
	totalenergy = []
	energyscin = []
	trueenergy = []
	countererror = 0
	CsignalE = []
	SsignalE = []
	MLFEM = 0
	mlemfraction = []
	differenceMLemfraction = []

	trainedlist = machinelearning.unpicklelistofclusters(answer.trainednamefile) #unpickle list of trained clusters in aswer

	for Event in range(int(answer.NofEvents)):

		tree.GetEntry(Event)

		#Set values of tree
		PrimaryParticleName = tree.PrimaryParticleName # MC truth: primary particle Geant4 name
		EnergyEscaped = tree.EscapedEnergy #Escaped energy
		PrimaryParticleEnergy = tree.PrimaryParticleEnergy - EnergyEscaped # MC truth: primary particle energy
		EnergyTot = tree.EnergyTot # Total energy deposited in calorimeter
		Energyem = tree.Energyem # Energy deposited by the em component
		EnergyScin = tree.EnergyScin # Energy deposited in Scin fibers (not Birk corrected)
		EnergyCher = tree.EnergyCher # Energy deposited in Cher fibers (not Birk corrected)
		NofCherenkovDetected = tree.NofCherenkovDetected # Total Cher p.e. detected
		VectorSignals = tree.VectorSignals # Vector of energy deposited in Scin fibers (Birk corrected)
		VectorSignalsCher = tree.VectorSignalsCher # Vector of Cher p.e. detected in Cher fibers

		#if doing training with jets primary particle energy is energy of a single track
		#I have to fix it
		if "jet" in answer.namefile: #if it is a jet file
			b = [s for s in answer.namefile[0:3] if s.isdigit()]
			c = ""
			for s in range(len(b)):
				c = c+b[s]
			d = int(c)
			PrimaryParticleEnergy = float(d*1000 - EnergyEscaped)
		#end of finxing jet energy problem

		FirstPrimaryParticleEnergy = tree.PrimaryParticleEnergy

		#fixing jet problem part 2
		if "jet" in answer.namefile:
			FirstPrimaryParticleEnergy = float(d*1000)

		if Event % 100 == 0:
			print "--------------------------------------------------------\n"
			print "Processing event " + str(Event) + " of " + str(NofEventsProcessed) + ": " + PrimaryParticleName + " energy " + str(PrimaryParticleEnergy) + " MeV" +"\n"

		#Set analysis and event parameters and create root histograms
		ScinTreshold = 0.03 #30 KeV, energy treshold for single photoelectron production
		CherTreshold = 0 #We assume to detect single Cherenkov p.e.
		ScinMaxFiber = max(list(VectorSignals)) #maximum energy deposit in scintillating fibers
		CherMaxFiber = max(list(VectorSignalsCher)) #maximum Cherenkov p.e. in clear fibers
		NumFibers = len(list(VectorSignals)) #number of fibers (Cherenkov+Scintillating) 
		NumModules = len(list(VectorSignals)) #number of modules

		#This version of the code deals with single particle events
		event = eventering.buildevent(VectorSignals, VectorSignalsCher)

		#Compute event informations (Energy Scin, Energy Cher)
		if event.chersignal == 0.0:
			print "Event with no Cherenkov signal"
			countererror = countererror +1
			continue
		if PrimaryParticleEnergy == 0.0:
			print "Event with wrong primary particle energy"
			countererror = countererror +1
			continue

		if answer.answerdualreadout == "y":
			event.computeEnergyScin(scincalibconstant) # Energy reconstructed with scin signal: scinsignal*calibconstantscin
			event.computeEnergyCher(chercalibconstant) # Energy reconstructed with cher signal: chersignal*calibconstantcher
			traditionalscinreconstructedenergy.append(event.EnergyScin)
			traditionalcherreconstructedenergy.append(event.EnergyCher)
		
			#FIND FAKE ID TO PERFORM TRADITIONAL CALORIMETRIC MEASUREMENTS
			#Perform traditional calorimetric measurement
			#Must be done after computeEnergyScin and computeEnergyCher
			if particletype == "electron": #here don't use particle type identification
				event.ID = "electron"
			else:
				event.ID = "hadron"
			truth_ID = event.ID
			event.computeEnergy(Chi)
			traditionalreconstructedenergy.append(event.Energy)
			if Event % 100 == 0:
				print "Traditional recon energy = " + str(event.Energy) 
		if answer.answermachinelearning == "y":
			'''
			#Find event ID with machine learning
			event.compute_ID(answer.trainednamefile, answer.Noftrainedevents) # "trainednamefile" is name of trained events

			#Compute correct and wrong ML IDs #still not used to perform particle ID capability
			if event.ID == truth_ID:
				correctML_ID += 1
			elif event.ID != truth_ID:
				wrongML_ID += 1
			'''
			#Find cluster ML energy with machine learning
			event.ID = str(particletype)
			MLFEM = event.compute_MLEnergy(trainedlist, answer.Noftrainedevents) # "namefile" is name of trained events
			machinelearningreconstructedenergy.append(event.MLEnergy) #not corrected for containment
			machinelearningscinreconstructedenergy.append(event.MLScinEnergy)
			machinelearningcherreconstructedenergy.append(event.MLCherEnergy)
			secondmlreconstructedenergy.append(event.secondMLEnergy) #not corrected for containment
			secondmlscinreconstructedenergy.append(event.secondMLScinEnergy) 
			secondmlcherreconstructedenergy.append(event.secondMLCherEnergy)
		
			if Event % 100 == 0:
				print "ML energy = " + str(event.MLEnergy)
				print "Second ML energy = " + str(event.secondMLEnergy)

		if answer.answerminimization == "y":
			#Find energy with parameters
			event.computeminEnergy(a, b)
			minEnergy.append(event.minEnergy) #not corrected for containment

			if Event % 100 == 0:
				print "Min Energy = " + str(event.minEnergy)
		
		if answer.answerdualreadout == "y":
			#DR Em fraction analysis
			trueemfraction.append(Energyem/PrimaryParticleEnergy) # Energy deposited by electrons and positrons / totalenergy in calorimeter
			e = traditionalcherreconstructedenergy[Event-countererror]/traditionalscinreconstructedenergy[Event-countererror] #used in line below
			dremfraction.append((-HEcher+e*HEscin)/(1-HEcher-e+e*HEscin)) # DR equation for em fraction
			differenceemfraction.append((trueemfraction[Event-countererror]-dremfraction[Event-countererror])/trueemfraction[Event-countererror])
			
			#ML Em fraction analysis
			mlemfraction.append(MLFEM)
			differenceMLemfraction.append((trueemfraction[Event-countererror]-mlemfraction[Event-countererror])/trueemfraction[Event-countererror])

		#Total energy, energy in scintillator and MC truth energy analysis
		totalenergy.append(EnergyTot)
		energyscin.append(EnergyScin) #to check if channeling is present
		trueenergy.append(PrimaryParticleEnergy)

		#Create C/E and S/E histograms
		if PrimaryParticleEnergy != 0:
			CsignalE.append(sum(VectorSignalsCher)/PrimaryParticleEnergy)
			SsignalE.append(sum(VectorSignals)/PrimaryParticleEnergy)

	#End of for loop

	if answer.answerdualreadout == "y" and answer.answermachinelearning == "y":
		#Combine DR Energy and 	ML Energy
		combinedenergy = []
		for Event in range(len(machinelearningreconstructedenergy)):
			combinedenergy.append((machinelearningreconstructedenergy[Event]+traditionalreconstructedenergy[Event])/2)

	if answer.answerdualreadout == "y" and answer.answermachinelearning == "y":
		#Combine DR Energy and ML Energy
		secondcombinedenergy = []
		for Event in range(len(secondmlreconstructedenergy)):
			secondcombinedenergy.append((secondmlreconstructedenergy[Event]+traditionalreconstructedenergy[Event])/2)

	if answer.answerdualreadout == "y" and answer.answermachinelearning == "y":
		#Print ROOT histograms of DR Energy anf ML Energy combined
		ROOTHistograms.create_thirdroothistograms(traditionalscinreconstructedenergy, machinelearningscinreconstructedenergy)
		ROOTHistograms.create_fourthroothistograms(traditionalcherreconstructedenergy, machinelearningcherreconstructedenergy)
		drmean, drmeanerror, drsigma, drsigmaerror, mlmean, mlmeanerror, mlsigma, mlsigmaerror = ROOTHistograms.create_secondroothistograms(traditionalreconstructedenergy, machinelearningreconstructedenergy)
		ROOTHistograms.create_fifthroothistograms(combinedenergy)
		ROOTHistograms.create_sixthroothistograms(secondcombinedenergy)

	if answer.answermachinelearning == "y":
		ROOTHistograms.create_seventhroothistograms(secondmlreconstructedenergy, secondmlscinreconstructedenergy, secondmlcherreconstructedenergy)
		ROOTHistograms.create_eigthroothistogram(machinelearningscinreconstructedenergy, machinelearningcherreconstructedenergy)
	
	if answer.answerdualreadout == "y":
		ROOTHistograms.create_ninethroothistogram(traditionalscinreconstructedenergy, traditionalcherreconstructedenergy)
	
	if answer.answerdualreadout == "y":
		
		#Print ROOT histograms of EM fraction DR and difference
		ROOTHistograms.create_fastroothistograms(trueemfraction, "True Fem", "Fem", "# Events", "TrueFem.pdf")
		ROOTHistograms.create_fastroothistograms(dremfraction, "Dr Fem", "Fem", "# Events", "DRFem.pdf")
		ROOTHistograms.create_fastroothistograms(differenceemfraction, "Difference Fem (DR)", "Fem", "# Events", "DiffDRFem.pdf")
		ROOTHistograms.create_emscatterplot(trueemfraction, differenceemfraction)
		
		#Print ROOT histograms of EM fraction ML and difference
		ROOTHistograms.create_fastroothistograms(mlemfraction, "ML Fem", "Fem", "# Events", "MLFem.pdf")
		ROOTHistograms.create_fastroothistograms(differenceMLemfraction, "Difference Fem (ML)", "Fem", "# Events", "DiffMLFem.pdf")
		ROOTHistograms.create_mlemscatterplot(trueemfraction, differenceMLemfraction)		

	#Print ROOT histograms of total energy, energy in scintillator and MC truth energy
	ROOTHistograms.create_fastroothistograms(totalenergy,"En Tot", "Energy (MeV)", "# Events", "EnergyTot.pdf")
	ROOTHistograms.create_fastroothistograms(energyscin,"En Scin", "Energy (MeV)", "# Events", "EnergyScin.pdf")
	ROOTHistograms.create_fastroothistograms(trueenergy,"MCt En", "Energy (MeV)", "# Events", "EnergyMCtruth.pdf")
	#Print C and S signals over E
	ROOTHistograms.create_fastroothistograms(CsignalE,"c/E", "s/E", "# Events", "CoverE.pdf")
	ROOTHistograms.create_fastroothistograms(SsignalE,"s/E", "c/E", "# Events", "SoverE.pdf")

	if answer.answerminimization == "y":
		#Print ROOT histograms of min energy (energy from minimization)
		minmean, minmeanerror, minsigma, minsigmaerror = ROOTHistograms.create_minenergyroothistograms(minEnergy,"Min En", "Energy (MeV)", "# Events", "minEnergy.pdf")

	
	#Print ML IDs counters
	print str(NofEventsProcessed) + " Events: " + str(correctML_ID) + " correct ID " + str(wrongML_ID) + " wrong ID."

	'''
	#Compute correlation and print into file
	print "Now computing correlations"
	correlationDRML = service.compute_correlation(traditionalreconstructedenergy, machinelearningreconstructedenergy)
	correlationScinCher = service.compute_correlation(traditionalscinreconstructedenergy, traditionalcherreconstructedenergy)
	correlationMLScinCher = service.compute_correlation(machinelearningscinreconstructedenergy, machinelearningcherreconstructedenergy)
	outputfile = open("fastresult.txt","w+")
	outputfile.write(str(PrimaryParticleEnergy) + " MeV\n")
	outputfile.write(str(PrimaryParticleEnergy) + " MeV\n")
	outputfile.write("Correlation DR Energy - ML Energy = " + str(correlationDRML) + " \n") 
	outputfile.write("Correlation DR Scin - DR Cher = " + str(correlationScinCher) + " \n")
	outputfile.write("Correlation ML Scin - ML Cher = " + str(correlationMLScinCher) + " \n")
	outputfile.write("Wrong and correct ID events from machine learning: \n")
	outputfile.write("correct ID: " + str(correctML_ID) + " \n")
	outputfile.write("wrong ID: " + str(wrongML_ID) + " \n")
	outputfile.close()
	'''

	#Create lego plot for first event
	tree.GetEntry(0)

	#Set values of tree
	PrimaryParticleName = tree.PrimaryParticleName # MC truth: primary particle Geant4 name
	PrimaryParticleEnergy = hadconteinment*tree.PrimaryParticleEnergy # MC truth: primary particle energy, 99% energy containment
	EnergyTot = tree.EnergyTot # Total energy deposited in calorimeter
	Energyem = tree.Energyem # Energy deposited by the em component
	EnergyScin = tree.EnergyScin # Energy deposited in Scin fibers (not Birk corrected)
	EnergyCher = tree.EnergyCher # Energy deposited in Cher fibers (not Birk corrected)
	NofCherenkovDetected = tree.NofCherenkovDetected # Total Cher p.e. detected
	VectorSignals = tree.VectorSignals # Vector of energy deposited in Scin fibers (Birk corrected)
	VectorSignalsCher = tree.VectorSignalsCher # Vector of Cher p.e. detected in Cher fibers

	#Create grouped vectors (1.2 x 1.2 cm^2)
	GroupedVectorSignals = mapgroup.group(VectorSignals)
	GroupedVectorSignalsCher = mapgroup.group(VectorSignalsCher)

	print "LEGO plot: " + PrimaryParticleName + " energy " + str(PrimaryParticleEnergy) + " MeV" +"\n"

	#Set analysis and event parameters and create root histograms
	ScinTreshold = 0.03 #30 KeV, energy treshold for single photoelectron production
	CherTreshold = 0 #We assume to detect single Cherenkov p.e.
	ScinMaxFiber = max(list(VectorSignals)) #maximum energy deposit in scintillating fibers
	CherMaxFiber = max(list(VectorSignalsCher)) #maximum Cherenkov p.e. in clear fibers
	NumFibers = len(list(VectorSignals)) #number of fibers (Cherenkov+Scintillating) 
	NumModules = len(list(GroupedVectorSignals)) #number of modules

	#Show lego plots ("images") of one event
	ROOTHistograms.create_firstroothistograms(PrimaryParticleName, VectorSignals, VectorSignalsCher, GroupedVectorSignals, 
	GroupedVectorSignalsCher, ScinTreshold, CherTreshold, ScinMaxFiber, CherMaxFiber, NumFibers, NumModules)
	
	gSystem.ProcessEvents()

	#Create folder result and move output file into
	if not os.path.exists(str(path)+"/"+"fastresult"):
		os.makedirs(str(path)+"/"+"fastresult")
	if not os.path.exists(str(path)+"/"+"fastresult/"+str(answer.absorbermaterial)):
		os.makedirs(str(path)+"/"+"fastresult/"+str(answer.absorbermaterial))
	if not os.path.exists(str(path)+"/"+"fastcalibration/"+str(answer.absorbermaterial)+"/"+str(endtime.strftime("%H-%M_%d_%m_%Y"))):
		os.makedirs(str(path)+"/"+"fastcalibration/"+str(answer.absorbermaterial)+"/"+str(endtime.strftime("%H-%M_%d_%m_%Y")))
	foldername = str(path)+"/"+"fastresult/"+str(answer.absorbermaterial)+"/"+str(endtime.strftime("%H-%M_%d_%m_%Y"))+"/"+answer.foldername # folder name asked at the beginning
	os.makedirs(str(foldername))
	for pdf_file in glob.iglob('*.pdf'):
		shutil.move(pdf_file, str(foldername))
	for txt_file in glob.iglob('*.txt'):
		shutil.move(txt_file, str(foldername))

	return drmean, drmeanerror, drsigma, drsigmaerror, (FirstPrimaryParticleEnergy*hadconteinment)/1000, mlmean, mlmeanerror, mlsigma, mlsigmaerror, minmean, minmeanerror, minsigma, minsigmaerror
#----------------------------------------------------------------------------------------------------
def resolution(nofpoints, vectorenergy, vectormean, vectormeanerror, vectorsigma, vectorsigmaerror, graphname):

	#Set parameters
	energy = 0
	mean = 0
	meanerror = 0
	sigma = 0
	sigmaerror = 0
	nofevents = 0
	energies = array('d')
	energieserrors = array('d')
	means = array('d')
	sigmasmeans = array('d')
	sigmasmeanserrors = array('d')
	a = 0
	b = 0
	mean_trueenergy = array('d')
	mean_trueenergyerrors = array('d')

	#For loop to fill graph
	for point in range(int(nofpoints)):
		energies.append(float(vectorenergy[point]))
		energieserrors.append(0.0)
		means.append(float(vectormean[point]))
		#mean_trueenergy.append(float(mean)/float(energy))
		#mean_trueenergyerrors.append((float(meanerror)/float(nofevents)**0.5)*(1/float(mean))*(float(mean)/float(energy)))
		sigmasmeans.append(float(vectorsigma[point])/float(vectormean[point]))
		sigmasmeanserrors.append((float(vectormeanerror[point])/float(vectormean[point])+float(vectorsigmaerror[point])/float(vectorsigma[point]))*(float(vectorsigma[point])/float(vectormean[point])))

	a, b = ROOTHistograms.create_resolutiongraph(nofpoints, energies, sigmasmeans, energieserrors, sigmasmeanserrors, graphname)

	print "Energy resolution: " + str(a*100) + "%" + "/sqrt(E) + " + str(b) 

	#Create linearity graph
	#ROOTHistograms.create_linearitygraph(n, energies, mean_trueenergy, energieserrors, mean_trueenergyerrors, graphname2)
	
	#Create folder result and move output file into
	os.chdir(str(path))
	if not os.path.exists(str(path)+"/"+"fastresolution"):
		os.makedirs(str(path)+"/"+"fastresolution")
	if not os.path.exists(str(path)+"/"+"fastresolution/"+str(absorbermaterial)):
		os.makedirs(str(path)+"/"+"fastresolution/"+str(absorbermaterial))
	if not os.path.exists(str(path)+"/"+"fastresolution/"+str(absorbermaterial)+"/"+str(endtime.strftime("%H-%M_%d_%m_%Y"))):
		os.makedirs(str(path)+"/"+"fastresolution/"+str(absorbermaterial)+"/"+str(endtime.strftime("%H-%M_%d_%m_%Y")))
	foldername = str(path)+"/"+"fastresolution/"+str(absorbermaterial)+"/"+str(endtime.strftime("%H-%M_%d_%m_%Y"))
	os.chdir("/home/lorenzo/cernbox/work/Git-to-Mac/AnalysisFullCalorimeter")
	for pdf_file in glob.iglob('*.pdf'):
		shutil.move(pdf_file, str(foldername))
#----------------------------------------------------------------------------------------------------
#Section used for the complete event analysis. Enable this part with Docopt.
if ToDo == "analysis":

	counter = 0

	answers = [] #list of answers
	Files = [] #list of all files to study
	Files = raw_input("List all the files to study: ").split()
	NofFiles = int(len(Files))

	for File in range(NofFiles):

		if counter == 0:
			absorbermaterial = raw_input("Choose absorber material: ")

		#Choose ROOT files
		gROOT.Reset()
		file = Files[int(counter)]
		filename = str(path)+"/"+"data-rotated/"+str(absorbermaterial)+"/"+str(file)+"/B4.root"
		inputFile = TFile(str(filename)) #root input file
		tree = TTree() #Tree name B4
		inputFile.GetObject("B4", tree)

		#Check if MATERIAL is correct, very important to not mix things!!!
		tree.GetEntry(0)
		geant4material = str(tree.AbsorberMaterial)
		print "Found Geant4 absorber material: "+str(geant4material)+" \n"
		geant4material = geant4material[0:len(absorbermaterial)]
		if geant4material != str(absorbermaterial):
			sys.exit()
		
		#Ask for how many events to be processed
		if counter == 0:
			NofEvents = tree.GetEntries()
			NofEventsProcessed = raw_input(str(NofEvents) + " Events, how many to process: ")

		inputFile.Close()

		if counter == 0:
			#Choose what analysis
			answerdualreadout = raw_input("Want to do dual-readout analysis? (y/n) ")
			answermachinelearning = raw_input("Want to do machine learning analysis? (y/n) ")
			answerminimization = raw_input("Want to do minimization analysis? (y/n) ")
			answermincalib = raw_input("Minimization with which calibration constants (20,40,60,80)? ")	

			particletype = raw_input("Input particle type (electron or hadron): ")

			namefile = raw_input("Input name of file for machine learning: ")
			namefile = str(path)+"/"+"fasttrained-rotated/"+str(absorbermaterial)+"/"+str(namefile)
			noftrained = raw_input("How many trained events to use? ")

		#Ask for folder name where to store results
		foldername = str(Files[counter])

		counter = counter + 1

		#save answers
		answer = questions.Analysisanswer(absorbermaterial, file, NofEventsProcessed, noftrained, particletype, str(answerdualreadout), str(answermachinelearning), str(answerminimization), foldername, str(answermincalib), namefile)
		answers.append(answer)

	#Now perform analysis on all files #resolution graphs obtained by default
	drmeanvector = []
	drmeanerrorvector = []
	drsigmavector = []
	drsigmaerrorvector = []
	primaryenergyvector = []
	mlmeanvector = []
	mlmeanerrorvector = []
	mlsigmavector = []
	mlsigmaerrorvector = []
	minmeanvector = []
	minmeanerrorvector = []
	minsigmavector = []
	minsigmaerrorvector = []
	for answer in answers:
		print "Analysis with: "+str(answer.namefile)+" \n"
		drmean, drmeanerror, drsigma, drsigmaerror, primaryenergy, mlmean, mlmeanerror, mlsigma, mlsigmaerror, minmean, minmeanerror, minsigma, minsigmaerror = analysis(answer)
		drmeanvector.append(drmean)
		drmeanerrorvector.append(drmeanerror)
		drsigmavector.append(drsigma)
		drsigmaerrorvector.append(drsigmaerror)
		primaryenergyvector.append(primaryenergy)
		mlmeanvector.append(mlmean)
		mlmeanerrorvector.append(mlmeanerror)
		mlsigmavector.append(mlsigma)
		mlsigmaerrorvector.append(mlsigmaerror)
		minmeanvector.append(minmean)
		minmeanerrorvector.append(minmeanerror)
		minsigmavector.append(minsigma)
		minsigmaerrorvector.append(minsigmaerror)
		

	resolution(int(NofFiles), primaryenergyvector, drmeanvector, drmeanerrorvector, drsigmavector, drsigmaerrorvector, "dualreadout.pdf")
	resolution(int(NofFiles), primaryenergyvector, mlmeanvector, mlmeanerrorvector, mlsigmavector, mlsigmaerrorvector, "machinelearning.pdf")
	resolution(int(NofFiles), primaryenergyvector, minmeanvector, minmeanerrorvector, minsigmavector, minsigmaerrorvector, "minimization.pdf")
#---------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------
if ToDo == "graph":

	nofpoints = raw_input("how many points? ")

	#Set parameters
	x = 0
	y = 0
	xs = array('d')
	ys = array('d')

	#For loop to fill graph
	for point in range(int(nofpoints)):
		x = raw_input("x: ")
		y = raw_input("y: ")
		xs.append(float(x))
		ys.append(float(y))

	ROOTHistograms.create_fastgraph2(xs, ys)
#---------------------------------------------------------------------------------------------------------
#If -do went wrong
if ToDo != "calibration" and ToDo != "training" and ToDo != "analysis" and ToDo != "minimization" and ToDo != "graph":
	print "Pass option: calibration, training, minimizarion, analysis, graph!"
#---------------------------------------------------------------------------------------------------------