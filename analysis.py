import random
import os
import shutil
import glob
from array import array
from ROOT import gROOT, TFile, TTree
import map
import mapgroup
import clustering
import machinelearning
import ROOTHistograms
import service

#Hardcoded calibration constants and Chi factor
#You can estimate them with dedicated part 
scincalibconstant = 23.92 #MeV/MeV
chercalibconstant = 23.00 #MeV/Cpe
Chi = 0.50 #to be better defined
HEcher = 0.360
HEscin = 0.679

#Set root file and tree to be analyzed 
gROOT.Reset()
file = raw_input("Insert namefile: ")
filename = str(file)+"B4.root"
inputFile = TFile(str(filename)) #root input file
tree = TTree() #Tree name B4
inputFile.GetObject("B4", tree)

#---------------------------------------------------------------------------------------------------
#Section used only if you want to estimate the calibration constants and Chi factor of a module.
#DREAM is calibrated with electrons. The input file must be from electron events with the same energy.
#To estimate h/e Cher and scin values and the Chi factor a second root file with pions (+ or -) is needed.

#To later perform h/e Cher and scin computation and Chi factor
print "You have given an electron file event, now pass a pion one for h/e estimation.\n"
namepionfile = raw_input("Insert name of ROOT pion (+ or -) file: ")

#Set parameters
NofEvents = tree.GetEntries()
fastchercalibconstant = 0 #Done without building clusters
fastscincalibconstant = 0 #Done without building clusters
fastChi = 0
chercalibconstant = 0
scincalibconstant = 0
Chi = 0
firstcounter = 0
secondcounter = 0

#Estimate Cherenkov and Scintillation calibration constants
for Event in range(int(NofEvents)):

	tree.GetEntry(Event)

	#Set values of tree
	PrimaryParticleName = tree.PrimaryParticleName # MC truth: primary particle Geant4 name
	PrimaryParticleEnergy = tree.PrimaryParticleEnergy # MC truth: primary particle energy
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

	#For output file
	FirstPrimaryParticleName = PrimaryParticleName
	FirstPrimaryParticleEnergy = PrimaryParticleEnergy

	#Estimating calibration constants using all the calorimeter signal
	if sum(GroupedVectorSignalsCher) == 0.0:
		continue
		print "One event with no Cherenkov signal"

	fastchercalibconstant += PrimaryParticleEnergy/sum(GroupedVectorSignalsCher)
	fastscincalibconstant += PrimaryParticleEnergy/sum(GroupedVectorSignals)

	firstcounter = firstcounter+1

	#Estimating calibration constants using cluster signals
	fixedindex = GroupedVectorSignals.index(max(GroupedVectorSignals))
	cluster = clustering.buildcluster(fixedindex, GroupedVectorSignals)
	cluster.find_modules(GroupedVectorSignals)
	cluster.compute_chersignal(GroupedVectorSignalsCher)

	if cluster.chersignal == 0.0:
		continue
		print "One cluster with no Cherenkov signal!"

	chercalibconstant += PrimaryParticleEnergy/cluster.chersignal
	scincalibconstant += PrimaryParticleEnergy/cluster.scinsignal

	secondcounter = secondcounter+1

#Finalize calibration constants
fastchercalibconstant = fastchercalibconstant/firstcounter
fastscincalibconstant = fastscincalibconstant/firstcounter
chercalibconstant = chercalibconstant/secondcounter
scincalibconstant = scincalibconstant/secondcounter

print "End of calibration constants -> now estimating h/e values."

#Estimating h/e Cher and scin values and Chi factor
#Set root file and tree
gROOT.Reset()
inputFile = TFile(namepionfile) #root input file
tree = TTree()
inputFile.GetObject("B4", tree)

#Set parameters and counters
NofEvents = tree.GetEntries()
Scinenergy1 = 0
Scinenergy2 = 0
Scinenergy3 = 0
Scinenergy4 = 0
Cherenergy1 = 0
Cherenergy2 = 0
Cherenergy3 = 0
Cherenergy4 = 0
FastScinenergy1 = 0
FastScinenergy2 = 0
FastScinenergy3 = 0
FastScinenergy4 = 0
FastCherenergy1 = 0
FastCherenergy2 = 0
FastCherenergy3 = 0
FastCherenergy4 = 0
counter1 = 0
counter2 = 0
counter3 = 0
counter4 = 0

for Event in range(int(NofEvents)):

	tree.GetEntry(Event)

	#Set values of tree
	PrimaryParticleName = tree.PrimaryParticleName # MC truth: primary particle Geant4 name
	PrimaryParticleEnergy = tree.PrimaryParticleEnergy # MC truth: primary particle energy
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

	#Create histograms
	Fem = Energyem/PrimaryParticleEnergy
	if Fem < 0.25 and Fem > 0.0:
		#Fast Part
		FastScinenergy1 += fastscincalibconstant*sum(GroupedVectorSignals)
		FastCherenergy1 += fastchercalibconstant*sum(GroupedVectorSignalsCher)
		#Cluster Part
		fixedindex = GroupedVectorSignals.index(max(GroupedVectorSignals))
		cluster = clustering.buildcluster(fixedindex, GroupedVectorSignals)
		if cluster.radius < 20:
			continue
		cluster.find_modules(GroupedVectorSignals)
		cluster.compute_chersignal(GroupedVectorSignalsCher)

		Scinenergy1 += scincalibconstant*cluster.scinsignal
		Cherenergy1 += chercalibconstant*cluster.chersignal

		counter1 += 1

	elif Fem < 0.50 and Fem > 0.25:
		#Fast Part
		FastScinenergy2 += fastscincalibconstant*sum(GroupedVectorSignals)
		FastCherenergy2 += fastchercalibconstant*sum(GroupedVectorSignalsCher)
		#Cluster Part
		fixedindex = GroupedVectorSignals.index(max(GroupedVectorSignals))
		cluster = clustering.buildcluster(fixedindex, GroupedVectorSignals)
		if cluster.radius < 20:
			continue
		cluster.find_modules(GroupedVectorSignals)
		cluster.compute_chersignal(GroupedVectorSignalsCher)

		Scinenergy2 += scincalibconstant*cluster.scinsignal
		Cherenergy2 += chercalibconstant*cluster.chersignal

		counter2 += 1

	elif Fem < 0.75 and Fem > 0.50:
		#Fast Part
		FastScinenergy3 += fastscincalibconstant*sum(GroupedVectorSignals)
		FastCherenergy3 += fastchercalibconstant*sum(GroupedVectorSignalsCher)
		#Cluster Part
		fixedindex = GroupedVectorSignals.index(max(GroupedVectorSignals))
		cluster = clustering.buildcluster(fixedindex, GroupedVectorSignals)
		if cluster.radius < 20:
			continue
		cluster.find_modules(GroupedVectorSignals)
		cluster.compute_chersignal(GroupedVectorSignalsCher)

		Scinenergy3 += scincalibconstant*cluster.scinsignal
		Cherenergy3 += chercalibconstant*cluster.chersignal

		counter3 += 1

	elif Fem < 1.0 and Fem > 0.75:
		#Fast Part
		FastScinenergy4 += fastscincalibconstant*sum(GroupedVectorSignals)
		FastCherenergy4 += fastchercalibconstant*sum(GroupedVectorSignalsCher)
		#Cluster Part
		fixedindex = GroupedVectorSignals.index(max(GroupedVectorSignals))
		cluster = clustering.buildcluster(fixedindex, GroupedVectorSignals)
		if cluster.radius < 20:
			continue
		cluster.find_modules(GroupedVectorSignals)
		cluster.compute_chersignal(GroupedVectorSignalsCher)

		Scinenergy4 += scincalibconstant*cluster.scinsignal
		Cherenergy4 += chercalibconstant*cluster.chersignal

		counter4 += 1

#Average energy reconstructed
FastScinenergy1 = FastScinenergy1/counter1
FastCherenergy1 = FastCherenergy1/counter1
Scinenergy1 = Scinenergy1/counter1
Cherenergy1 = Cherenergy1/counter1
FastScinenergy2 = FastScinenergy2/counter2
FastCherenergy2 = FastCherenergy2/counter2
Scinenergy2 = Scinenergy2/counter2
Cherenergy2 = Cherenergy2/counter2
FastScinenergy3 = FastScinenergy3/counter3
FastCherenergy3 = FastCherenergy3/counter3
Scinenergy3 = Scinenergy3/counter3
Cherenergy3 = Cherenergy3/counter3
FastScinenergy4 = FastScinenergy4/counter4
FastCherenergy4 = FastCherenergy4/counter4
Scinenergy4 = Scinenergy4/counter4
Cherenergy4 = Cherenergy4/counter4

#Values grouped in list for ROOT graphs
FastScinenergy = array('d',[FastScinenergy1, FastScinenergy2, FastScinenergy3, FastScinenergy4])
FastCherenergy = array('d',[FastCherenergy1, FastCherenergy2, FastCherenergy3, FastCherenergy4])
Scinenergy = array('d',[Scinenergy1, Scinenergy2, Scinenergy3, Scinenergy4])
Cherenergy = array('d',[Cherenergy1, Cherenergy2, Cherenergy3, Cherenergy4])
fem = array('d',[0.25/2, 0.5/2, 0.75/2, 1.0/2])

#Create graphs, compute he Cher and scin values and Chi factor
Fasthescin, Fasthecher, hescin, hecher = ROOTHistograms.create_graph(FastScinenergy, FastCherenergy, Scinenergy, Cherenergy, fem, PrimaryParticleEnergy)
fastChi = (1 - Fasthescin)/(1 - Fasthecher)
Chi = (1 - hescin)/(1 - hecher)

#Print calibration results on file
outputfile = open("calibration.txt","w+")
outputfile.write("Calibration constants for electron + hadron of energies: \n")
outputfile.write(str(FirstPrimaryParticleEnergy) + " MeV\n")
outputfile.write(str(PrimaryParticleEnergy) + " MeV\n")
outputfile.write("-------------------------FAST ANALYSIS---------------------------------------------------------\n")
outputfile.write("fastchercalibconstant = " + str(fastchercalibconstant) + "MeV/Cpe" + " fastscincalibconstant = " + str(fastscincalibconstant) + " MeV/MeV\n")
outputfile.write("fasthescin = " + str(Fasthescin) + " fasthecher = " + str(Fasthecher) + " fastChi = " + str(fastChi) + " \n")
outputfile.write("-------------------------CLUSTER ANALYSIS------------------------------------------------------\n")
outputfile.write("chercalibconstant = " + str(chercalibconstant) + " MeV/Cpe " + "scincalibconstant = " + str(scincalibconstant) + " MeV/MeV\n")
outputfile.write("hescin = " + str(hescin) + " hecher = " + str(hecher) + " Chi = " + str(Chi))
outputfile.close()

#Create folder calibration and move output file into
if not os.path.exists("calibration"):
    os.makedirs("calibration")
foldername = raw_input("Insert folder name (take care folder does not exist!): ")
foldername = "calibration/"+foldername
os.makedirs(str(foldername))
for eps_file in glob.iglob('*.eps'):
    shutil.move(eps_file, str(foldername))
for txt_file in glob.iglob('*.txt'):
	shutil.move(txt_file, str(foldername))
print "End of h/e values estimation."
#---------------------------------------------------------------------------------------------------	
'''
#---------------------------------------------------------------------------------------------------
#Section used only if you want to ceate a trained list of clusters to later perform machine learning
#ID and energy reconstruction. ROOT events must be single particle events to be classified in as "electron"
#or "hadron" with a given MC truth energy. Enable this part with Docopt.

#Set name of pickle file containing trained clusters
namefile = raw_input("Input pickle namefile (with .p extension): ")
particletype = raw_input("Input particle type (electron or hadron): ")

#Ask for how many events to be processed
NofEvents = tree.GetEntries()
NofEventsProcessed = raw_input(str(NofEvents) + "Events,how many to process: ")

for Event in range(int(NofEventsProcessed)):
	
	tree.GetEntry(Event)

	#Set values of tree
	PrimaryParticleName = tree.PrimaryParticleName # MC truth: primary particle Geant4 name
	PrimaryParticleEnergy = tree.PrimaryParticleEnergy # MC truth: primary particle energy
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

	print "--------------------------------------------------------\n"
	print "Processing event " + str(Event) + ": " + PrimaryParticleName + " energy " + str(PrimaryParticleEnergy) + " MeV" +"\n"

	#Create cluster
	fixedindex = GroupedVectorSignals.index(max(GroupedVectorSignals))
	cluster = clustering.buildcluster(fixedindex, GroupedVectorSignals)
	cluster.find_modules(GroupedVectorSignals)
	cluster.compute_chersignal(GroupedVectorSignalsCher)
	cluster.computeEnergyScin(scincalibconstant)
	cluster.computeEnergyCher(chercalibconstant)
	if cluster.chersignal == 0.0:
		print "One cluster with no Cherenkov signal"
		continue
	cluster.Energy = PrimaryParticleEnergy
	if particletype == "electron":
		cluster.ID = "electron"
	else:
		cluster.ID = "hadron"

	#Creat list of trained clusters, random shuffle it and save in pickle file
	listoftrainedclusters.append(cluster)

random.shuffle(listoftrainedclusters)
machinelearning.picklelistofclusters(listoftrainedclusters, namefile)
print "Pickled " + str(len(listoftrainedclusters)) + " trained clusters in " + namefile +"."
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
#Section used for the complete cluster analysis. Enable this part with Docopt.

#Set name of pickle file containing trained clusters
namefile = raw_input("Input name of file for machine learning: ")
particletype = raw_input("Input particle type (electron or hadron): ")

#Ask for how many events to be processed
NofEvents = tree.GetEntries()
NofEventsProcessed = raw_input(str(NofEvents) + "Events,how many to process: ")

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
correctML_ID = 0
wrongML_ID = 0
trueemfraction = []
dremfraction = []
differenceemfraction = []
totalenergy = []
energyscin = []
clusterradius = []

for Event in range(int(NofEventsProcessed)):

	tree.GetEntry(Event)

	#Set values of tree
	PrimaryParticleName = tree.PrimaryParticleName # MC truth: primary particle Geant4 name
	PrimaryParticleEnergy = tree.PrimaryParticleEnergy # MC truth: primary particle energy
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

	print "--------------------------------------------------------\n"
	print "Processing event " + str(Event) + " of " + str(NofEventsProcessed) + ": " + PrimaryParticleName + " energy " + str(PrimaryParticleEnergy) + " MeV" +"\n"

	#Set analysis and event parameters and create root histograms
	ScinTreshold = 0.03 #30 KeV, energy treshold for single photoelectron production
	CherTreshold = 0 #We assume to detect single Cherenkov p.e.
	ScinMaxFiber = max(list(VectorSignals)) #maximum energy deposit in scintillating fibers
	CherMaxFiber = max(list(VectorSignalsCher)) #maximum Cherenkov p.e. in clear fibers
	NumFibers = len(list(VectorSignals)) #number of fibers (Cherenkov+Scintillating) 
	NumModules = len(list(GroupedVectorSignals)) #number of modules

	#This version of the code deals with single particle events
	#Find a single maximum and build a cluster around it
	fixedindex = GroupedVectorSignals.index(max(GroupedVectorSignals))
	cluster = clustering.buildcluster(fixedindex, GroupedVectorSignals)
	#Fill cluster radius container
	clusterradius.append(cluster.radius)

	#Compute cluster informations (Cherenkov signal, Energy Scin, Energy Cher)
	cluster.find_modules(GroupedVectorSignals)
	cluster.compute_chersignal(GroupedVectorSignalsCher)
	if cluster.chersignal == 0.0:
		print "Cluster with no Cherenkov signal"
		continue
	cluster.computeEnergyScin(scincalibconstant) # Energy reconstructed with scin signal: scinsignal*calibconstantscin
	cluster.computeEnergyCher(chercalibconstant) # Energy reconstructed with cher signal: chersignal*calibconstantcher
	traditionalscinreconstructedenergy.append(cluster.EnergyScin)
	traditionalcherreconstructedenergy.append(cluster.EnergyCher)
	
	#FIND FAKE ID TO PERFORM TRADITIONAL CALORIMETRIC MEASUREMENTS
	#Perform traditional calorimetric measurement
	#Must be done after compute_chersignal, computeEnergyScin and computeEnergyCher
	if particletype == "electron":
		cluster.ID = "electron"
	else:
		cluster.ID = "hadron"
	truth_ID = cluster.ID
	cluster.computeEnergy(Chi)
	traditionalreconstructedenergy.append(cluster.Energy)
	print "Traditional energy = " + str(cluster.Energy) 

	#Find cluster ID with machine learning
	cluster.compute_ID(namefile) # "namefile" is name of trained events

	#Compute correct and wrong ML IDs
	if cluster.ID == truth_ID:
		correctML_ID += 1
	elif cluster.ID != truth_ID:
		wrongML_ID += 1

	#Find cluster ML energy with machine learning
	#Must be done after cluster.compute_ID
	#Here I assign by default the correct ID
	cluster.ID = str(particletype)
	cluster.compute_MLEnergy(namefile) # "namefile" is name of trained events
	machinelearningreconstructedenergy.append(cluster.MLEnergy)
	machinelearningscinreconstructedenergy.append(cluster.MLScinEnergy)
	machinelearningcherreconstructedenergy.append(cluster.MLCherEnergy)
	secondmlreconstructedenergy.append(cluster.secondMLEnergy)
	secondmlscinreconstructedenergy.append(cluster.secondMLScinEnergy)
	secondmlcherreconstructedenergy.append(cluster.secondMLCherEnergy)
	print "ML energy = " + str(cluster.MLEnergy)
	print "Second ML energy = " + str(cluster.secondMLEnergy)


	#Em fraction analysis
	trueemfraction.append(Energyem/PrimaryParticleEnergy) # Energy deposited by electrons and positrons / totalenergy
	e = traditionalcherreconstructedenergy[Event]/traditionalscinreconstructedenergy[Event]
	dremfraction.append((-e*HEcher+HEscin)/(e-e*HEcher-1+HEscin)) # DR equation for em fraction
	differenceemfraction.append(trueemfraction[Event]-dremfraction[Event])

	#Total energy and energy in scintillator analysis
	totalenergy.append(EnergyTot)
	energyscin.append(EnergyScin)

#End of machine learning analysis	
print "End of " + str(NofEventsProcessed) + " events analysis"

#Combine DR Energy and 	ML Energy
combinedenergy = []
for Event in range(len(machinelearningreconstructedenergy)):
	combinedenergy.append((machinelearningreconstructedenergy[Event]+traditionalreconstructedenergy[Event])/2)

#Combine DR Energy and second ML Energy
secondcombinedenergy = []
for Event in range(len(secondmlreconstructedenergy)):
	secondcombinedenergy.append((secondmlreconstructedenergy[Event]+traditionalreconstructedenergy[Event])/2)

#Print ROOT histograms of DR Energy anf ML Energy, Scin energy and ML Scin energy, Cher energy and ML Cher energy
ROOTHistograms.create_thirdroothistograms(traditionalscinreconstructedenergy, machinelearningscinreconstructedenergy)
ROOTHistograms.create_fourthroothistograms(traditionalcherreconstructedenergy, machinelearningcherreconstructedenergy)
ROOTHistograms.create_secondroothistograms(traditionalreconstructedenergy, machinelearningreconstructedenergy)
ROOTHistograms.create_fifthroothistograms(combinedenergy)
ROOTHistograms.create_sixthroothistograms(secondcombinedenergy)
ROOTHistograms.create_seventhroothistograms(secondmlreconstructedenergy, secondmlscinreconstructedenergy, secondmlcherreconstructedenergy)
ROOTHistograms.create_eigthroothistogram(machinelearningscinreconstructedenergy, machinelearningcherreconstructedenergy)
ROOTHistograms.create_ninethroothistogram(traditionalscinreconstructedenergy, traditionalcherreconstructedenergy)

#Print ROOT histograms of EM fraction DR, ML and difference
ROOTHistograms.create_fastroothistograms(trueemfraction, "True Fem", "Fem", "# Events", "TrueFem.eps")
ROOTHistograms.create_fastroothistograms(dremfraction, "Dr Fem", "Fem", "# Events", "DRFem.eps")
ROOTHistograms.create_fastroothistograms(differenceemfraction, "Difference Fem", "Fem", "# Events", "DiffFem.eps")

#Print ROOT histograms of total energy and energy in scintillator
ROOTHistograms.create_fastroothistograms(totalenergy,"En Tot", "Energy (MeV)", "# Events", "EnergyTot.eps")
ROOTHistograms.create_fastroothistograms(energyscin,"En Scin", "Energy (MeV)", "# Events", "EnergyScin.eps")

#Print ROOTHistogram of cluster radius
ROOTHistograms.create_fastroothistograms(clusterradius,"Radius", "cm", "# Events", "Radius.eps")

#Print ML IDs counters
print str(NofEventsProcessed) + " Events: " + str(correctML_ID) + " correct ID " + str(wrongML_ID) + " wrong ID."

#Compute correlation and print into file
print "Now computing correlations"
correlationDRML = service.compute_correlation(traditionalreconstructedenergy, machinelearningreconstructedenergy)
correlationScinCher = service.compute_correlation(traditionalscinreconstructedenergy, traditionalcherreconstructedenergy)
correlationMLScinCher = service.compute_correlation(machinelearningscinreconstructedenergy, machinelearningcherreconstructedenergy)
outputfile = open("result.txt","w+")
outputfile.write(str(PrimaryParticleEnergy) + " MeV\n")
outputfile.write(str(PrimaryParticleEnergy) + " MeV\n")
outputfile.write("Correlation DR Energy - ML Energy = " + str(correlationDRML) + " \n") 
outputfile.write("Correlation DR Scin - DR Cher = " + str(correlationScinCher) + " \n")
outputfile.write("Correlation ML Scin - ML Cher = " + str(correlationMLScinCher) + " \n")
outputfile.write("Wrong and correct ID events from machine learning: \n")
outputfile.write("correct ID: " + str(correctML_ID) + " \n")
outputfile.write("wrong ID: " + str(wrongML_ID) + " \n")
outputfile.close()

#Create lego plot for first event
tree.GetEntry(0)

#Set values of tree
PrimaryParticleName = tree.PrimaryParticleName # MC truth: primary particle Geant4 name
PrimaryParticleEnergy = tree.PrimaryParticleEnergy # MC truth: primary particle energy
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

print "--------------------------------------------------------\n"
print "Processing event: " + PrimaryParticleName + " energy " + str(PrimaryParticleEnergy) + " MeV" +"\n"

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

#Create folder result and move output file into
if not os.path.exists("result"):
    os.makedirs("result")
foldername = raw_input("Insert folder name (take care folder does not exist yet!): ")
foldername = "result/"+foldername
os.makedirs(str(foldername))
for eps_file in glob.iglob('*.eps'):
    shutil.move(eps_file, str(foldername))
for txt_file in glob.iglob('*.txt'):
	shutil.move(txt_file, str(foldername))
#---------------------------------------------------------------------------------------------------------
'''