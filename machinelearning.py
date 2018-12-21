import clustering
import cPickle as pickle
import glob

def picklesinglecluster(vectorsignals, chervectorsignals, ID, MCEnergy, namefile):
	"""Function to pickle a cluster, vector signals scintillation 
	and cherenkov must be provided. Provide also known ID, MC truth energy and 
	name of pickle file"""
	startindex = vectorsignals.index(max(vectorsignals))
	cluster = clustering.buildcluster(startindex, vectorsignals)
	cluster.find_modules(vectorsignals)
	cluster.compute_chersignal(chervectorsignals)
	cluster.ID = str(ID)
	cluster.Energy = MCEnergy
	print(cluster.modules, cluster.chersignal, cluster.ID)
	file = open(str(namefile),"wb")
	pickle.dump(cluster, file)
	file.close()

def unpicklesinglecluster(namefile):
	"""Unpickle signle cluster"""
	file = open(str(namefile), "rb")
	cluster = pickle.load(file)
	print(cluster.modules, cluster.chersignal, cluster.ID)

def picklelistofclusters(listofclusters, namefile):
	"""Function to pickle list of clusters trained for 
	machine learning"""
	file = open(str(namefile),"wb")
	pickle.dump(listofclusters, file)
	file.close()

def unpicklelistofclusters(namefile):
	"""Function to unpickle list of clusters trained for
	machine learning"""
	file = open(str(namefile), "rb")
	listofclusters = pickle.load(file)
	return listofclusters

def differencebetweenclusters(cluster, trainedcluster):
	"""Function to compute difference between clusters"""
	difference1 = (trainedcluster.chersignal/trainedcluster.scinsignal) - (cluster.chersignal/cluster.scinsignal)
	#difference2 = trainedcluster.radius - cluster.radius
	return abs(difference1) # Not using (+ difference2) #return absolute value of difference

#def seconddifferencebetweenclusters(cluster, trainedcluster): # Not used -> going to cancel it
#	"""Function to perform difference on Energy in scintillator"""
#	differece = trainedcluster.EnergyScin - cluster.EnergyScin
#	return abs(differece) #return absolute value of difference

def find_ID(cluster, namefile, noftrained):
	"""Function to find cluster ID with machine learning"""
	trainedlistofclusters = unpicklelistofclusters(namefile)
	listofdifferences = []
	electroncounter = 0
	hadroncounter = 0
	for trainedcluster in trainedlistofclusters:
		listofdifferences.append(differencebetweenclusters(cluster, trainedcluster))
	trainedlistofclusters.sort(key=dict(zip(trainedlistofclusters, listofdifferences)).get) #sort list of trained cluster according to list of differences
	del trainedlistofclusters[0]
	for i in range(int(noftrained)):
		if trainedlistofclusters[i].ID == "electron":
			electroncounter = electroncounter + 1
		else:
			hadroncounter = hadroncounter + 1
	if electroncounter >= noftrained/2:
		cluster.ID = "electron"
	else:
		cluster.ID = "hadron"
	return cluster.ID

def find_energy(cluster, trainedlist, noftrained):
	"""Function to compute by machine learning the energy of 
	a cluster, the cluster ID must me prior computed"""
	selectedlistofclusters = []
	listofdifferences = []
	secondlistofdifferences = []
	chercalibconstant = 0
	scincalibconstant = 0
	secondchercalibconstant = 0
	secondscincalibconstant = 0
	counter = 0
	MLFEM = 0
	weights = []
	''' I don't want this part on ID identification, going to cancel it
	for index in range(len(trainedlistofclusters)):
		if trainedlistofclusters[index].ID == cluster.ID:
			selectedlistofclusters.append(trainedlistofclusters[index]) # selected list of clusters with same ID of cluster to analyze
	'''
	selectedlistofclusters = trainedlist
	for trainedcluster in selectedlistofclusters:
		listofdifferences.append(differencebetweenclusters(cluster, trainedcluster))
	selectedlistofclusters.sort(key=dict(zip(selectedlistofclusters, listofdifferences)).get)
	del selectedlistofclusters[0] #cancel first element (equal to one under reconstruction)
	#print "Doing machine learning with trainedcluster C/S: " + str(selectedlistofclusters[0].scinsignal/selectedlistofclusters[0].chersignal) + "\n"
	#print "                       with cluster C/S: " + str(cluster.scinsignal/cluster.chersignal) + "\n"
	#print "                       abs difference: " + str(abs(selectedlistofclusters[0].scinsignal/selectedlistofclusters[0].chersignal - cluster.scinsignal/cluster.chersignal))
	#print listofdifferences[0:10]
	#listofdifferences.sort()
	#print listofdifferences[0:10]
	#for i in range(int(noftrained)):
		#weights for weighted sum
		#weights.append((1/(cluster.chersignal/cluster.scinsignal - selectedlistofclusters[i].chersignal/selectedlistofclusters[i].scinsignal))**2)
	for i in range(int(noftrained)):
		#chercalibconstant += (selectedlistofclusters[i].Energy/selectedlistofclusters[i].chersignal)*weights[i] #weighted sum
		#scincalibconstant += (selectedlistofclusters[i].Energy/selectedlistofclusters[i].scinsignal)*weights[i]
		#MLFEM += selectedlistofclusters[i].FEM*weights[i] #weighted sum
		chercalibconstant += (selectedlistofclusters[i].Energy/selectedlistofclusters[i].chersignal)
		scincalibconstant += (selectedlistofclusters[i].Energy/selectedlistofclusters[i].scinsignal)
		MLFEM += selectedlistofclusters[i].FEM
		#print "doing ML with energy cluster: " + str(selectedlistofclusters[i].Energy) + " \n"
	#chercalibconstant = chercalibconstant/sum(weights)
	#scincalibconstant = scincalibconstant/sum(weights)
	#MLFEM = MLFEM/sum(weights)
	chercalibconstant = chercalibconstant/int(noftrained)
	scincalibconstant = scincalibconstant/int(noftrained)
	MLFEM = MLFEM/int(noftrained)
	cluster.MLEnergy = (cluster.chersignal*chercalibconstant + cluster.scinsignal*scincalibconstant)/2
	cluster.MLScinEnergy = cluster.scinsignal*scincalibconstant
	cluster.MLCherEnergy = cluster.chersignal*chercalibconstant
	#Now find second ML energiess
	secondlistofclusters = selectedlistofclusters[0:int(float(noftrained)/2)] # half events
	weights = weights[0:(int(noftrained)/2)] #half weights
	for i in range(int(float(noftrained)/2)):
		#secondchercalibconstant += (secondlistofclusters[i].Energy/secondlistofclusters[i].chersignal)*weights[i]
		#secondscincalibconstant += (secondlistofclusters[i].Energy/secondlistofclusters[i].scinsignal)*weights[i]
		secondchercalibconstant += (secondlistofclusters[i].Energy/secondlistofclusters[i].chersignal)
		secondscincalibconstant += (secondlistofclusters[i].Energy/secondlistofclusters[i].scinsignal)
	#secondchercalibconstant = secondchercalibconstant/sum(weights)
	#secondscincalibconstant = secondscincalibconstant/sum(weights)
	secondchercalibconstant = secondchercalibconstant/(int(noftrained)/2)
	secondscincalibconstant = secondscincalibconstant/(int(noftrained)/2)
	cluster.secondMLEnergy = (cluster.chersignal*secondchercalibconstant + cluster.scinsignal*secondscincalibconstant)/2
	cluster.secondMLScinEnergy = cluster.scinsignal*secondscincalibconstant
	cluster.secondMLCherEnergy = cluster.chersignal*secondchercalibconstant
	#return all energies
	return cluster.MLEnergy, cluster.MLScinEnergy, cluster.MLCherEnergy, cluster.secondMLEnergy, cluster.secondMLScinEnergy, cluster.secondMLCherEnergy, MLFEM

'''#not using it -> going to cancel it
def differencefile(cluster, trainedfile):
	"""Function used to find the trainedfile belonging
	to hadrons or electrons with same closest energy to custer"""
	listofcluster = unpicklelistofclusters(trainedfile)
	trainedcluster = listofcluster[0]
	difference = cluster.scinsignal - trainedcluster.scinsignal
	return abs(difference)

def find_twostepenergyfile(cluster):
	"""Find energy with trained file with closest energy"""
	listofscinsignals = []
	listoffile = glob.glob("../fasttrained-rotated/*.p")
	for trainedfile in listoffile:
		listofscinsignals.append(differencefile(cluster, trainedfile))
	listoffile.sort(key=dict(zip(listoffile, listofscinsignals)).get)
	namefile = str(listoffile[0])
	print "Doing two step analysis with: " + str(namefile)
	return namefile
'''




