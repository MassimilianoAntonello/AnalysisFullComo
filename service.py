def compute_average(distribution):
	"""Function to compute average of a ditribution,
	numbers stored in vectors"""
	average = 0
	for entry in range(len(distribution)):
		average = average + distribution[entry]
	average = average/len(distribution)
	return average

def compute_std(distribution):
	"""Function to compute std of a distribution,
	numbers store in vectors"""
	average = compute_average(distribution)
	std = 0
	for entry in range(len(distribution)):
		std = std + (distribution[entry] - average)**2
	std = (std/(len(distribution)-1))**0.5
	return std

def compute_correlation(distribution1, distribution2):
	"""Funcition to compute correlation factor between 
	two distributions, numbers stored in vectors"""
	average1 = compute_average(distribution1)
	average2 = compute_average(distribution2)
	std1 = compute_std(distribution1)
	std2 = compute_std(distribution2)
	newdistribution = []
	for entry in range(len(distribution1)):
		newdistribution.append((distribution1[entry]-average1)*(distribution2[entry]-average2))
	newaverage = compute_average(newdistribution)
	correlation = newaverage/(std1*std2)
	return correlation
