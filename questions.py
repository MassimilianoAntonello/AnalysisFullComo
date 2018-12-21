class Calibrationanswer:
	"""Class with all calibration answers"""
	def __init__(self, absorbermaterial, electronnamefile, electronNofEvents, pionnamefile, pionNofEvents, FolderName):
		self.absorbermaterial = absorbermaterial
		self.electronnamefile = electronnamefile
		self.pionnamefile = pionnamefile
		self.pionNofEvents = pionNofEvents
		self.FolderName = FolderName
		self.electronNofEvents = electronNofEvents

class Traininganswer:
	"""Class with all training answers"""
	def __init__(self, absorbermaterial, namefile, trainednamefile, particletype, NofEvents):
		self.absorbermaterial = absorbermaterial
		self.trainednamefile = trainednamefile
		self.NofEvents = NofEvents
		self.particletype = particletype #electron or hadron
		self.namefile = namefile

class Minimizationanswer:
	"""Class with all minimization answers"""
	def __init__(self, absorbermaterial, namefile, NofEvents, FolderName):
		self.absorbermaterial = absorbermaterial
		self.NofEvents = NofEvents
		self.FolderName = FolderName
		self.namefile = namefile

class Analysisanswer:
	"""Class with all analysis answers"""
	def __init__(self, absorbermaterial, namefile, NofEvents, Noftrainedevents, particletype, answerdualreadout, answermachinelearning, answerminimization,  foldername, answermincalib, trainednamefile):
		self.absorbermaterial = absorbermaterial
		self.namefile = namefile
		self.NofEvents = NofEvents
		self.trainednamefile = trainednamefile
		self.answerdualreadout = answerdualreadout
		self.answerminimization = answerminimization
		self.answermachinelearning = answermachinelearning
		self.answermincalib = answermincalib
		self.particletype = particletype
		self.Noftrainedevents = Noftrainedevents
		self.foldername = foldername