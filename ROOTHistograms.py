from ROOT import gStyle, TCanvas, TH1F, TH2F, TF1, gPad, TGraph, Fit, TGraphErrors
import map
import mapgroup
from array import array
import numpy as np

def create_firstroothistograms(PrimaryParticleName, VectorSignals, VectorSignalsCher, 
	GroupedVectorSignals, GroupedVectorSignalsCher, ScinTreshold, 
	CherTreshold, ScinMaxFiber, CherMaxFiber, NumFibers, NumModules):
	"""Function to perform ROOT histograms"""
	
	#Set ROOT histograms
	TH2Signals = TH2F("ScatterplotSignals",PrimaryParticleName,111*8,1.2*0,1.2*111,111*8,1.2*0,1.2*111)
	TH2SignalsGrouped = TH2F("ScatterplotSignalsGrouped",PrimaryParticleName,111,1.2*0,1.2*111,111,1.2*0,1.2*111)
	TH2SignalsCher = TH2F("ScatterplotSignalsCher",PrimaryParticleName,111*8,0.0*111,1.2*111,111*8,0.0*111,1.2*111)
	TH2SignalsCherGrouped = TH2F("ScatterplotSignalsCherGrouped",PrimaryParticleName,111,1.2*0,1.2*111,111,1.2*0,1.2*111)
	TH1Signals = TH1F("Scintillation",PrimaryParticleName,100,0.0,ScinMaxFiber+200.0)
	TH1SignalsCher = TH1F("Cherenkov",PrimaryParticleName,100,0.0,CherMaxFiber+5)

	#Fill histograms in for loop
	for fiberindex in range(NumFibers):
		X,Y = map.mapXY(fiberindex)
		if VectorSignals[fiberindex] > ScinTreshold:
			TH2Signals.Fill(X,Y,VectorSignals[fiberindex])
			TH1Signals.Fill(VectorSignals[fiberindex])
		if VectorSignalsCher[fiberindex] > CherTreshold:
			TH2SignalsCher.Fill(X,Y,VectorSignalsCher[fiberindex])
			TH1SignalsCher.Fill(VectorSignalsCher[fiberindex])

	for moduleindex in range(NumModules):
		X,Y = mapgroup.mapgroupedXY(moduleindex)
		TH2SignalsGrouped.Fill(X,Y,GroupedVectorSignals[moduleindex])
		TH2SignalsCherGrouped.Fill(X,Y,GroupedVectorSignalsCher[moduleindex])
		
	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetPalette(1) #Root palette style
	Style.SetOptStat(0) #Do not show statistics
	TH2Signals.SetLineWidth(0) #TH2Signals #No line width
	TH2Signals.SetLineColor(2)
	#TH2Signals.SetFillColorAlpha(2, 0.)
	XAxis = TH2Signals.GetXaxis()
	XAxis.SetTitle("x (cm)")
	XAxis.CenterTitle()
	XAxis.SetTitleOffset(1.8)
	YAxis = TH2Signals.GetYaxis()
	YAxis.SetTitle("y (cm)")
	YAxis.CenterTitle()
	YAxis.SetTitleOffset(1.8)
	ZAxis = TH2Signals.GetZaxis()
	ZAxis.SetTitle("Energy (MeV)")
	ZAxis.SetTitleOffset(1.4)
	TH2Signals.Draw("LEGO2Z 0 FB")
	gPad.SaveAs("ImageScintillation.pdf")
	TH2SignalsGrouped.SetLineWidth(0) #TH2GroupedSignals #No line width
	TH2SignalsGrouped.SetLineColor(2)
	#TH2Signals.SetFillColorAlpha(2, 0.)
	XAxis = TH2SignalsGrouped.GetXaxis()
	XAxis.SetTitle("x (cm)")
	XAxis.CenterTitle()
	XAxis.SetTitleOffset(1.8)
	YAxis = TH2SignalsGrouped.GetYaxis()
	YAxis.SetTitle("y (cm)")
	YAxis.CenterTitle()
	YAxis.SetTitleOffset(1.8)
	ZAxis = TH2SignalsGrouped.GetZaxis()
	ZAxis.SetTitle("Energy (MeV)")
	ZAxis.SetTitleOffset(1.4)
	TH2SignalsGrouped.Draw("LEGO2Z 0 FB")
	gPad.SaveAs("ImageScintillationGrouped.pdf")
	TH2SignalsCherGrouped.SetLineWidth(0) #TH2GroupedCherSignals #No line width
	TH2SignalsCherGrouped.SetLineColor(4)
	#TH2Signals.SetFillColorAlpha(2, 0.)
	XAxis = TH2SignalsCherGrouped.GetXaxis()
	XAxis.SetTitle("x (cm)")
	XAxis.CenterTitle()
	XAxis.SetTitleOffset(1.8)
	YAxis = TH2SignalsCherGrouped.GetYaxis()
	YAxis.SetTitle("y (cm)")
	YAxis.CenterTitle()
	YAxis.SetTitleOffset(1.8)
	ZAxis = TH2SignalsCherGrouped.GetZaxis()
	ZAxis.SetTitle("Energy (MeV)")
	ZAxis.SetTitleOffset(1.4)
	TH2SignalsCherGrouped.Draw("LEGO2Z 0 FB")
	gPad.SaveAs("ImageCherenkovGrouped.pdf")
	TH2SignalsCher.SetLineWidth(0) #TH2SignalsCher #No line width
	TH2SignalsCher.SetLineColor(4)
	XAxis = TH2SignalsCher.GetXaxis()
	XAxis.SetTitle("x (cm)")
	XAxis.CenterTitle()
	XAxis.SetTitleOffset(1.8)
	YAxis = TH2SignalsCher.GetYaxis()
	YAxis.SetTitle("y (cm)")
	YAxis.CenterTitle()
	YAxis.SetTitleOffset(1.8)
	ZAxis = TH2SignalsCher.GetZaxis()
	ZAxis.SetTitle("Energy (MeV)")
	ZAxis.SetTitleOffset(1.4)
	TH2SignalsCher.Draw("LEGO2Z FB 0")
	gPad.SaveAs("ImageCherenkov.pdf")
	Style.SetLineWidth(1) #TH1Signals
	Style.SetOptStat(1) #Show statistics
	gPad.SetLogy()
	gPad.SetLogx()
	XAxis = TH1Signals.GetXaxis()
	XAxis.SetTitle("Energy (MeV)")
	XAxis.SetTitleOffset(1.2)
	YAxis = TH1Signals.GetYaxis()
	YAxis.SetTitle("# fibers")
	TH1Signals.Draw()
	gPad.SaveAs("EnergyFibers.pdf")
	XAxis = TH1SignalsCher.GetXaxis() #TH1SignalsCher
	XAxis.SetTitle("# Cher p.e.")
	XAxis.SetTitleOffset(1.2)
	YAxis = TH1SignalsCher.GetYaxis()
	YAxis.SetTitle("# fibers")
	TH1SignalsCher.Draw()
	gPad.SaveAs("CherpeFibers.pdf")
	gPad.Close()

def create_secondroothistograms(traditionalreconstructedenergy, machinelearningreconstructedenergy):
	"""Function to perform ROOT histograms"""
	
	#Set ROOT histograms
	TH1TraditionalEnergy = TH1F("DR Energy","",100,0.0, np.mean(traditionalreconstructedenergy)*2)
	TH1MLEnergy = TH1F("ML Energy","",100,0.0, np.mean(machinelearningreconstructedenergy)*2)

	#Set ROOT 2D histogram
	TH2FScatter = TH2F("", "", 100, 0.0, np.mean(traditionalreconstructedenergy)*2, 100, 0.0, np.mean(machinelearningreconstructedenergy)*2) 

	#Fill histograms in for loop
	for Event in range(len(traditionalreconstructedenergy)):
		TH1TraditionalEnergy.Fill(traditionalreconstructedenergy[Event])
		TH1MLEnergy.Fill(machinelearningreconstructedenergy[Event])
		TH2FScatter.Fill(traditionalreconstructedenergy[Event], machinelearningreconstructedenergy[Event])

	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetLineWidth(1) #TH1TraditionalEnergy
	Style.SetOptStat(1) #Show statistics
	XAxis = TH1TraditionalEnergy.GetXaxis()
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TH1TraditionalEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TH1TraditionalEnergy.Fit("gaus")
	myfit = TH1TraditionalEnergy.GetFunction("gaus")
	drmean = myfit.GetParameter(1)
	drmeanerror = myfit.GetParError(1)
	drsigma = myfit.GetParameter(2)
	drsigmaerror = myfit.GetParError(2)
	Style.SetOptFit()
	TH1TraditionalEnergy.Draw()
	gPad.SaveAs("DREnergy.pdf")
	gPad.Close()
	XAxis = TH1MLEnergy.GetXaxis() #TH1SignalsCher
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TH1MLEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TH1MLEnergy.Fit("gaus")
	mymlfit = TH1MLEnergy.GetFunction("gaus")
	mlmean = mymlfit.GetParameter(1)
	mlmeanerror = mymlfit.GetParError(1)
	mlsigma = mymlfit.GetParameter(2)
	mlsigmaerror = mymlfit.GetParError(2)
	Style.SetOptFit()
	TH1MLEnergy.Draw()
	gPad.SaveAs("MLEnergy.pdf")	
	gPad.Close()
	Style.SetOptStat(0) #Dont shown statistics
	XAxis = TH2FScatter.GetXaxis() #TH2FScatter
	XAxis.SetTitle("DR Energy (MeV)")
	YAxis = TH2FScatter.GetYaxis()
	YAxis.SetTitle("ML Energy (MeV)")
	TH2FScatter.Draw("COLZ")
	gPad.SaveAs("ScatterEnergies.pdf")
	gPad.Close()
	return drmean, drmeanerror, drsigma, drsigmaerror, mlmean, mlmeanerror, mlsigma, mlsigmaerror 

def create_emscatterplot(trueem, differenceem):
	"""Function to perform ROOT histograms"""
	
	#Set ROOT 2D histogram
	TH2FScatter = TH2F("", "Fem", 100, 0.0, np.mean(trueem)*2, 100, 0.0, np.mean(differenceem)*2) 

	#Fill histograms in for loop
	for Event in range(len(trueem)):
		TH2FScatter.Fill(trueem[Event], differenceem[Event])

	#Draw + DrawOptions histograms		
	Style = gStyle	
	Style.SetOptStat(0) 
	XAxis = TH2FScatter.GetXaxis() #TH2FScatter
	XAxis.SetTitle("True fem")
	YAxis = TH2FScatter.GetYaxis()
	YAxis.SetTitle("Difference em (DR)")
	TH2FScatter.Draw("COLZ")
	gPad.SaveAs("ScatterfemDR.pdf")
	gPad.Close()

def create_mlemscatterplot(trueem, differencemlem):
	"""Function to perform ROOT histograms"""
	
	#Set ROOT 2D histogram
	TH2FScatter = TH2F("", "Fem", 100, 0.0, np.mean(trueem)*2, 100, 0.0, np.mean(differencemlem)*2) 

	#Fill histograms in for loop
	for Event in range(len(trueem)):
		TH2FScatter.Fill(trueem[Event], differencemlem[Event])

	#Draw + DrawOptions histograms		
	Style = gStyle	
	Style.SetOptStat(0) 
	XAxis = TH2FScatter.GetXaxis() #TH2FScatter
	XAxis.SetTitle("True fem")
	YAxis = TH2FScatter.GetYaxis()
	YAxis.SetTitle("Difference em (ML)")
	TH2FScatter.Draw("COLZ")
	gPad.SaveAs("ScatterfemML.pdf")
	gPad.Close()

def create_thirdroothistograms(traditionalscinreconstructedenergy, machinelearningscinreconstructedenergy):
	"""Function to perform ROOT histograms"""
	
	#Set ROOT histograms
	TH1TraditionalScinEnergy = TH1F("Scin Energy","",100,0.0, np.mean(traditionalscinreconstructedenergy)*2)
	TH1MLScinEnergy = TH1F("ML Scin Energy","",100,0.0, np.mean(machinelearningscinreconstructedenergy)*2)

	#Fill histograms in for loop
	for Event in range(len(traditionalscinreconstructedenergy)):
		TH1TraditionalScinEnergy.Fill(traditionalscinreconstructedenergy[Event])
		TH1MLScinEnergy.Fill(machinelearningscinreconstructedenergy[Event])

	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetLineWidth(1) #TH1TraditionalEnergy
	Style.SetOptStat(1) #Show statistics
	XAxis = TH1TraditionalScinEnergy.GetXaxis()
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TH1TraditionalScinEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TH1TraditionalScinEnergy.Draw()
	gPad.SaveAs("TraditionalScinEnergy.pdf")
	gPad.Close()	
	XAxis = TH1MLScinEnergy.GetXaxis() #TH1SignalsCher
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TH1MLScinEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TH1MLScinEnergy.Fit("gaus")
	Style.SetOptFit()
	TH1MLScinEnergy.Draw()
	gPad.SaveAs("MLScinEnergy.pdf")	
	gPad.Close()

def create_fourthroothistograms(traditionalcherreconstructedenergy, machinelearningcherreconstructedenergy):
	"""Function to perform ROOT histograms"""
	
	#Set ROOT histograms
	TH1TraditionalCherEnergy = TH1F("Cher Energy","",100,0.0, np.mean(traditionalcherreconstructedenergy)*2)
	TH1MLCherEnergy = TH1F("ML Cher Energy","",100,0.0, np.mean(machinelearningcherreconstructedenergy)*2)

	#Fill histograms in for loop
	for Event in range(len(traditionalcherreconstructedenergy)):
		TH1TraditionalCherEnergy.Fill(traditionalcherreconstructedenergy[Event])
		TH1MLCherEnergy.Fill(machinelearningcherreconstructedenergy[Event])

	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetLineWidth(1) #TH1TraditionalCherEnergy
	Style.SetOptStat(1) #Show statistics
	XAxis = TH1TraditionalCherEnergy.GetXaxis()
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TH1TraditionalCherEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TH1TraditionalCherEnergy.Draw()
	gPad.SaveAs("TraditionalCherEnergy.pdf")
	gPad.Close()
	XAxis = TH1MLCherEnergy.GetXaxis() #TH1SignalsCher
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TH1MLCherEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	Style.SetOptFit()
	TH1MLCherEnergy.Fit("gaus")
	TH1MLCherEnergy.Draw()
	gPad.SaveAs("MLCherEnergy.pdf")	
	gPad.Close()

def create_fifthroothistograms(CombinedEnergy):
	"""Function to perform ROOT histograms"""
	
	#Set ROOT histograms
	TH1CombinedEnergy = TH1F("Comb Energy","",100,0.0, np.mean(CombinedEnergy)*2)

	#Fill histograms in for loop
	for Event in range(len(CombinedEnergy)):
		TH1CombinedEnergy.Fill(CombinedEnergy[Event])

	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetLineWidth(1) #TH1TraditionalCherEnergy
	Style.SetOptStat(1) #Show statistics
	XAxis = TH1CombinedEnergy.GetXaxis()
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TH1CombinedEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TH1CombinedEnergy.Draw()
	gPad.SaveAs("CombinedEnergy.pdf")
	gPad.Close()

def create_sixthroothistograms(SecondCombinedEnergy):
	"""Function to perform ROOT histograms"""
	
	#Set ROOT histograms
	TH1SecondCombinedEnergy = TH1F("Second Comb Energy","",100,0.0, np.mean(SecondCombinedEnergy)*2)

	#Fill histograms in for loop
	for Event in range(len(SecondCombinedEnergy)):
		TH1SecondCombinedEnergy.Fill(SecondCombinedEnergy[Event])

	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetLineWidth(1) #TH1TraditionalCherEnergy
	Style.SetOptStat(1) #Show statistics
	XAxis = TH1SecondCombinedEnergy.GetXaxis()
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TH1SecondCombinedEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TH1SecondCombinedEnergy.Draw()
	gPad.SaveAs("SecondCombinedEnergy.pdf")
	gPad.Close()

def create_seventhroothistograms(secondmlreconstructedenergy, secondmlscinreconstructedenergy, secondmlcherreconstructedenergy):
	"""Function to perform ROOT histograms"""
	
	#Set ROOT histograms
	TH1SMLEnergy = TH1F("Sec ML Energy","",100,0.0, np.mean(secondmlreconstructedenergy)*2)
	TH1SMLScinEnergy = TH1F("Sec ML Scin Energy","",100,0.0, np.mean(secondmlscinreconstructedenergy)*2)
	TH1SMLCherEnergy = TH1F("Sec ML Cher Energy","",100,0.0, np.mean(secondmlcherreconstructedenergy)*2)

	#Fill histograms in for loop
	for Event in range(len(secondmlreconstructedenergy)):
		TH1SMLEnergy.Fill(secondmlreconstructedenergy[Event])
		TH1SMLScinEnergy.Fill(secondmlscinreconstructedenergy[Event])
		TH1SMLCherEnergy.Fill(secondmlcherreconstructedenergy[Event])

	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetLineWidth(1) #TH1SMLEnergy
	Style.SetOptStat(1) #Show statistics
	XAxis = TH1SMLEnergy.GetXaxis()
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TH1SMLEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TH1SMLEnergy.Fit("gaus")
	Style.SetOptFit()
	TH1SMLEnergy.Draw()
	gPad.SaveAs("SecondMLEnergy.pdf")
	gPad.Close()
	XAxis = TH1SMLScinEnergy.GetXaxis() #TH1SMLScinEnergy
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TH1SMLScinEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TH1SMLScinEnergy.Fit("gaus")
	Style.SetOptFit()
	TH1SMLScinEnergy.Draw()
	gPad.SaveAs("SecondMLScinEnergy.pdf")
	gPad.Close()	
	XAxis = TH1SMLCherEnergy.GetXaxis() #TH1SMLCherEnergy
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TH1SMLCherEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TH1SMLCherEnergy.Fit("gaus")
	Style.SetOptFit()
	TH1SMLCherEnergy.Draw()
	gPad.SaveAs("SecondMLCherEnergy.pdf")
	gPad.Close()
''' #not using it -> going to cancel it
def create_Tsecondroothistograms(Tsecondmlreconstructedenergy, Tsecondmlscinreconstructedenergy, Tsecondmlcherreconstructedenergy):
	"""Function to perform ROOT histograms"""
	
	#Set ROOT histograms
	TTH1SMLEnergy = TH1F("TSec ML Energy","",100,0.0, max(Tsecondmlreconstructedenergy)+5000)
	TTH1SMLScinEnergy = TH1F("TSec ML Scin Energy","",100,0.0, max(Tsecondmlscinreconstructedenergy)+5000)
	TTH1SMLCherEnergy = TH1F("TSec ML Cher Energy","",100,0.0,max(Tsecondmlcherreconstructedenergy)+5000)

	#Fill histograms in for loop
	for Event in range(len(Tsecondmlreconstructedenergy)):
		TTH1SMLEnergy.Fill(Tsecondmlreconstructedenergy[Event])
		TTH1SMLScinEnergy.Fill(Tsecondmlscinreconstructedenergy[Event])
		TTH1SMLCherEnergy.Fill(Tsecondmlcherreconstructedenergy[Event])

	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetLineWidth(1) #TH1SMLEnergy
	Style.SetOptStat(1) #Show statistics
	XAxis = TTH1SMLEnergy.GetXaxis()
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TTH1SMLEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TTH1SMLEnergy.Fit("gaus")
	Style.SetOptFit()
	TTH1SMLEnergy.Draw()
	gPad.SaveAs("TSecondMLEnergy.pdf")
	XAxis = TTH1SMLScinEnergy.GetXaxis() #TH1SMLScinEnergy
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TTH1SMLScinEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TTH1SMLScinEnergy.Fit("gaus")
	Style.SetOptFit()
	TTH1SMLScinEnergy.Draw()
	gPad.SaveAs("TSecondMLScinEnergy.pdf")	
	XAxis = TTH1SMLCherEnergy.GetXaxis() #TH1SMLCherEnergy
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TTH1SMLCherEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TTH1SMLCherEnergy.Fit("gaus")
	Style.SetOptFit()
	TTH1SMLCherEnergy.Draw()
	gPad.SaveAs("TSecondMLCherEnergy.pdf")

def create_Troothistograms(Tsecondmlreconstructedenergy, Tsecondmlscinreconstructedenergy, Tsecondmlcherreconstructedenergy):
	"""Function to perform ROOT histograms"""
	
	#Set ROOT histograms
	TTH1SMLEnergy = TH1F("TSec ML Energy","",100,0.0, max(Tsecondmlreconstructedenergy)+5000)
	TTH1SMLScinEnergy = TH1F("TSec ML Scin Energy","",100,0.0, max(Tsecondmlscinreconstructedenergy)+5000)
	TTH1SMLCherEnergy = TH1F("TSec ML Cher Energy","",100,0.0,max(Tsecondmlcherreconstructedenergy)+5000)

	#Fill histograms in for loop
	for Event in range(len(Tsecondmlreconstructedenergy)):
		TTH1SMLEnergy.Fill(Tsecondmlreconstructedenergy[Event])
		TTH1SMLScinEnergy.Fill(Tsecondmlscinreconstructedenergy[Event])
		TTH1SMLCherEnergy.Fill(Tsecondmlcherreconstructedenergy[Event])

	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetLineWidth(1) #TH1SMLEnergy
	Style.SetOptStat(1) #Show statistics
	XAxis = TTH1SMLEnergy.GetXaxis()
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TTH1SMLEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TTH1SMLEnergy.Fit("gaus")
	Style.SetOptFit()
	TTH1SMLEnergy.Draw()
	gPad.SaveAs("TMLEnergy.pdf")
	XAxis = TTH1SMLScinEnergy.GetXaxis() #TH1SMLScinEnergy
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TTH1SMLScinEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TTH1SMLScinEnergy.Fit("gaus")
	Style.SetOptFit()
	TTH1SMLScinEnergy.Draw()
	gPad.SaveAs("TMLScinEnergy.pdf")	
	XAxis = TTH1SMLCherEnergy.GetXaxis() #TH1SMLCherEnergy
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TTH1SMLCherEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TTH1SMLCherEnergy.Fit("gaus")
	Style.SetOptFit()
	TTH1SMLCherEnergy.Draw()
	gPad.SaveAs("TMLCherEnergy.pdf")

def create_T2roothistograms(Tsecondmlreconstructedenergy, Tsecondmlscinreconstructedenergy, Tsecondmlcherreconstructedenergy):
	"""Function to perform ROOT histograms"""
	
	#Set ROOT histograms
	TTH1SMLEnergy = TH1F("TSec ML Energy","",100,0.0, max(Tsecondmlreconstructedenergy)+5000)
	TTH1SMLScinEnergy = TH1F("TSec ML Scin Energy","",100,0.0, max(Tsecondmlscinreconstructedenergy)+5000)
	TTH1SMLCherEnergy = TH1F("TSec ML Cher Energy","",100,0.0,max(Tsecondmlcherreconstructedenergy)+5000)

	#Fill histograms in for loop
	for Event in range(len(Tsecondmlreconstructedenergy)):
		TTH1SMLEnergy.Fill(Tsecondmlreconstructedenergy[Event])
		TTH1SMLScinEnergy.Fill(Tsecondmlscinreconstructedenergy[Event])
		TTH1SMLCherEnergy.Fill(Tsecondmlcherreconstructedenergy[Event])

	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetLineWidth(1) #TH1SMLEnergy
	Style.SetOptStat(1) #Show statistics
	XAxis = TTH1SMLEnergy.GetXaxis()
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TTH1SMLEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TTH1SMLEnergy.Fit("gaus")
	Style.SetOptFit()
	TTH1SMLEnergy.Draw()
	gPad.SaveAs("T2SecondMLEnergy.pdf")
	XAxis = TTH1SMLScinEnergy.GetXaxis() #TH1SMLScinEnergy
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TTH1SMLScinEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TTH1SMLScinEnergy.Fit("gaus")
	Style.SetOptFit()
	TTH1SMLScinEnergy.Draw()
	gPad.SaveAs("T2SecondMLScinEnergy.pdf")	
	XAxis = TTH1SMLCherEnergy.GetXaxis() #TH1SMLCherEnergy
	XAxis.SetTitle("Energy (MeV)")
	YAxis = TTH1SMLCherEnergy.GetYaxis()
	YAxis.SetTitle("# events")
	TTH1SMLCherEnergy.Fit("gaus")
	Style.SetOptFit()
	TTH1SMLCherEnergy.Draw()
	gPad.SaveAs("T2SecondMLCherEnergy.pdf")
'''
def create_eigthroothistogram(machinelearningscinreconstructedenergy, machinelearningcherreconstructedenergy):
	"""Function to perform ROOT histograms"""

	#Set ROOT 2D histogram
	TH2FMLScinCherEnergy = TH2F("", "", 100, 0.0, np.mean(machinelearningscinreconstructedenergy)*2, 100, 0.0, np.mean(machinelearningcherreconstructedenergy)*2) 

	#Fill histogram in for loop
	for event in range(len(machinelearningscinreconstructedenergy)):
		TH2FMLScinCherEnergy.Fill(machinelearningscinreconstructedenergy[event], machinelearningcherreconstructedenergy[event])

	#Draw + Draw Options
	XAxis = TH2FMLScinCherEnergy.GetXaxis() #TH2FMLScinCherEnergy
	XAxis.SetTitle("ML Scin Energy (MeV)")
	YAxis = TH2FMLScinCherEnergy.GetYaxis()
	YAxis.SetTitle("ML Cher Energy (MeV)")
	TH2FMLScinCherEnergy.Draw("COLZ")
	gPad.SaveAs("ScatterMLScinCherEnergies.pdf")
	gPad.Close()

def create_ninethroothistogram(traditionalscinreconstructedenergy, traditionalcherreconstructedenergy):
	"""Function to perform ROOT histograms"""

	#Set ROOT 2D histogram
	TH2FScinCherEnergy = TH2F("", "", 100, 0.0, np.mean(traditionalscinreconstructedenergy)*2, 100, 0.0, np.mean(traditionalcherreconstructedenergy)*2) 

	#Fill histogram in for loop
	for event in range(len(traditionalscinreconstructedenergy)):
		TH2FScinCherEnergy.Fill(traditionalscinreconstructedenergy[event], traditionalcherreconstructedenergy[event])

	#Draw + Draw Options
	XAxis = TH2FScinCherEnergy.GetXaxis() #TH2FScinCherEnergy
	XAxis.SetTitle("Scin Energy (MeV)")
	YAxis = TH2FScinCherEnergy.GetYaxis()
	YAxis.SetTitle("Cher Energy (MeV)")
	TH2FScinCherEnergy.Draw("COLZ")
	gPad.SaveAs("ScatterScinCherEnergies.pdf")
	gPad.Close()

def create_graph(Fastscinenergy, Fastcherenergy, Scinenergy, Cherenergy, fem, PrimaryParticleEnergy):
	"""Function to perform ROOT graphs and compute he values"""

	#How many graph points
	n = len(Fastscinenergy)

	TGraphfasthescin = TGraph(n, fem, Fastscinenergy)
	TGraphfasthecher = TGraph(n, fem, Fastcherenergy)
	TGraphhescin = TGraph(n, fem, Scinenergy)
	TGraphhecher = TGraph(n, fem, Cherenergy)

	#Draw + DrawOptions, Fit + parameter estimation
	Style = gStyle
	XAxis = TGraphfasthescin.GetXaxis() #TGraphfasthescin
	TGraphfasthescin.SetMarkerColor(4)
	TGraphfasthescin.SetMarkerStyle(20)
	TGraphfasthescin.SetMarkerSize(1)
	XAxis.SetTitle("fem")
	YAxis = TGraphfasthescin.GetYaxis()
	YAxis.SetTitle("Energy scin (MeV)")
	TGraphfasthescin.Fit("pol1")
	myfit = TGraphfasthescin.GetFunction("pol1")
	Fasthescin = myfit.GetParameter(0)/PrimaryParticleEnergy             
	TGraphfasthescin.Draw("AP")
	gPad.SaveAs("Fasthescin.pdf")
	gPad.Close()
	XAxis = TGraphfasthecher.GetXaxis() #TGraphfasthecher
	TGraphfasthecher.SetMarkerColor(4)
	TGraphfasthecher.SetMarkerStyle(20)
	TGraphfasthecher.SetMarkerSize(1)
	XAxis.SetTitle("fem")
	YAxis = TGraphfasthecher.GetYaxis()
	YAxis.SetTitle("Energy Cher (MeV)")
	TGraphfasthecher.Fit("pol1")
	myfit = TGraphfasthecher.GetFunction("pol1")
	Fasthecher = myfit.GetParameter(0)/PrimaryParticleEnergy  
	TGraphfasthecher.Draw("AP")
	gPad.SaveAs("Fasthecher.pdf")
	gPad.Close()
	XAxis = TGraphhescin.GetXaxis() #TGraphhescin
	TGraphhescin.SetMarkerColor(4)
	TGraphhescin.SetMarkerStyle(20)
	TGraphhescin.SetMarkerSize(1)
	XAxis.SetTitle("fem")
	YAxis = TGraphhescin.GetYaxis()
	YAxis.SetTitle("Energy scin (MeV)")
	TGraphhescin.Fit("pol1")
	myfit = TGraphhescin.GetFunction("pol1")
	hescin = myfit.GetParameter(0)/PrimaryParticleEnergy  
	TGraphhescin.Draw("AP")
	gPad.SaveAs("hescin.pdf")
	gPad.Close()
	XAxis = TGraphhecher.GetXaxis() #TGraphhecher
	TGraphhecher.SetMarkerColor(4)
	TGraphhecher.SetMarkerStyle(20)
	TGraphhecher.SetMarkerSize(1)
	XAxis.SetTitle("fem")
	YAxis = TGraphhecher.GetYaxis()
	YAxis.SetTitle("Energy cher (MeV)")
	TGraphhecher.Fit("pol1")
	myfit = TGraphhecher.GetFunction("pol1")
	hecher = myfit.GetParameter(0)/PrimaryParticleEnergy  
	TGraphhecher.Draw("AP")
	gPad.SaveAs("hecher.pdf")
	gPad.Close()
	return Fasthescin, Fasthecher, hescin, hecher

def create_fastgraph(Fastscinenergy, Fastcherenergy, fem, PrimaryParticleEnergy):
	"""Function to perform ROOT graphs and compute he values"""

	#How many graph points
	n = len(Fastscinenergy)

	TGraphfasthescin = TGraph(n, fem, Fastscinenergy)
	TGraphfasthecher = TGraph(n, fem, Fastcherenergy)
	
	#Draw + DrawOptions, Fit + parameter estimation
	Style = gStyle
	XAxis = TGraphfasthescin.GetXaxis() #TGraphfasthescin
	TGraphfasthescin.SetMarkerColor(4)
	TGraphfasthescin.SetMarkerStyle(1)
	TGraphfasthescin.SetMarkerSize(1)
	XAxis.SetTitle("fem")
	YAxis = TGraphfasthescin.GetYaxis()
	YAxis.SetTitle("Energy scin (MeV)")
	TGraphfasthescin.Fit("pol1")
	myfit = TGraphfasthescin.GetFunction("pol1")
	Fasthescin = myfit.GetParameter(0)/PrimaryParticleEnergy             
	TGraphfasthescin.Draw("AP")
	gPad.SaveAs("Fasthescin.pdf")
	gPad.Close()
	XAxis = TGraphfasthecher.GetXaxis() #TGraphfasthecher
	TGraphfasthecher.SetMarkerColor(4)
	TGraphfasthecher.SetMarkerStyle(1)
	TGraphfasthecher.SetMarkerSize(1)
	XAxis.SetTitle("fem")
	YAxis = TGraphfasthecher.GetYaxis()
	YAxis.SetTitle("Energy Cher (MeV)")
	TGraphfasthecher.Fit("pol1")
	myfit = TGraphfasthecher.GetFunction("pol1")
	Fasthecher = myfit.GetParameter(0)/PrimaryParticleEnergy  
	TGraphfasthecher.Draw("AP")
	gPad.SaveAs("Fasthecher.pdf")
	gPad.Close()
	return Fasthescin, Fasthecher

def create_resolutiongraph(n, energies, sigmasmeans, energieserrors, sigmasmeanserrors, graphname):
	"""Function to perform ROOT graphs of resolutions"""
	#How many points
	n = int(n)

	TGraphresolution = TGraphErrors(n, energies, sigmasmeans, energieserrors, sigmasmeanserrors)
	
	#Draw + DrawOptions, Fit + parameter estimation
	Style = gStyle
	Style.SetOptFit()
	XAxis = TGraphresolution.GetXaxis() #TGraphresolution
	TGraphresolution.SetMarkerColor(4)
	TGraphresolution.SetMarkerStyle(20)
	TGraphresolution.SetMarkerSize(2)
	XAxis.SetTitle("Energy (GeV)")
	YAxis = TGraphresolution.GetYaxis()
	YAxis.SetTitle("Sigma/Mean")
	resolutionfit = TF1("resolutionfit", '([0]/((x)**0.5))+[1]', 0, max(energies)) #somma non quadratura
	TGraphresolution.Fit("resolutionfit")
	a = resolutionfit.GetParameter(0)
	b = resolutionfit.GetParameter(1)             
	TGraphresolution.Draw("AP")
	gPad.SaveAs(graphname)
	gPad.Close()
	return a, b

def create_linearitygraph(n, energies, energieserrors, means, sigmameans, graphname):
	"""Function to perform ROOT graphs of resolutions"""
	#How many points
	n = int(n)

	TGraphlinearity = TGraphErrors(n, energies, means, energieserrors, sigmameans)
	
	#Draw + DrawOptions, Fit + parameter estimation
	Style = gStyle
	Style.SetOptFit()
	XAxis = TGraphlinearity.GetXaxis() #TGraphresolution
	TGraphlinearity.SetMarkerColor(4)
	TGraphlinearity.SetMarkerStyle(20)
	TGraphlinearity.SetMarkerSize(2)
	XAxis.SetTitle("Energy (GeV)")
	YAxis = TGraphlinearity.GetYaxis()
	YAxis.SetTitle("Mean/TrueEnergy")
	TGraphlinearity.Draw("AP")
	gPad.SaveAs(graphname)
	gPad.Close()

def create_fastgrapherror(Fastscinenergy, Fastcherenergy, rmsScin, rmsCher, fem, rmsfem, PrimaryParticleEnergy):
	"""Function to perform ROOT graphs and compute he values"""

	#How many graph points
	n = len(Fastscinenergy)

	TGraphfasthescin = TGraphErrors(n, fem, Fastscinenergy, rmsfem, rmsScin)
	TGraphfasthecher = TGraphErrors(n, fem, Fastcherenergy, rmsfem, rmsCher)
	
	#Draw + DrawOptions, Fit + parameter estimation
	Style = gStyle
	Style.SetOptFit(1)
	Style.SetOptStat(0) #Do not show statistics
	XAxis = TGraphfasthescin.GetXaxis() #TGraphfasthescin
	TGraphfasthescin.SetMarkerColor(4)
	TGraphfasthescin.SetMarkerStyle(1)
	TGraphfasthescin.SetMarkerSize(1)
	XAxis.SetTitle("fem")
	XAxis.SetLimits(0.0,1.0)
	YAxis = TGraphfasthescin.GetYaxis()
	YAxis.SetTitle("Energy scin (MeV)")
	YAxis.SetRange(int(0.0),int(float(PrimaryParticleEnergy)+10000.0))
	TGraphfasthescin.Fit("pol1")
	myfit = TGraphfasthescin.GetFunction("pol1")
	Fasthescin = myfit.GetParameter(0)/PrimaryParticleEnergy    
	Fasthescinerror = myfit.GetParError(0)/PrimaryParticleEnergy 
	Fasthescin2 = (PrimaryParticleEnergy - myfit.GetParameter(1))/PrimaryParticleEnergy
	TGraphfasthescin.SetTitle("")
	TGraphfasthescin.Draw("AP")
	gPad.SaveAs("Fasthescin-pol1fit.pdf")
	gPad.Close()
	combinedfit = TF1("combinedfit", '(x)*([1]-[0]*[1])+[0]*[1]', 0, max(fem))
	combinedfit.FixParameter(1, PrimaryParticleEnergy)
	TGraphfasthescin.Fit("combinedfit")
	hes = combinedfit.GetParameter(0)   
	heserror = combinedfit.GetParError(0)          
	TGraphfasthescin.Draw("AP")
	gPad.SaveAs("Fasthescin.pdf")
	gPad.Close()
	XAxis = TGraphfasthecher.GetXaxis() #TGraphfasthecher
	TGraphfasthecher.SetMarkerColor(4)
	TGraphfasthecher.SetMarkerStyle(1)
	TGraphfasthecher.SetMarkerSize(1)
	XAxis.SetTitle("fem")
	XAxis.SetLimits(0.0,1.0)
	YAxis = TGraphfasthecher.GetYaxis()
	YAxis.SetTitle("Energy Cher (MeV)")
	YAxis.SetRange(int(0.0),int(float(PrimaryParticleEnergy)+10000.0))
	TGraphfasthecher.Fit("pol1")
	myfit = TGraphfasthecher.GetFunction("pol1")
	Fasthecher = myfit.GetParameter(0)/PrimaryParticleEnergy
	Fasthechererror = myfit.GetParError(0)/PrimaryParticleEnergy
	Fasthecher2 = (PrimaryParticleEnergy - myfit.GetParameter(1))/PrimaryParticleEnergy
	TGraphfasthecher.SetTitle("")
	TGraphfasthecher.Draw("AP")
	gPad.SaveAs("Fasthecher-pol1fit.pdf")
	gPad.Close()  
	TGraphfasthecher.Fit("combinedfit")
	hec = combinedfit.GetParameter(0)  
	hecerror = combinedfit.GetParError(0)           
	TGraphfasthecher.Draw("AP")
	gPad.SaveAs("Fasthecher.pdf")
	gPad.Close()
	return Fasthescin, Fasthecher, Fasthescinerror, Fasthechererror, Fasthescin2, Fasthecher2, hes, hec, heserror, hecerror

def create_fastgraph2(xs, ys):
	"""Function to perform ROOT graph"""

	#How many graph points
	n = len(xs)

	TGraph2 = TGraph(n, xs, ys)
	
	#Draw + DrawOptions, Fit + parameter estimation
	Style = gStyle
	XAxis = TGraph2.GetXaxis() #TGraphfasthescin
	TGraph2.SetMarkerColor(4)
	TGraph2.SetMarkerStyle(20)
	TGraph2.SetMarkerSize(1)
	XAxis.SetLimits(0.,0.6)
	XAxis.SetTitle("")
	YAxis = TGraph2.GetYaxis()
	YAxis.SetTitle("")
	TGraph2.Draw("ACP")
	gPad.SaveAs("Fastgraph.pdf")
	gPad.Close()

def create_fastroothistograms(vector, histogramtitle, xtitle, ytitle, histogramname):
	"""Function to perform ROOT histograms"""

	rms = 0
	
	#Set ROOT histograms
	TH1Hist = TH1F(histogramtitle,"",100,0.0, np.mean(vector)*2)

	#Fill histograms in for loop
	for entry in range(len(vector)):
		TH1Hist.Fill(vector[entry])

	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetLineWidth(1) #TH1Hist
	Style.SetOptStat(1) #Show statistics
	XAxis = TH1Hist.GetXaxis()
	XAxis.SetTitle(xtitle)
	YAxis = TH1Hist.GetYaxis()
	YAxis.SetTitle(ytitle)
	TH1Hist.Draw()
	mean = TH1Hist.GetMean()
	rms = TH1Hist.GetRMS()
	Entries = TH1Hist.GetEntries()
	gPad.SaveAs(histogramname)
	gPad.Close()
	return mean, rms, Entries

def create_fastroothistogramsgaus(vector, histogramtitle, xtitle, ytitle, histogramname):
	"""Function to perform ROOT histograms"""

	rms = 0
	
	#Set ROOT histograms
	TH1Hist = TH1F(histogramtitle,"",100,0.0, np.mean(vector)*2)

	#Fill histograms in for loop
	for entry in range(len(vector)):
		TH1Hist.Fill(vector[entry])

	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetLineWidth(1) #TH1Hist
	Style.SetOptStat(1) #Show statistics
	XAxis = TH1Hist.GetXaxis()
	XAxis.SetTitle(xtitle)
	YAxis = TH1Hist.GetYaxis()
	YAxis.SetTitle(ytitle)
	TH1Hist.Fit("gaus")
	Style.SetOptFit()
	TH1Hist.Draw()
	myfit = TH1Hist.GetFunction("gaus")
	mean = myfit.GetParameter(1)
	rms = myfit.GetParameter(2)
	Entries = TH1Hist.GetEntries()
	gPad.SaveAs(histogramname)
	gPad.Close()
	return mean, rms, Entries


def create_minenergyroothistograms(vector, histogramtitle, xtitle, ytitle, histogramname):
	"""Function to perform ROOT histograms"""

	rms = 0
	
	#Set ROOT histograms
	TH1Hist = TH1F(histogramtitle,"",100,0.0, np.mean(vector)*2)

	#Fill histograms in for loop
	for entry in range(len(vector)):
		TH1Hist.Fill(vector[entry])

	#Draw + DrawOptions histograms		
	Style = gStyle
	Style.SetLineWidth(1) #TH1Hist
	Style.SetOptStat(1) #Show statistics
	XAxis = TH1Hist.GetXaxis()
	XAxis.SetTitle(xtitle)
	YAxis = TH1Hist.GetYaxis()
	YAxis.SetTitle(ytitle)
	TH1Hist.Fit("gaus")
	myminfit = TH1Hist.GetFunction("gaus")
	minmean = myminfit.GetParameter(1)
	minmeanerror = myminfit.GetParError(1)
	minsigma = myminfit.GetParameter(2)
	minsigmaerror = myminfit.GetParError(2)
	Style.SetOptFit()
	TH1Hist.Draw()
	mean = TH1Hist.GetMean()
	rms = TH1Hist.GetRMS()
	Entries = TH1Hist.GetEntries()
	gPad.SaveAs(histogramname)
	gPad.Close()
	return minmean, minmeanerror, minsigma, minsigmaerror



	

