from ROOT import TFile
from ROOT import TDirectory
from ROOT import gDirectory
from ROOT import TH1F
from ROOT import TH2F
from math import sqrt
#
# Setup
#

#Working Directory Location
workDirLoc = '/home/cranelli/HO_Muon/Fall13/RAW/CMSSW_6_2_9/src/'

#File with Weights and Data Set Prefixes
inFileLoc = workDirLoc + 'Analysis/hoTriggerAnalyzer/test/WeightingFiles/Weights.txt'
#Location of Histograms made from the Data Sets
rootFileDir = "/data/users/cranelli/HOL1Muon/HOL1Muon_Histograms/QCD/"
version = "/Version_1_3/"
rootFileName = "/L1MuonHistogram.root"
histDirName = "demo/"


#
# Begining of Script
#

outRootFileLoc = rootFileDir + version +"WeightedMerger/" + rootFileName
outRootFile = TFile(outRootFileLoc, "RECREATE")
histDir = outRootFile.mkdir(histDirName)

def weightedMerger(histogramName):

    inFile = open(inFileLoc, 'r')
    
    
    hMerge = TH1F()
    #print type(hMerge)

    #
    # Loop Over Data Sets in the Weights File 
    #

    firstDataset = True
    for line in inFile:
        entries = line.split()

        if entries[0] == "DataSet": continue #Skip Header

        dataset = entries[0]
        weight = float(entries[1])

       
        rootFileLoc = rootFileDir + version + dataset+rootFileName
        rootFile = TFile(rootFileLoc, "READ")
        # h1 = TH1F()
        h1 = rootFile.Get("demo/"+histogramName)
    
        
        if firstDataset :
            firstDataset = False
            outRootFile.cd()
            hMerge = h1.Clone(histogramName)
            hMerge.Reset("CE")


        h1.Sumw2()
        hMerge.Add(h1, weight)


        # Check Errors, Looks Correct
        #print dataset 
        #for i in range(0,h1.GetNbinsX()):
            #print h1.GetBinCenter(i), h1.GetBinContent(i)
            #print h1.GetBinError(i), sqrt(h1.GetBinContent(i))
            #print hMerge.GetBinError(i), h1.GetBinError(i)*weight
        #print type(hMerge)
  
    rootFile.Close()

    # outRootFile.Open()
    print type(hMerge)

    histDir.cd()
    hMerge.Write()
    # outRootFile.Write() 
    inFile.close()
    

# Select Which Histograms to Merge
    
weightedMerger("L1MuonBarrel_Pt")
weightedMerger("L1MuonBarrelwithMipMatch_Pt")
weightedMerger("L1MuonBarrelwithMipHLTMatch_Pt")
weightedMerger("L1MuonBarrel_Pt")
weightedMerger("L1MuonBarrelwithHLTMatch_Pt")
#weightedMerger("L1MuonBarrelandHOReco_DeltaEtaDeltaPhi")
weightedMerger("L1MuonBarrelandHLTMuon_DeltaEta")
weightedMerger("L1MuonBarrelandHLTMuon_DeltaPhi")
weightedMerger("L1MuonBarrel_Eta")
weightedMerger("L1MuonBarrel_Phi")
weightedMerger("hltMu5_Eta")
weightedMerger("hltMu5_Phi")

outRootFile.Close()

