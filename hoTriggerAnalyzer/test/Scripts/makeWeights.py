from ROOT import TFile
from ROOT import gDirectory
from ROOT import TH1F

#
# Setup
#

#Working Directory Location
workDirLoc = '/home/cranelli/HO_Muon/Fall13/RAW/CMSSW_6_2_9/src/'

#File with Cross Section and Data Set Prefixes
inFileLoc = workDirLoc + 'Analysis/hoTriggerAnalyzer/test/WeightingFiles/CrossSection.txt'
#Location of Histograms made from the Data Sets
rootFileDir = "/data/users/cranelli/HOL1Muon/HOL1Muon_Histograms/QCD/"
version = "/Version_1_3/"
rootFileName = "/L1MuonHistogram.root"

#File where Weights and Data Set Prefixes are Written To
outFileLoc = workDirLoc + 'Analysis/hoTriggerAnalyzer/test/WeightingFiles/Weights.txt'


#
# Begining of Script
#

inFile = open(inFileLoc, 'r')
outFile = open(outFileLoc, 'w')
outFile.write("DataSet".ljust(20)+ "Weight"+ "\n")
#maxlength = 8; #For formating the Output File

#
#Loop Over Data Sets in the Cross Section File 
#

for line in inFile:
    entries = line.split()
    if entries[0] == "DataSet": continue #Skip Header
    dataset = entries[0]
    crossSection = float(entries[1])

    #if len(dataset) > maxlength: maxlength = len(dataset)
    #print dataset
    #print crossSection

    #
    # Get the number of events from the root file
    # There is an event counter called eventCountHistogram
    #

    rootFileLoc = rootFileDir + version + dataset+rootFileName
    rootFile = TFile(rootFileLoc)
    eventCountHistogram = rootFile.Get("demo/Events_Count")
    numEvents = eventCountHistogram.GetBinContent(2)
    rootFile.Close()
    #print numEvents
    #directory = file.Get("demo")
    #print eventCountHistogram


    #Calculate Weight Cross Section Divided By Num Events
    weight = crossSection/numEvents


    #Print out weights and Save to File. 
    print "Weight for " + dataset + " : ", "%.9f" %weight
    outFile.write(dataset.ljust(20) + ("%.9f" %weight).rjust(15)+'\n')
                  #repr(weight).rjust(0)+'\n')

inFile.close()
outFile.close()
    
    



