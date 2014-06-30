import FWCore.ParameterSet.Config as cms

process = cms.Process("demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string('HOMuonTree.root')
                                   )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/work/c/cranelli/'
    'public/HO_Muon/Samples/'
    '0009CDC8-CCD6-E311-A9AD-90E6BAE8CC37.root'

    #'root://xrootd.unl.edu//store//mc/'
    #'Fall13dr/QCD_Pt-300to470_Tune4C_13TeV_pythia8/GEN-SIM-RAW/'
    #'castor_tsg_PU40bx25_POSTLS162_V2-v1/'
    #'00000/001C52C5-2EA4-E311-AA23-003048678F9C.root'
    #'file:/data/users/cranelli/HOL1Muon/HOL1Muon_Samples/'
    #'Fall13dr/QCD_Pt-300to470_Tune4C_13TeV_pythia8/GEN-SIM-RAW/'
    #'RAW_QCD_Pt-300to470_PU40bx25_POSTLS162_V2.root'
    )
       
)

# RawToDigi and Necessary Configuartion Files
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
#Global Tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag='PostLS170_V5::All'
#RawToDigi
process.load('Configuration.StandardSequences.RawToDigi_cff')

#turn off HO ZS
#process.hcalRawData.HO = cms.untracked.InputTag("simHcalUnsuppressedDigis", "", "")

#L1Extra
#process.load('L1Trigger.Configuration.L1Extra_cff')

#horeco
#process.load('Configuration.StandardSequences.Reconstruction_cff')


process.demo = cms.EDAnalyzer(
    'hoMuonTreeBuilder',
    genInfoSrc = cms.InputTag("generator"),
    genSrc = cms.InputTag("genParticles"),
    l1MuonSrc=cms.InputTag("l1extraParticles"),
    #stdMuSrc = cms.InputTag("standAloneMuons"),
    horecoSrc = cms.InputTag("horeco"),
    hltSumAODSrc = cms.InputTag("hltTriggerSummaryAOD")
    #L1GtTmLInputTag = cms.InputTag("l1GtTriggerMenuLite")
    )

#Path definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.l1extra_step = cms.Path(process.L1Extra)
#process.horeco_step = cms.Path(process.horeco)
process.demo_step = cms.Path(process.demo)


#For the HLT                                                                                                             
# customisation of the process.                                                                                          

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC                        
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC
#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC               
process = customizeHLTforMC(process)
# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1
#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs     
process = customisePostLS1(process)
# End of customisation functions     

#Schedule Definition

process.schedule = cms.Schedule(process.demo_step)

#process.schedule = cms.Schedule(process.raw2digi_step, process.l1extra_step,
#                                process.horeco_step, process.demo_step)

#process.p = cms.Path(process.RawToDigi)

#process.p = cms.Path(process.RawToDigi * process.L1Extra*process.demo)
