import FWCore.ParameterSet.Config as cms

process = cms.Process("Tester")

process.load("FWCore.MessageService.MessageLogger_cfi")


process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string('testerOutput.root')
                                   )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring('file:RootFiles/RAW2DIGI_L1Extra.root')
                            )

# RawToDigi and Necessary Configuartion Files
process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('Configuration.EventContent.EventContent_cff')
#process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
#process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')

#Global Tag
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag='PostLS162_V1::All'

process.tester = cms.EDAnalyzer(
    'tester',
    genSrc = cms.InputTag("genParticles")
    #l1GlobalTriggerReadoutSrc = cms.InputTag("gtDigis")
    #l1MuonSrc=cms.InputTag("l1extraParticles"),
    #stdMuSrc = cms.InputTag("standAloneMuons"),
    #horecoSrc = cms.InputTag("horeco"),
    #L1GtTmLInputTag = cms.InputTag("l1GtTriggerMenuLite")
    )

#Path definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.l1extra_step = cms.Path(process.L1Extra)
#process.horeco_step = cms.Path(process.horeco)
#process.tester_step = cms.Path(process.tester)

#Schedule Definition
#process.schedule = cms.Schedule(process.tester)

process.p = cms.Path(process.tester)

#process.p = cms.Path(process.RawToDigi * process.L1Extra*process.demo)
