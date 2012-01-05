import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("DQMServices.Core.DQMStore_cfg")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("rfio:/castor/cern.ch/user/j/jalimena/176304/HT/out_176304_HT_0.root")
                            )

process.load("HLTrigger.HLTanalyzers.hltOfflineReproducibilityDQM_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cfi')
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['hltonline']
process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')

process.DQMoutput = cms.OutputModule("PoolOutputModule",
                                     splitLevel = cms.untracked.int32(0),
                                     outputCommands = process.DQMEventContent.outputCommands,
                                     fileName = cms.untracked.string('MyFirstDQMExample.root'),
                                     dataset = cms.untracked.PSet(
    filterName = cms.untracked.string(''),
    dataTier = cms.untracked.string('')
    )
                                     )

process.p_step = cms.Path(process.hltOfflineReproducibilityDQM)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

process.schedule = cms.Schedule(process.p_step,process.endjob_step,process.DQMoutput_step)
