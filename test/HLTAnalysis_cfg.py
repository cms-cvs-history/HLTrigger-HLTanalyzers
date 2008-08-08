import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/relval/2008/6/20/RelVal-RelValZTT-1213920853/0000/38CB1BCE-863E-DD11-A69B-001617DBD332.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)
process.hltanalysis = cms.EDAnalyzer("HLTAnalyzer",
    mctruth = cms.InputTag('genParticles'),
    genEventScale = cms.InputTag('genEventScale'),
    recjets = cms.InputTag('MCJetCorJetIcone5'),
    genjets = cms.InputTag('iterativeCone5GenJets'),
    recmet = cms.InputTag('met'),
    genmet = cms.InputTag('genMet'),
    calotowers = cms.InputTag('towerMaker'),
    muon = cms.InputTag('muons'),
    Photon = cms.InputTag('correctedPhotons'),
    Electron = cms.InputTag('pixelMatchGsfElectrons'),
    ht = cms.InputTag('htMet'),
    l1GtObjectMapRecord = cms.InputTag('hltL1GtObjectMap'),
    l1GtReadoutRecord = cms.InputTag('hltGtDigis'),
    l1GctCounts = cms.InputTag("hltGctDigis"),
    l1extramc = cms.string('hltL1extraParticles'),
    hltresults = cms.InputTag('TriggerResults'),
    RunParameters = cms.PSet(
        GenJetMin = cms.double(10.0),
        Monte = cms.bool(True),
        CalJetMin = cms.double(10.0),
        HistogramFile = cms.string('TEST_HLTAnalyzer.root'),
        EtaMin = cms.double(-5.2),
        Debug = cms.bool(False),
        EtaMax = cms.double(5.2)
    )
)

process.p1 = cms.Path(process.hltanalysis)


