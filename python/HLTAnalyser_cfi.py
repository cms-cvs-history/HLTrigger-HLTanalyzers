import FWCore.ParameterSet.Config as cms

hltanalysis = cms.EDAnalyzer("HLTAnalyzer",
    mctruth = cms.InputTag("genParticleCandidates"),
    genEventScale = cms.InputTag("genEventScale"),
    genjets = cms.InputTag("iterativeCone5GenJets"),
    genmet = cms.InputTag("genMet"),
    recjets = cms.InputTag("MCJetCorJetIcone5"),
    recmet = cms.InputTag("met"),
    ht = cms.InputTag("htMet"),
    calotowers = cms.InputTag("towerMaker"),
    muon = cms.InputTag("muons"),
    Electron = cms.InputTag("pixelMatchGsfElectrons"),
    Photon = cms.InputTag("correctedPhotons"),
    l1GtObjectMapRecord = cms.InputTag("l1GtEmulDigis"),
#    l1GctCounts = cms.InputTag("l1GctEmulDigis"),
    l1GctCounts = cms.InputTag("hltGctDigis"),
    l1GtReadoutRecord = cms.InputTag("l1GmtEmulDigis"),
    l1extramc = cms.string('hltL1extraParticles'),
    hltresults = cms.InputTag("TriggerResults"),
    MuCandTag2 = cms.InputTag("hltL2MuonCandidates"),
    MuCandTag3 = cms.InputTag("hltL3MuonCandidates"),
    MuIsolTag3 = cms.InputTag("hltL3MuonIsolations"),
    MuIsolTag2 = cms.InputTag("hltL2MuonIsolations"),
    MuLinkTag = cms.InputTag("hltL3Muons"),
    RunParameters = cms.PSet(
        GenJetMin = cms.double(0.0),
        Monte = cms.bool(True),
        CalJetMin = cms.double(0.0),
        HistogramFile = cms.string('TEST.root'),
        EtaMin = cms.double(-5.2),
        Debug = cms.bool(False),
        EtaMax = cms.double(5.2)
    )
)



