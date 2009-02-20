import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_0_0_pre6/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_30X_v1/0005/0C40A6C6-D4DD-DD11-AE88-0019B9F6C674.root',
    '/store/relval/CMSSW_3_0_0_pre6/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_30X_v1/0005/1240DB12-D5DD-DD11-84FB-001D09F29619.root',
    '/store/relval/CMSSW_3_0_0_pre6/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_30X_v1/0005/32F858F2-D6DD-DD11-9A51-000423D6CA72.root'

    ),
    secondaryFileNames = cms.untracked.vstring(
    )
)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 100 )
)

process.load("Configuration.StandardSequences.GeometryPilot2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# Conditions: fake or frontier
# process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'STARTUP_30X::All'

process.load("Configuration.StandardSequences.L1Emulator_cff")
# Choose a menu/prescale/mask from one of the choices
# in L1TriggerConfig.L1GtConfigProducers.Luminosity
#process.load("Configuration.StandardSequences.L1TriggerDefaultMenu_cff")

# Run HLT
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.load("HLTrigger.Configuration.HLT_2E30_cff")
#process.schedule = process.HLTSchedule

#process.hltL1gtTrigReport = cms.EDAnalyzer( "L1GtTrigReport",
#    UseL1GlobalTriggerRecord = cms.bool( False ),
#    L1GtRecordInputTag = cms.InputTag( "hltGtDigis" )
#)
#process.hltTrigReport = cms.EDAnalyzer( "HLTrigReport",
#    HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT' )
#)
#process.HLTAnalyzerEndpath = cms.EndPath( process.hltL1gtTrigReport + process.hltTrigReport )
#process.schedule.append(process.HLTAnalyzerEndpath)

# OpenHLT specificss
# Define the HLT reco paths
process.load("HLTrigger.HLTanalyzers.HLTopen_cff")

# AlCa OpenHLT specific settings
process.hltEcalRegionalRestFEDs.Pi0ListToIgnore =  cms.InputTag("hltEcalRegionalPi0FEDs")
process.hltEcalRegionalJetsFEDs.Pi0ListToIgnore =  cms.InputTag("hltEcalRegionalPi0FEDs")
process.hltEcalRegionalEgammaFEDs.Pi0ListToIgnore =  cms.InputTag("hltEcalRegionalPi0FEDs")
process.hltEcalRegionalMuonsFEDs.Pi0ListToIgnore =  cms.InputTag("hltEcalRegionalPi0FEDs")
#process.hltEcalRegionalTausFEDs.Pi0ListToIgnore =  cms.InputTag("hltEcalRegionalPi0FEDs")

# Define the analyzer modules
process.load("HLTrigger.HLTanalyzers.HLTAnalyser_cfi")
process.analyzeThis = cms.Path( process.hltanalysis )

# Schedule the whole thing
process.schedule = cms.Schedule( 
    process.DoHltMuon, 
    process.DoHLTJets, 
    process.DoHLTPhoton, 
    process.DoHLTElectron, 
    process.DoHLTElectronStartUpWindows, 
    process.DoHLTElectronLargeWindows, 
#    process.DoHLTTau, 
    process.DoHLTBTag,
    process.DoHLTAlCaECALPhiSym,
    process.DoHLTAlCaPi0,
    process.analyzeThis )
