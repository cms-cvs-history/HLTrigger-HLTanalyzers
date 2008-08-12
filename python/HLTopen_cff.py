import FWCore.ParameterSet.Config as cms
 
# import the whole HLT menu
from HLTrigger.Configuration.HLT_2E30_cff import *

# create the muon HLT reco path
DoHltMuon = cms.Path( HLTL1muonrecoSequence + HLTL2muonrecoSequence + HLTL2muonisorecoSequence + HLTL3muonrecoSequence + HLTL3muonisorecoSequence )

# create the tau HLT reco path

# ...

