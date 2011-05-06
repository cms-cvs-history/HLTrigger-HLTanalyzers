#include "HLTMessages.h"

const char * kSubEventMap                 = "subevent map";
const char * kHLTjets                     = "uncorrected HLT jets";
const char * kHLTCorjets                  = "corrected HLT jets";
const char * kRecjets                     = "uncorrected reconstructed jets";
const char * kRecCorjets                  = "corrected reconstructed jets";
const char * kGenjets                     = "generator jets";
const char * kRecmet                      = "reconstructed MET";
const char * kGenmet                      = "generator MET";
const char * kCaloTowers                  = "calo towers";
const char * kHt                          = "HT object";
const char * kRecoPFJets                  = "reco particle flow jets"; 
const char * kMuon                        = "muon candidates";
const char * kTaus                        = "tau candidates";
const char * kPFTaus                      = "pftau candidates";
const char * kPFTausTightCone             = "pftau candidates tight cone";
const char * kPFJets                      = "particle flow jets";
const char * kRecoPFTaus                            = "reco pftau candidates";
const char * ktheRecoPFTauDiscrByTanCOnePercent     = "reco PFTau Discriminator By TanC One Percent"; 
const char * ktheRecoPFTauDiscrByTanCHalfPercent    = "reco PFTau Discriminator By TanC Half Percent"; 
const char * ktheRecoPFTauDiscrByTanCQuarterPercent = "reco PFTau Discriminator By TanC Quarter Percent"; 
const char * ktheRecoPFTauDiscrByTanCTenthPercent   = "reco PFTau Discriminator By TanC Tenth Percent"; 
const char * ktheRecoPFTauDiscrByIsolation          = "reco PFTau Discriminator By Isolation";
const char * ktheRecoPFTauDiscrAgainstMuon          = "reco PFTau Discriminator Against Muon";
const char * ktheRecoPFTauDiscrAgainstElec          = "reco PFTau Discriminator Against Elec"; 
const char * kHltresults                  = "HLT results";
const char * kL1extemi                    = "L1 isolated EM objects";
const char * kL1extemn                    = "L1 non isolated EM objects";
const char * kL1extmu                     = "L1 muon objects";
const char * kL1extjetc                   = "L1 central jet objects";
const char * kL1extjetf                   = "L1 forward jet objects";
const char * kL1exttaujet                 = "L1 tau jet objects";
const char * kL1extmet                    = "L1 EtMiss object";
const char * kL1extmht                    = "L1 HtMiss object";
const char * kL1GtRR                      = "L1 GT readout record";
const char * kL1GtOMRec                   = "L1 GT object map";
const char * kL1GctBitCounts              = "L1 GCT HF bit counts";
const char * kL1GctRingSums               = "L1 GCT HF ring sums";
const char * kMctruth                     = "generator particles";
const char * kSimhit                      = "SimHit information";
const char * kGenEventInfo                = "generator information";
const char * kMucands2                    = "L2 muon candidates";
const char * kMucands3                    = "L3 muon candidates";
const char * kMunovtxcands2               = "L2 no-vertex muon candidates"; 
const char * kIsoMap2                     = "L2 muon isolation map";
const char * kIsoMap3                     = "L3 muon isolation map";
const char * kMulinks                     = "L3 muon link";
const char * kOniaPixelCands              = "Pixel track candidates in resonance with a L3 muon";
const char * kOniaTrackCands              = "Strip track candidates in resonance with a L3 muon";
const char * kDimuvtxcands3               = "L3 dimuon vertex";

const char * kBTagJets                    = "L2 b-jet collection";
const char * kBTagCorrectedJets           = "L2 calibrated b-jet collection";
const char * kBTagLifetimeBJetsL25        = "L2.5 b-jet lifetime tags";
const char * kBTagLifetimeBJetsL3         = "L3 b-jet lifetime tags";
const char * kBTagLifetimeBJetsL25SingleTrack = "L2.5 b-jet lifetime tags (SingleTrack)";
const char * kBTagLifetimeBJetsL3SingleTrack  = "L3 b-jet lifetime tags (SingleTrack)";
const char * kBTagSoftmuonBJetsL25        = "L2.5 b-jet soft muon tags";
const char * kBTagSoftmuonBJetsL3         = "L3 b-jet soft muon tags";
const char * kBTagPerformanceBJetsL25     = "L2.5 b-jet perf. meas. tag";
const char * kBTagPerformanceBJetsL3      = "L3 b-jet perf. meas. tag";

const char * kElectrons                   = "electron candidates";
const char * kPhotons                     = "photon candidates";
const char * kCandIso                     = "isol eg candidate";
const char * kCandNonIso                  = "non-isol eg candidate";
const char * kEcalIso                     = "Ecal isol map";
const char * kEcalNonIso                  = "Ecal non-isol map";
const char * kHcalIsoPho                  = "Hcal isol photon map";
const char * kHcalNonIsoPho               = "Hcal non-isol photon map";
const char * kIsoPhoTrackIsol             = "Track isol photon map";
const char * kNonIsoPhoTrackIsol          = "Track non-isol photon map";
const char * kHFECALClusters              = "HF ECAL clusters";   
const char * kHFElectrons                 = "HF Electrons"; 
const char * kIsoElectron                 = "isol electron";
const char * kNonIsoElectron              = "isol electron";
const char * kIsoEleHcal                  = "isol Hcal electron";
const char * kNonIsoEleHcal               = "isol Hcal electron"; 
const char * kIsoEleTrackIsol             = "isol Track electron";
const char * kNonIsoEleTrackIsol          = "isol Track electron";
const char * kL1IsoPixelSeeds             = "pixelSeed-SC association map for electron";
const char * kL1NonIsoPixelSeeds          = "pixelSeed-SC for electron";
const char * kNonIsoR9                    = "Spike-cleaning";
const char * kIsoR9                       = "Spike-cleaning"; 
const char * kNonIsoR9ID                  = "isol R9 ID";
const char * kIsoR9ID                     = "non-isol R9 ID";
const char * kIsoHoverEH                  = "H for H/E isol photon map";
const char * kNonIsoHoverEH               = "H for H/E non-isol photon map";

const char * kEErechits                   = "ECAL Endcap RecHits";
const char * kEBrechits                   = "ECAL Barrel RecHits"; 
const char * kHBHErechits                 = "HCAL Endcap-Barrel RecHits"; 
const char * kHOrechits                   = "HCAL HO RecHits";  
const char * kHFrechits                   = "HCAL HF RecHits"; 
const char * kpi0EErechits                = "ECAL pi0 Endcap RecHits"; 
const char * kpi0EBrechits                = "ECAL pi0 Barrel RecHits";  
const char * kIsoPixelTracksL3            = "L3 Iso Pixel Tracks"; 
const char * kIsoPixelTracksL2            = "L2 Iso Pixel Tracks";
const char * kIsoPixelTrackVertices       = "Pixel Vertices";
const char * kPixelTracksL3               = "L3 Pixel Tracks"; 
const char * kRecoVerticesHLT             = "Reconstructed vertices, HLT"; 
const char * kRecoVerticesOffline0        = "Reconstructed vertices, Offline0";
const char * kRecoVerticesOffline1        = "Reconstructed vertices, Offline1";
