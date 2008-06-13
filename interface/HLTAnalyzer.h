#include <iostream>

#include "HLTrigger/HLTanalyzers/interface/HLTEgamma.h"
#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"
#include "HLTrigger/HLTanalyzers/interface/HLTJets.h"
#include "HLTrigger/HLTanalyzers/interface/HLTBJet.h"
#include "HLTrigger/HLTanalyzers/interface/HLTMCtruth.h"
#include "HLTrigger/HLTanalyzers/interface/HLTMuon.h"
#include "HLTrigger/HLTanalyzers/interface/EventHeader.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
/* #include "DataFormats/L1GlobalTrigger/interface/L1GtLogicParser.h" */

/** \class HLTAnalyzer
  *  
  * $Date: November 2006
  * $Revision: 
  * \author P. Bargassa - Rice U.
  */

class HLTAnalyzer : public edm::EDAnalyzer {
public:
  explicit HLTAnalyzer(edm::ParameterSet const& conf);
  virtual void analyze(edm::Event const& e, edm::EventSetup const& iSetup);
  virtual void endJob();

  // Analysis tree to be filled
  TTree *HltTree;

private:
  // variables persistent across events should be declared here.
  //
  ///Default analyses
 
  EventHeader evt_header_;
  HLTInfo     hlt_analysis_;
  HLTEgamma   elm_analysis_;
  HLTJets     jet_analysis_;
  HLTBJet     bjet_analysis_;
  HLTMCtruth  mct_analysis_;
  HLTMuon     muon_analysis_;

  edm::InputTag recjets_,genjets_,recmet_,genmet_,ht_, calotowers_,hltobj_,hltresults_,genEventScale_;
  edm::InputTag Electron_,Photon_,muon_;
  std::string l1extramc_;
  edm::InputTag particleMapSource_,mctruth_; 
/*   std::string ecalDigisLabel_,hcalDigisLabel_; */

  edm::InputTag MuCandTag2_,MuIsolTag2_,MuCandTag3_,MuIsolTag3_,MuLinkTag_;
  edm::InputTag myHLT1Tau_/*,myHLT2Tau_*/;
  edm::InputTag gtReadoutRecord_,gtObjectMap_; 
  edm::InputTag gctCounts_;
  edm::InputTag lifetimeBjetL2_, lifetimeBjetL25_, lifetimeBjetL3_;
  edm::InputTag lifetimeBjetPixelTracks_, lifetimeBjetRegionalTracks_;
  edm::InputTag softmuonBjetL2_, softmuonBjetL25_, softmuonBjetL3_;
  edm::InputTag performanceBjetL2_, performanceBjetL25_, performanceBjetL3_;

  int errCnt;
  const int errMax() { return 100; }

  string _HistName; // Name of histogram file
  double _EtaMin,_EtaMax;
  TFile* m_file; // pointer to Histogram file

};
