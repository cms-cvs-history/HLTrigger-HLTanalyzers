#ifndef HLTMUON_H
#define HLTMUON_H

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TNamed.h"
#include <vector>
#include <map>
#include "TROOT.h"
#include "TChain.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/MuonReco/interface/MuIsoDeposit.h"

#include "HLTrigger/HLTanalyzers/interface/JetUtil.h"
#include "HLTrigger/HLTanalyzers/interface/CaloTowerBoundries.h"

#include "DataFormats/METReco/interface/CaloMETCollection.h"

typedef std::vector<std::string> MyStrings;

/** \class HLTMuon
  *  
  * $Date: November 2006
  * $Revision: 
  * \author P. Bargassa - Rice U.
  */
	// Revised V. Rekovic - U Minnesota, Date: December 2007
	// (added handles muonl2pterr, muonl3pterr (relative errors on track pT)
class HLTMuon {
public:
  HLTMuon(); 

  void setup(const edm::ParameterSet& pSet, TTree* tree);

  /** Analyze the Data */
  void analyze(const MuonCollection& muon,
	       const RecoChargedCandidateCollection& mucands2,
	       const reco::MuIsoAssociationMap& isoMap2,
	       const RecoChargedCandidateCollection& mucands3,
	       const reco::MuIsoAssociationMap& isoMap3,
	       const CaloGeometry& geom,
	       TTree* tree);


private:

  // Tree variables
  float *muonpt, *muonphi, *muoneta, *muonet, *muone; 
  float *muonl2pt, *muonl2eta, *muonl2phi, *muonl2dr, *muonl2dz;
  float *muonl3pt, *muonl3eta, *muonl3phi, *muonl3dr, *muonl3dz;
  float *muonl2pterr, *muonl3pterr;
  int nmuon, nmu2cand, nmu3cand;
  int *muonl2chg, *muonl2iso, *muonl3chg, *muonl3iso;
	

  // input variables
  bool _Monte,_Debug;

  int evtCounter;

  const float etaBarrel() {return 1.4;}

};

#endif
