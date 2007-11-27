#ifndef HLTEGAMMA_H
#define HLTEGAMMA_H

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TNamed.h"
#include <vector>
#include <map>
#include "TROOT.h"
#include "TChain.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "HLTrigger/HLTanalyzers/interface/JetUtil.h"
#include "HLTrigger/HLTanalyzers/interface/CaloTowerBoundries.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "L1TriggerConfig/L1Geometry/interface/L1CaloGeometry.h"
#include "L1TriggerConfig/L1Geometry/interface/L1CaloGeometryRecord.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SeedSuperClusterAssociation.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"
#include "L1TriggerConfig/L1Geometry/interface/L1CaloGeometry.h"
#include "L1TriggerConfig/L1Geometry/interface/L1CaloGeometryRecord.h"
#include "TTree.h"
#include "TFile.h"
#include <vector>
#include <algorithm>
#include <memory>

typedef std::vector<std::string> MyStrings;

/** \class HLTEgamma
  *  
  * $Date: November 2006
  * $Revision: 
  * \author P. Bargassa - Rice U.
  */
class HLTEgamma {
public:
  HLTEgamma(); 

  void setup(const edm::ParameterSet& pSet, TTree* tree);

  /** Analyze the Data */
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup,
	       const ElectronCollection& electron,
	       const PhotonCollection& photon,
	       const CaloGeometry& geom,
	       TTree* tree);

  void MakeL1IsolatedPhotons(edm::Event const& e, edm::EventSetup const& iSetup);
  void MakeL1NonIsolatedPhotons(edm::Event const& e, edm::EventSetup const& iSetup);
  void MakeL1IsolatedElectrons(edm::Event const& e, edm::EventSetup const& iSetup);
  void MakeL1NonIsolatedElectrons(edm::Event const& e, edm::EventSetup const& iSetup);

private:

  // Tree variables
  float *elpt, *elphi, *eleta, *elet, *ele; 
  float *photonpt, *photonphi, *photoneta, *photonet, *photone; 
  float *hphotet, *hphoteta, *hphotphi, *hphoteiso, *hphothiso, *hphottiso;
  float *heleet,*heleeta,*helephi,*heleE,*helep,*helehiso,*heletiso;
  int *hphotl1iso,*helel1iso;

  int nele,nphoton,nhltgam,nhltele;

  std::string CandIso_,CandNonIso_,EcalNonIso_,EcalIso_,HcalIsoPho_,HcalNonIsoPho_,IsoPhoTrackIsol_,NonIsoPhoTrackIsol_;
  std::string IsoEleHcalTag_,NonIsoEleHcalTag_,IsoElectronTag_,NonIsoElectronTag_,IsoEleTrackIsolTag_,NonIsoEleTrackIsolTag_;
  
  class myHLTPhoton {
  public:
    float Et;
    float eta;
    float phi;
    float ecalIsol;
    float hcalIsol;
    float trackIsol;
    bool L1Isolated;

    float et() const { return Et; } // Function defined as such to be compatible with EtGreater()
  };
  std::vector<myHLTPhoton> theHLTPhotons;

  class myHLTElectron {
  public:
    float Et;
    float eta;
    float phi;
    float E;
    float p;
    float hcalIsol;
    float trackIsol;
    bool L1Isolated;

    float et() const {return Et;}
  };
  std::vector<myHLTElectron>  theHLTElectrons;

// input variables
  bool _Monte,_Debug;
  int evtCounter;
  const float etaBarrel() {return 1.4;}

};

#endif
