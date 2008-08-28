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
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "HLTrigger/HLTanalyzers/interface/JetUtil.h"
#include "HLTrigger/HLTanalyzers/interface/CaloTowerBoundries.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"

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
  void analyze(const reco::PixelMatchGsfElectronCollection& electron,
	       const reco::PhotonCollection& photon,
	       TTree* tree);


private:

  // Tree variables
  int nele, nphoton;
  float *elpt, *elphi, *eleta, *elet, *ele; 
  float *photonpt, *photonphi, *photoneta, *photonet, *photone; 

  //get hold of the pixel seed - supercluster association map
  
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
    int pixelSeeds;
    float et() const {return Et;}
    bool newSC;
  };
// input variables
  bool _Monte,_Debug;
  int evtCounter;
  const float etaBarrel() {return 1.4;}

};

#endif
