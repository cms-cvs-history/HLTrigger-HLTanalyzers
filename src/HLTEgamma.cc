#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <cstdlib>
#include <cstring>

#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeedFwd.h"
#include "HLTrigger/HLTanalyzers/interface/HLTEgamma.h"

HLTEgamma::HLTEgamma() {
  evtCounter = 0;
}

/*  Setup the analysis to put the branch-variables into the tree. */
void HLTEgamma::setup(const edm::ParameterSet& pSet, TTree* HltTree) {

  CandIso_                              = pSet.getParameter<edm::InputTag> ("CandIso");
  CandNonIso_                           = pSet.getParameter<edm::InputTag> ("CandNonIso");
  EcalIso_                              = pSet.getParameter<edm::InputTag> ("EcalIso");
  EcalNonIso_                           = pSet.getParameter<edm::InputTag> ("EcalNonIso");
  HcalIsoPho_                           = pSet.getParameter<edm::InputTag> ("HcalIsoPho");
  HcalNonIsoPho_                        = pSet.getParameter<edm::InputTag> ("HcalNonIsoPho");
  IsoPhoTrackIsol_                      = pSet.getParameter<edm::InputTag> ("IsoPhoTrackIsol");
  NonIsoPhoTrackIsol_                   = pSet.getParameter<edm::InputTag> ("NonIsoPhoTrackIsol");
  IsoElectronTag_                       = pSet.getParameter<edm::InputTag> ("IsoElectrons");
  NonIsoElectronTag_                    = pSet.getParameter<edm::InputTag> ("NonIsoElectrons");
  IsoEleHcalTag_                        = pSet.getParameter<edm::InputTag> ("HcalIsoEle");
  NonIsoEleHcalTag_                     = pSet.getParameter<edm::InputTag> ("HcalNonIsoEle");
  IsoEleTrackIsolTag_                   = pSet.getParameter<edm::InputTag> ("IsoEleTrackIsol");
  NonIsoEleTrackIsolTag_                = pSet.getParameter<edm::InputTag> ("NonIsoEleTrackIsol");
  IsoElectronLargeWindowsTag_           = pSet.getParameter<edm::InputTag> ("IsoElectronsLargeWindows");
  NonIsoElectronLargeWindowsTag_        = pSet.getParameter<edm::InputTag> ("NonIsoElectronsLargeWindows");
  IsoEleTrackIsolLargeWindowsTag_       = pSet.getParameter<edm::InputTag> ("IsoEleTrackIsolLargeWindows");
  NonIsoEleTrackIsolLargeWindowsTag_    = pSet.getParameter<edm::InputTag> ("NonIsoEleTrackIsolLargeWindows");
  L1IsoPixelSeedsTag_                   = pSet.getParameter<edm::InputTag> ("PixelSeedL1Iso");
  L1NonIsoPixelSeedsTag_                = pSet.getParameter<edm::InputTag> ("PixelSeedL1NonIso");
  L1IsoPixelSeedsLargeWindowsTag_       = pSet.getParameter<edm::InputTag> ("PixelSeedL1IsoLargeWindows");
  L1NonIsoPixelSeedsLargeWindowsTag_    = pSet.getParameter<edm::InputTag> ("PixelSeedL1NonIsoLargeWindows");

  const int kMaxEl = 10000;
  elpt              = new float[kMaxEl];
  elphi             = new float[kMaxEl];
  eleta             = new float[kMaxEl];
  elet              = new float[kMaxEl];
  ele               = new float[kMaxEl];
  const int kMaxPhot = 10000;
  photonpt          = new float[kMaxPhot];
  photonphi         = new float[kMaxPhot];
  photoneta         = new float[kMaxPhot];
  photonet          = new float[kMaxPhot];
  photone           = new float[kMaxPhot];
  const int kMaxhPhot = 500;
  hphotet           = new float[kMaxhPhot];
  hphoteta          = new float[kMaxhPhot];
  hphotphi          = new float[kMaxhPhot];
  hphoteiso         = new float[kMaxhPhot];
  hphothiso         = new float[kMaxhPhot];
  hphottiso         = new float[kMaxhPhot];
  hphotl1iso        = new int[kMaxhPhot];
  const int kMaxhEle = 500;
  heleet            = new float[kMaxhEle];
  heleeta           = new float[kMaxhEle];
  helephi           = new float[kMaxhEle];
  heleE             = new float[kMaxhEle];
  helep             = new float[kMaxhEle];
  helehiso          = new float[kMaxhEle];
  heletiso          = new float[kMaxhEle];
  helel1iso         = new int[kMaxhEle];
  helePixelSeeds    = new int[kMaxhEle];
  heleNewSC         = new int[kMaxhEle];
  const int kMaxhEleLW = 500;
  heleetLW          = new float[kMaxhEleLW];
  heleetaLW         = new float[kMaxhEleLW];
  helephiLW         = new float[kMaxhEleLW];
  heleELW           = new float[kMaxhEleLW];
  helepLW           = new float[kMaxhEleLW];
  helehisoLW        = new float[kMaxhEleLW];
  heletisoLW        = new float[kMaxhEleLW];
  helel1isoLW       = new int[kMaxhEleLW];
  helePixelSeedsLW  = new int[kMaxhEleLW];
  heleNewSCLW       = new int[kMaxhEleLW];

  // Egamma-specific branches of the tree
  HltTree->Branch("NrecoElec",          & nele,             "NrecoElec/I");
  HltTree->Branch("recoElecPt",         elpt,               "recoElecPt[NrecoElec]/F");
  HltTree->Branch("recoElecPhi",        elphi,              "recoElecPhi[NrecoElec]/F");
  HltTree->Branch("recoElecEta",        eleta,              "recoElecEta[NrecoElec]/F");
  HltTree->Branch("recoElecEt",         elet,               "recoElecEt[NrecoElec]/F");
  HltTree->Branch("recoElecE",          ele,                "recoElecE[NrecoElec]/F");
  HltTree->Branch("NrecoPhot",          &nphoton,           "NrecoPhot/I");
  HltTree->Branch("recoPhotPt",         photonpt,           "recoPhotPt[NrecoPhot]/F");
  HltTree->Branch("recoPhotPhi",        photonphi,          "recoPhotPhi[NrecoPhot]/F");
  HltTree->Branch("recoPhotEta",        photoneta,          "recoPhotEta[NrecoPhot]/F");
  HltTree->Branch("recoPhotEt",         photonet,           "recoPhotEt[NrecoPhot]/F");
  HltTree->Branch("recoPhotE",          photone,            "recoPhotE[NrecoPhot]/F");
  HltTree->Branch("NohPhot",            & nhltgam,          "NohPhot/I");
  HltTree->Branch("ohPhotEt",           hphotet,            "ohPhotEt[NohPhot]/F");
  HltTree->Branch("ohPhotEta",          hphoteta,           "ohPhotEta[NohPhot]/F");
  HltTree->Branch("ohPhotPhi",          hphotphi,           "ohPhotPhi[NohPhot]/F");
  HltTree->Branch("ohPhotEiso",         hphoteiso,          "ohPhotEiso[NohPhot]/F");
  HltTree->Branch("ohPhotHiso",         hphothiso,          "ohPhotHiso[NohPhot]/F");
  HltTree->Branch("ohPhotTiso",         hphottiso,          "ohPhotTiso[NohPhot]/F");
  HltTree->Branch("ohPhotL1iso",        hphotl1iso,         "ohPhotL1iso[NohPhot]/I");
  HltTree->Branch("NohEle",             & nhltele,          "NohEle/I");
  HltTree->Branch("ohEleEt",            heleet,             "ohEleEt[NohEle]/F");
  HltTree->Branch("ohEleEta",           heleeta,            "ohEleEta[NohEle]/F");
  HltTree->Branch("ohElePhi",           helephi,            "ohElePhi[NohEle]/F");
  HltTree->Branch("ohEleE",             heleE,              "ohEleE[NohEle]/F");
  HltTree->Branch("ohEleP",             helep,              "ohEleP[NohEle]/F");
  HltTree->Branch("ohEleHiso",          helehiso,           "ohEleHiso[NohEle]/F");
  HltTree->Branch("ohEleTiso",          heletiso,           "ohEleTiso[NohEle]/F");
  HltTree->Branch("ohEleL1iso",         helel1iso,          "ohEleLiso[NohEle]/I");
  HltTree->Branch("ohElePixelSeeds",    helePixelSeeds,     "ohElePixelSeeds[NohEle]/I");
  HltTree->Branch("ohEleNewSC",         heleNewSC,          "ohEleNewSC[NohEle]/I");
  HltTree->Branch("NohEleLW",           & nhlteleLW,        "NohEleLW/I");
  HltTree->Branch("ohEleEtLW",          heleetLW,           "ohEleEtLW[NohEleLW]/F");
  HltTree->Branch("ohEleEtaLW",         heleetaLW,          "ohEleEtaLW[NohEleLW]/F");
  HltTree->Branch("ohElePhiLW",         helephiLW,          "ohElePhiLW[NohEleLW]/F");
  HltTree->Branch("ohEleELW",           heleELW,            "ohEleELW[NohEleLW]/F");
  HltTree->Branch("ohElePLW",           helepLW,            "ohElePLW[NohEleLW]/F");
  HltTree->Branch("ohEleHisoLW",        helehisoLW,         "ohEleHisoLW[NohEleLW]/F");
  HltTree->Branch("ohEleTisoLW",        heletisoLW,         "ohEleTisoLW[NohEleLW]/F");
  HltTree->Branch("ohEleL1isoLW",       helel1isoLW,        "ohEleLisoLW[NohEleLW]/I");
  HltTree->Branch("ohElePixelSeedsLW",  helePixelSeedsLW,   "ohElePixelSeedsLW[NohEleLW]/I");
  HltTree->Branch("ohEleNewSCLW",       heleNewSCLW,        "ohEleNewSCLW[NohEleLW]/I");
}

/* **Analyze the event** */
void HLTEgamma::analyze(edm::Event const& event, edm::EventSetup const& iSetup,
                        const reco::GsfElectronCollection * electrons,
                        const reco::PhotonCollection      * photons,
                        TTree* HltTree) {

  edm::Handle<reco::ElectronCollection>                 electronIsoHandle;
  edm::Handle<reco::ElectronCollection>                 electronIsoHandleLW;
  edm::Handle<reco::ElectronCollection>                 electronNonIsoHandle;
  edm::Handle<reco::ElectronCollection>                 electronNonIsoHandleLW;
  edm::Handle<reco::ElectronIsolationMap>               NonIsoTrackEleIsolMap;
  edm::Handle<reco::ElectronIsolationMap>               NonIsoTrackEleIsolMapLW;
  edm::Handle<reco::ElectronIsolationMap>               TrackEleIsolMap;
  edm::Handle<reco::ElectronIsolationMap>               TrackEleIsolMapLW;
  edm::Handle<reco::ElectronPixelSeedCollection>        L1IsoPixelSeedsMap;
  edm::Handle<reco::ElectronPixelSeedCollection>        L1IsoPixelSeedsMapLW;
  edm::Handle<reco::ElectronPixelSeedCollection>        L1NonIsoPixelSeedsMap;
  edm::Handle<reco::ElectronPixelSeedCollection>        L1NonIsoPixelSeedsMapLW;
  edm::Handle<reco::RecoEcalCandidateCollection>        recoIsolecalcands;
  edm::Handle<reco::RecoEcalCandidateCollection>        recoNonIsolecalcands;
  edm::Handle<reco::RecoEcalCandidateIsolationMap>      EcalIsolMap;
  edm::Handle<reco::RecoEcalCandidateIsolationMap>      EcalNonIsolMap;
  edm::Handle<reco::RecoEcalCandidateIsolationMap>      HcalEleIsolMap;
  edm::Handle<reco::RecoEcalCandidateIsolationMap>      HcalIsolMap;
  edm::Handle<reco::RecoEcalCandidateIsolationMap>      HcalNonIsolMap;
  edm::Handle<reco::RecoEcalCandidateIsolationMap>      TrackIsolMap;
  edm::Handle<reco::RecoEcalCandidateIsolationMap>      TrackNonIsolMap;

  event.getByLabel(CandIso_,                           recoIsolecalcands);
  event.getByLabel(CandNonIso_,                        recoNonIsolecalcands);
  event.getByLabel(EcalIso_,                           EcalIsolMap);
  event.getByLabel(EcalNonIso_,                        EcalNonIsolMap);
  event.getByLabel(HcalIsoPho_,                        HcalIsolMap);
  event.getByLabel(HcalNonIsoPho_,                     HcalNonIsolMap);
  event.getByLabel(IsoElectronLargeWindowsTag_,        electronIsoHandleLW);
  event.getByLabel(IsoElectronTag_,                    electronIsoHandle);
  event.getByLabel(IsoEleHcalTag_,                     HcalEleIsolMap);
  event.getByLabel(IsoEleTrackIsolLargeWindowsTag_,    TrackEleIsolMapLW);
  event.getByLabel(IsoEleTrackIsolTag_,                TrackEleIsolMap);
  event.getByLabel(IsoPhoTrackIsol_,                   TrackIsolMap);
  event.getByLabel(L1IsoPixelSeedsLargeWindowsTag_,    L1IsoPixelSeedsMapLW);
  event.getByLabel(L1IsoPixelSeedsTag_,                L1IsoPixelSeedsMap);
  event.getByLabel(L1NonIsoPixelSeedsLargeWindowsTag_, L1NonIsoPixelSeedsMapLW);
  event.getByLabel(L1NonIsoPixelSeedsTag_,             L1NonIsoPixelSeedsMap);
  event.getByLabel(NonIsoElectronLargeWindowsTag_,     electronNonIsoHandleLW);
  event.getByLabel(NonIsoElectronTag_,                 electronNonIsoHandle);
  event.getByLabel(NonIsoEleHcalTag_,                  HcalEleIsolMap);
  event.getByLabel(NonIsoEleTrackIsolLargeWindowsTag_, NonIsoTrackEleIsolMapLW);
  event.getByLabel(NonIsoEleTrackIsolTag_,             NonIsoTrackEleIsolMap);
  event.getByLabel(NonIsoPhoTrackIsol_,                TrackNonIsolMap);

  if (electrons) {
    GsfElectronCollection myelectrons( electrons->begin(), electrons->end() );
    nele = myelectrons.size();
    std::sort(myelectrons.begin(), myelectrons.end(), EtGreater());
    int iel = 0;
    for (GsfElectronCollection::const_iterator i = myelectrons.begin(); i != myelectrons.end(); i++) {
      elpt[iel]  = i->pt();
      elphi[iel] = i->phi();
      eleta[iel] = i->eta();
      elet[iel]  = i->et();
      ele[iel]   = i->energy();
      iel++;
    }
  }
  else {nele = 0;}

  if (photons) {
    PhotonCollection myphotons(* photons);
    nphoton = myphotons.size();
    std::sort(myphotons.begin(), myphotons.end(), EtGreater());
    int ipho = 0;
    for (PhotonCollection::const_iterator i = myphotons.begin(); i!= myphotons.end(); i++) {
      photonpt[ipho] = i->pt();
      photonphi[ipho] = i->phi();
      photoneta[ipho] = i->eta();
      photonet[ipho] = i->et();
      photone[ipho] = i->energy();
      ipho++;
    }
  }
  else {nphoton = 0;}

  /////////////////////////////// Open-HLT Egammas ///////////////////////////////

  theHLTPhotons.clear();
  MakeL1IsolatedPhotons(
      recoIsolecalcands.product(), 
      EcalIsolMap.product(), 
      HcalIsolMap.product(), 
      TrackIsolMap.product());
  MakeL1NonIsolatedPhotons(
      recoNonIsolecalcands.product(), 
      EcalNonIsolMap.product(), 
      HcalNonIsolMap.product(), 
      TrackNonIsolMap.product());
  std::sort(theHLTPhotons.begin(), theHLTPhotons.end(), EtGreater());
  nhltgam = theHLTPhotons.size();
  for (int u = 0; u < nhltgam; u++) {
    hphotet[u]    = theHLTPhotons[u].Et;
    hphoteta[u]   = theHLTPhotons[u].eta;
    hphotphi[u]   = theHLTPhotons[u].phi;
    hphoteiso[u]  = theHLTPhotons[u].ecalIsol;
    hphothiso[u]  = theHLTPhotons[u].hcalIsol;
    hphottiso[u]  = theHLTPhotons[u].trackIsol;
    hphotl1iso[u] = theHLTPhotons[u].L1Isolated;
  }

  theHLTElectrons.clear();
  MakeL1IsolatedElectrons(
      electronIsoHandle.product(), 
      recoIsolecalcands.product(), 
      HcalEleIsolMap.product(), 
      L1IsoPixelSeedsMap.product(), 
      TrackEleIsolMap.product());
  MakeL1NonIsolatedElectrons(
      electronNonIsoHandle.product(), 
      recoNonIsolecalcands.product(), 
      HcalEleIsolMap.product(), 
      L1NonIsoPixelSeedsMap.product(), 
      NonIsoTrackEleIsolMap.product());
  std::sort(theHLTElectrons.begin(), theHLTElectrons.end(), EtGreater());
  nhltele = theHLTElectrons.size();
  for (int u = 0; u < nhltele; u++) {
    heleet[u]         = theHLTElectrons[u].Et;
    heleeta[u]        = theHLTElectrons[u].eta;
    helephi[u]        = theHLTElectrons[u].phi;
    heleE[u]          = theHLTElectrons[u].E;
    helep[u]          = theHLTElectrons[u].p;
    helehiso[u]       = theHLTElectrons[u].hcalIsol;
    helePixelSeeds[u] = theHLTElectrons[u].pixelSeeds;
    heletiso[u]       = theHLTElectrons[u].trackIsol;
    helel1iso[u]      = theHLTElectrons[u].L1Isolated;
    heleNewSC[u]      = theHLTElectrons[u].newSC;
  }

  theHLTElectronsLargeWindows.clear();
  MakeL1IsolatedElectronsLargeWindows(
      electronIsoHandleLW.product(), 
      recoIsolecalcands.product(), 
      HcalEleIsolMap.product(), 
      L1IsoPixelSeedsMapLW.product(), 
      TrackEleIsolMapLW.product());
  MakeL1NonIsolatedElectronsLargeWindows(
      electronNonIsoHandleLW.product(), 
      recoNonIsolecalcands.product(), 
      HcalEleIsolMap.product(), 
      L1NonIsoPixelSeedsMapLW.product(), 
      NonIsoTrackEleIsolMapLW.product());
  std::sort(theHLTElectronsLargeWindows.begin(), theHLTElectronsLargeWindows.end(), EtGreater());
  nhlteleLW = theHLTElectronsLargeWindows.size();
  for (int u = 0; u < nhltele; u++) {
    heleetLW[u]         = theHLTElectronsLargeWindows[u].Et;
    heleetaLW[u]        = theHLTElectronsLargeWindows[u].eta;
    helephiLW[u]        = theHLTElectronsLargeWindows[u].phi;
    heleELW[u]          = theHLTElectronsLargeWindows[u].E;
    helepLW[u]          = theHLTElectronsLargeWindows[u].p;
    helehisoLW[u]       = theHLTElectronsLargeWindows[u].hcalIsol;
    helePixelSeedsLW[u] = theHLTElectronsLargeWindows[u].pixelSeeds;
    heletisoLW[u]       = theHLTElectronsLargeWindows[u].trackIsol;
    helel1isoLW[u]      = theHLTElectronsLargeWindows[u].L1Isolated;
    heleNewSCLW[u]      = theHLTElectronsLargeWindows[u].newSC;
  }

}


void HLTEgamma::MakeL1IsolatedPhotons(
    const reco::RecoEcalCandidateCollection   * recoIsolecalcands,
    const reco::RecoEcalCandidateIsolationMap * EcalIsolMap,
    const reco::RecoEcalCandidateIsolationMap * HcalIsolMap,
    const reco::RecoEcalCandidateIsolationMap * TrackIsolMap)
{
  string EGerrMsg("");
  bool foundCand     = true;
  bool foundEcalIMap = true;
  bool foundHcalIMap = true;
  bool foundTckIMap  = true;

  if (not recoIsolecalcands) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol eg candidate";
    foundCand = false;
  }
  if (not EcalIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No Ecal isol map";
    foundEcalIMap = false;
  }
  if (not HcalIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No Hcal isol photon map";
    foundHcalIMap = false;
  }
  if (not TrackIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No Track isol photon map";
    foundTckIMap = false;
  }

  // Iterator to the isolation-map
  reco::RecoEcalCandidateIsolationMap::const_iterator mapi;

  if (foundCand) {
    // loop over SuperCluster and fill the HLTPhotons
    for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand = recoIsolecalcands->begin();
         recoecalcand!= recoIsolecalcands->end(); recoecalcand++) {

      myHLTPhoton pho;
      pho.ecalIsol   = -999;
      pho.hcalIsol   = -999;
      pho.trackIsol  = -999;
      pho.L1Isolated = true;
      pho.Et         = recoecalcand->et();
      pho.eta        = recoecalcand->eta();
      pho.phi        = recoecalcand->phi();

      // Method to get the reference to the candidate
      reco::RecoEcalCandidateRef ref = reco::RecoEcalCandidateRef(recoIsolecalcands, distance(recoIsolecalcands->begin(), recoecalcand));

      // First/Second member of the Map: Ref-to-Candidate(mapi)/Isolation(->val)
      // fill the ecal Isolation
      if (foundEcalIMap) {
        mapi = (*EcalIsolMap).find(ref);
        if (mapi !=(*EcalIsolMap).end()) { pho.ecalIsol = mapi->val;}
      }
      // fill the hcal Isolation
      if (foundHcalIMap) {
        mapi = (*HcalIsolMap).find(ref);
        if (mapi !=(*HcalIsolMap).end()) { pho.hcalIsol = mapi->val;}
      }
      // fill the track Isolation
      if (foundTckIMap) {
        mapi = (*TrackIsolMap).find(ref);
        if (mapi !=(*TrackIsolMap).end()) { pho.trackIsol = mapi->val;}
      }

      // store the photon into the vector
      theHLTPhotons.push_back(pho);
    }
  }
}

void HLTEgamma::MakeL1NonIsolatedPhotons(
    const reco::RecoEcalCandidateCollection   * recoNonIsolecalcands,
    const reco::RecoEcalCandidateIsolationMap * EcalNonIsolMap,
    const reco::RecoEcalCandidateIsolationMap * HcalNonIsolMap, 
    const reco::RecoEcalCandidateIsolationMap * TrackNonIsolMap)
{
  string EGerrMsg("");
  bool foundCand     = true;
  bool foundEcalIMap = true;
  bool foundHcalIMap = true;
  bool foundTckIMap  = true;

  if (not recoNonIsolecalcands) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No non-isol eg candidate";
    foundCand = false;
  }
  if (not EcalNonIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No Ecal non-isol map";
    foundEcalIMap = false;
  }
  if (not HcalNonIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No Hcal non-isol photon map";
    foundHcalIMap = false;
  }
  if (not TrackNonIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No Track non-isol photon map";
    foundTckIMap = false;
  }

  reco::RecoEcalCandidateIsolationMap::const_iterator mapi;

  if (foundCand) {
    for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand = recoNonIsolecalcands->begin();
         recoecalcand!= recoNonIsolecalcands->end(); recoecalcand++) { 
      // loop over SuperCluster and fill the HLTPhotons
      myHLTPhoton pho;
      pho.ecalIsol   = -999;
      pho.hcalIsol   = -999;
      pho.trackIsol  = -999;
      pho.L1Isolated = false;
      pho.Et         = recoecalcand->et();
      pho.eta        = recoecalcand->eta();
      pho.phi        = recoecalcand->phi();

      reco::RecoEcalCandidateRef ref = reco::RecoEcalCandidateRef(recoNonIsolecalcands, distance(recoNonIsolecalcands->begin(), recoecalcand));

      // fill the ecal Isolation
      if (foundEcalIMap) {
        mapi = (*EcalNonIsolMap).find(ref);
        if (mapi !=(*EcalNonIsolMap).end()) { pho.ecalIsol = mapi->val;}
      }
      // fill the hcal Isolation
      if (foundHcalIMap) {
        mapi = (*HcalNonIsolMap).find(ref);
        if (mapi !=(*HcalNonIsolMap).end()) { pho.hcalIsol = mapi->val;}
      }
      // fill the track Isolation
      if (foundTckIMap) {
        mapi = (*TrackNonIsolMap).find(ref);
        if (mapi !=(*TrackNonIsolMap).end()) { pho.trackIsol = mapi->val;}
      }

      // store the photon into the vector
      theHLTPhotons.push_back(pho);
    }
  }
}

void HLTEgamma::MakeL1IsolatedElectrons(
    const reco::ElectronCollection            * electronIsoHandle,
    const reco::RecoEcalCandidateCollection   * recoIsolecalcands,
    const reco::RecoEcalCandidateIsolationMap * HcalEleIsolMap,
    const reco::ElectronPixelSeedCollection   * L1IsoPixelSeedsMap,
    const reco::ElectronIsolationMap          * TrackEleIsolMap)
{
  // if there are electrons, then the isolation maps and the SC should be in the event; if not it is an error

  string EGerrMsg("");
  bool foundCand       = true;
  bool foundHandle     = true;
  bool foundHcalIMap   = true;
  bool foundPixelSCMap = true;
  bool foundTckIMap    = true;

  if (not recoIsolecalcands) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol eg candidate";
    foundCand = false;
  }
  if (not electronIsoHandle) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol electron";
    foundHandle = false;
  }
  if (not HcalEleIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol Hcal electron";
    foundHcalIMap = false;
  }
  if (not L1IsoPixelSeedsMap) {
    foundPixelSCMap = false;
    EGerrMsg = EGerrMsg + "  HLTEgamma: No pixelSeed-SC association map for electron";
  }
  if (not TrackEleIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol Track electron";
    foundTckIMap = false;
  }

  if (foundCand) {
    for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand = recoIsolecalcands->begin();
         recoecalcand!= recoIsolecalcands->end(); recoecalcand++) {
      // get the ref to the SC:
      reco::RecoEcalCandidateRef ref = reco::RecoEcalCandidateRef(recoIsolecalcands, distance(recoIsolecalcands->begin(), recoecalcand));
      reco::SuperClusterRef recrSC = ref->superCluster();
      //reco::SuperClusterRef recrSC = recoecalcand->superCluster();

      myHLTElectron ele;
      ele.hcalIsol   = -999;
      ele.trackIsol  = -999;
      ele.L1Isolated = true;
      ele.p          = -999;
      ele.pixelSeeds = -999;
      ele.newSC      = true;
      ele.Et         = recoecalcand->et();
      ele.eta        = recoecalcand->eta();
      ele.phi        = recoecalcand->phi();
      ele.E          = recrSC->energy();

      // fill the hcal Isolation
      if (foundHcalIMap) {
        //reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( reco::RecoEcalCandidateRef(recoIsolecalcands, distance(recoIsolecalcands->begin(), recoecalcand)) );
        reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( ref );
        if (mapi !=(*HcalEleIsolMap).end()) { ele.hcalIsol = mapi->val; }
      }
      // look if the SC has associated pixelSeeds
      int nmatch = 0;

      if (foundPixelSCMap) {
        for (reco::ElectronPixelSeedCollection::const_iterator it = L1IsoPixelSeedsMap->begin();
             it != L1IsoPixelSeedsMap->end(); it++) {
          const reco::SuperClusterRef & scRef = it->superCluster();
          if (&(*recrSC) ==  &(*scRef)) { nmatch++; }
        }
      }

      ele.pixelSeeds = nmatch;

      // look if the SC was promoted to an electron:
      if (foundHandle) {
        bool FirstElectron = true;
        reco::ElectronRef electronref;
        for (reco::ElectronCollection::const_iterator iElectron = electronIsoHandle->begin();
             iElectron != electronIsoHandle->end(); iElectron++) {
          // 1) find the SC from the electron
          electronref = reco::ElectronRef(electronIsoHandle, iElectron - electronIsoHandle->begin());
          const reco::SuperClusterRef theClus = electronref->superCluster(); // SC from the electron;
          if (&(*recrSC) ==  &(*theClus)) {     // ref is the RecoEcalCandidateRef corresponding to the electron
            if (FirstElectron) {                // the first electron is stored in ele, keeping the ele.newSC = true
              FirstElectron = false;
              ele.p = electronref->track()->momentum().R();
              // fill the track Isolation
              if (foundTckIMap) {
                reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find(electronref);
                if (mapTr != (*TrackEleIsolMap).end()) { ele.trackIsol = mapTr->val; }
              }
            }
            else {
              // FirstElectron is false, i.e. the SC of this electron is common to another electron.
              // A new  myHLTElectron is inserted in the theHLTElectrons vector setting newSC = false
              myHLTElectron ele2;
              ele2.hcalIsol  = ele.hcalIsol;
              ele2.trackIsol = -999;
              ele2.Et  = ele.Et;
              ele2.eta = ele.eta;
              ele2.phi = ele.phi;
              ele2.E   = ele.E;
              ele2.L1Isolated = ele.L1Isolated;
              ele2.pixelSeeds = ele.pixelSeeds;
              ele2.newSC = false;
              ele2.p = electronref->track()->momentum().R();
              // fill the track Isolation
              if (foundTckIMap) {
                reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find( electronref);
                if (mapTr !=(*TrackEleIsolMap).end()) { ele2.trackIsol = mapTr->val;}
              }
              theHLTElectrons.push_back(ele2);
            }
          }
        } // end of loop over electrons
      } // end of if (foundHandle) {

      //store the electron into the vector
      theHLTElectrons.push_back(ele);
    } // end of loop over ecalCandidates
  } // end of if (foundCand) {
}


void HLTEgamma::MakeL1NonIsolatedElectrons(
    const reco::ElectronCollection            * electronNonIsoHandle,
    const reco::RecoEcalCandidateCollection   * recoNonIsolecalcands,
    const reco::RecoEcalCandidateIsolationMap * HcalEleIsolMap,
    const reco::ElectronPixelSeedCollection   * L1NonIsoPixelSeedsMap,
    const reco::ElectronIsolationMap          * TrackEleIsolMap)
{
  // if there are electrons, then the isolation maps and the SC should be in the event; if not it is an error

  string EGerrMsg("");
  bool foundCand       = true;
  bool foundHandle     = true;
  bool foundHcalIMap   = true;
  bool foundPixelSCMap = true;
  bool foundTckIMap    = true;

  if (not recoNonIsolecalcands) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol eg candidate";
    foundCand = false;
  }
  if (not electronNonIsoHandle) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol electron";
    foundHandle = false;
  }
  if (not HcalEleIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol Hcal electron";
    foundHcalIMap = false;
  }
  if (not L1NonIsoPixelSeedsMap) {
    foundPixelSCMap = false;
    EGerrMsg = EGerrMsg + "  HLTEgamma: No pixelSeed-SC for electron";
  }
  if (not TrackEleIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol Track electron";
    foundTckIMap = false;
  }

  if (foundCand) {
    for(reco::RecoEcalCandidateCollection::const_iterator recoecalcand = recoNonIsolecalcands->begin();
                                        recoecalcand!= recoNonIsolecalcands->end(); recoecalcand++) {
      //get the ref to the SC:
      reco::RecoEcalCandidateRef ref = reco::RecoEcalCandidateRef(recoNonIsolecalcands, distance(recoNonIsolecalcands->begin(), recoecalcand));
      reco::SuperClusterRef recrSC = ref->superCluster();
      //reco::SuperClusterRef recrSC = recoecalcand->superCluster();

      myHLTElectron ele;
      ele.hcalIsol   = -999; 
      ele.trackIsol  = -999;
      ele.L1Isolated = false; 
      ele.p          = -999;
      ele.pixelSeeds = -999; 
      ele.newSC      = true;
      ele.Et         = recoecalcand->et();
      ele.eta        = recoecalcand->eta();
      ele.phi        = recoecalcand->phi();
      ele.E          = recrSC->energy();

      // fill the hcal Isolation
      if (foundHcalIMap) {
        // reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( reco::RecoEcalCandidateRef(recoNonIsolecalcands, distance(recoNonIsolecalcands->begin(), recoecalcand)) );
        reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( ref );
        if (mapi !=(*HcalEleIsolMap).end()) {ele.hcalIsol = mapi->val;}
      }
      // look if the SC has associated pixelSeeds
      int nmatch = 0;

      if (foundPixelSCMap) {
        for (reco::ElectronPixelSeedCollection::const_iterator it = L1NonIsoPixelSeedsMap->begin();
             it != L1NonIsoPixelSeedsMap->end(); it++) {
          const reco::SuperClusterRef & scRef = it->superCluster();
          if (&(*recrSC) == &(*scRef)) { nmatch++;}
        }
      }

      ele.pixelSeeds = nmatch;

      // look if the SC was promoted to an electron:
      if (foundHandle) {
        bool FirstElectron = true;
        reco::ElectronRef electronref;
        for(reco::ElectronCollection::const_iterator iElectron = electronNonIsoHandle->begin(); iElectron !=
              electronNonIsoHandle->end();iElectron++) {
          // 1) find the SC from the electron
          electronref = reco::ElectronRef(electronNonIsoHandle, iElectron - electronNonIsoHandle->begin());
          const reco::SuperClusterRef theClus = electronref->superCluster(); //SC from the electron;
          if (&(*recrSC) ==  &(*theClus)) { // ref is the RecoEcalCandidateRef corresponding to the electron
            if (FirstElectron) { //the first electron is stored in ele, keeping the ele.newSC = true
              FirstElectron = false;
              ele.p = electronref->track()->momentum().R();
              // fill the track Isolation
              if (foundTckIMap) {
                reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find( electronref);
                if (mapTr !=(*TrackEleIsolMap).end()) { ele.trackIsol = mapTr->val;}
              }
            } else {
              // FirstElectron is false, i.e. the SC of this electron is common to another electron.
              // A new myHLTElectron is inserted in the theHLTElectrons vector setting newSC = false
              myHLTElectron ele2;
              ele2.hcalIsol   = ele.hcalIsol;
              ele2.trackIsol  =-999;
              ele2.Et         = ele.Et;
              ele2.eta        = ele.eta;
              ele2.phi        = ele.phi;
              ele2.E          = ele.E;
              ele2.L1Isolated = ele.L1Isolated;
              ele2.pixelSeeds = ele.pixelSeeds;
              ele2.newSC      = false;
              ele2.p          = electronref->track()->momentum().R();
              // fill the track Isolation
              if (foundTckIMap) {
                reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find( electronref);
                if (mapTr !=(*TrackEleIsolMap).end()) { ele2.trackIsol = mapTr->val;}
              }
              theHLTElectrons.push_back(ele2);
            }
          }
        } // end of loop over electrons
      } // end of if (foundHandle) {

      // store the electron into the vector
      theHLTElectrons.push_back(ele);
    } // end of loop over ecalCandidates
  } // end of if (foundCand) {

}


void HLTEgamma::MakeL1IsolatedElectronsLargeWindows(
    const reco::ElectronCollection            * electronIsoHandle,
    const reco::RecoEcalCandidateCollection   * recoIsolecalcands,
    const reco::RecoEcalCandidateIsolationMap * HcalEleIsolMap,
    const reco::ElectronPixelSeedCollection   * L1IsoPixelSeedsMap,
    const reco::ElectronIsolationMap          * TrackEleIsolMap)
{
  // if there are electrons, then the isolation maps and the SC should be in the event; if not it is an error

  string EGerrMsg("");
  bool foundCand       = true;
  bool foundHandle     = true;
  bool foundHcalIMap   = true;
  bool foundPixelSCMap = true; 
  bool foundTckIMap    = true;

  if (not recoIsolecalcands) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol eg candidate";
    foundCand = false;
  }
  if (not electronIsoHandle) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol electron";
    foundHandle = false;
  }
  if (not HcalEleIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol Hcal electron";
    foundHcalIMap = false;
  }
  if (not L1IsoPixelSeedsMap) {
    foundPixelSCMap = false;
    EGerrMsg = EGerrMsg + "  HLTEgamma: No pixelSeed-SC for electron";
  }
  if (not TrackEleIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol Track electron";
    foundTckIMap = false;
  }

  if (foundCand) {
    for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand = recoIsolecalcands->begin();
         recoecalcand!= recoIsolecalcands->end(); recoecalcand++) {
      // get the ref to the SC:
      reco::RecoEcalCandidateRef ref = reco::RecoEcalCandidateRef(recoIsolecalcands, distance(recoIsolecalcands->begin(), recoecalcand));
      reco::SuperClusterRef recrSC = ref->superCluster();
      //reco::SuperClusterRef recrSC = recoecalcand->superCluster();

      myHLTElectron ele;
      ele.hcalIsol   = -999; 
      ele.trackIsol  = -999;
      ele.L1Isolated = true; 
      ele.p          = -999;
      ele.pixelSeeds = -999; 
      ele.newSC      = true;
      ele.Et         = recoecalcand->et();
      ele.eta        = recoecalcand->eta();
      ele.phi        = recoecalcand->phi();
      ele.E          = recrSC->energy();

      // fill the hcal Isolation
      if (foundHcalIMap) {
        //reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( reco::RecoEcalCandidateRef(recoIsolecalcands, distance(recoIsolecalcands->begin(), recoecalcand)) );
        reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( ref );
        if (mapi !=(*HcalEleIsolMap).end()) {ele.hcalIsol = mapi->val;}
      }
      // look if the SC has associated pixelSeeds
      int nmatch = 0;

      if (foundPixelSCMap) {
        for(reco::ElectronPixelSeedCollection::const_iterator it = L1IsoPixelSeedsMap->begin();
            it != L1IsoPixelSeedsMap->end(); it++) {
          const reco::SuperClusterRef & scRef = it->superCluster();
          if (&(*recrSC) == &(*scRef)) { nmatch++;}
        }
      }

      ele.pixelSeeds = nmatch;

      // look if the SC was promoted to an electron:
      if (foundHandle) {
        bool FirstElectron = true;
        reco::ElectronRef electronref;
        for(reco::ElectronCollection::const_iterator iElectron = electronIsoHandle->begin(); iElectron !=
              electronIsoHandle->end();iElectron++) {
          // 1) find the SC from the electron
          electronref = reco::ElectronRef(electronIsoHandle, iElectron - electronIsoHandle->begin());
          const reco::SuperClusterRef theClus = electronref->superCluster(); //SC from the electron;
          if (&(*recrSC) == &(*theClus)) { // ref is the RecoEcalCandidateRef corresponding to the electron
            if (FirstElectron) { //the first electron is stored in ele, keeping the ele.newSC = true
              FirstElectron = false;
              ele.p = electronref->track()->momentum().R();
              // fill the track Isolation
              if (foundTckIMap) {
                reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find( electronref);
                if (mapTr !=(*TrackEleIsolMap).end()) { ele.trackIsol = mapTr->val;}
              }
            }
            else {
              // FirstElectron is false, i.e. the SC of this electron is common to another electron.
              // A new  myHLTElectron is inserted in the theHLTElectrons vector setting newSC = false
              myHLTElectron ele2;
              ele2.hcalIsol   = ele.hcalIsol;
              ele2.trackIsol  = -999;
              ele2.Et         = ele.Et;
              ele2.eta        = ele.eta;
              ele2.phi        = ele.phi;
              ele2.E          = ele.E;
              ele2.L1Isolated = ele.L1Isolated;
              ele2.pixelSeeds = ele.pixelSeeds;
              ele2.newSC      = false;
              ele2.p          = electronref->track()->momentum().R();
              // fill the track Isolation
              if (foundTckIMap) {
                reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find( electronref);
                if (mapTr !=(*TrackEleIsolMap).end()) { ele2.trackIsol = mapTr->val;}
              }
              theHLTElectronsLargeWindows.push_back(ele2);
            }
          }
        } // end of loop over electrons
      } // end of if (foundHandle) {

      // store the electron into the vector
      theHLTElectronsLargeWindows.push_back(ele);
    } // end of loop over ecalCandidates
  } // end of if (foundCand) {
}


void HLTEgamma::MakeL1NonIsolatedElectronsLargeWindows(
    const reco::ElectronCollection            * electronNonIsoHandle,
    const reco::RecoEcalCandidateCollection   * recoNonIsolecalcands,
    const reco::RecoEcalCandidateIsolationMap * HcalEleIsolMap,
    const reco::ElectronPixelSeedCollection   * L1NonIsoPixelSeedsMap,
    const reco::ElectronIsolationMap          * TrackEleIsolMap)
{
  // if there are electrons, then the isolation maps and the SC should be in the event; if not it is an error

  string EGerrMsg("");
  bool foundCand       = true;
  bool foundHandle     = true;
  bool foundHcalIMap   = true;
  bool foundPixelSCMap = true;
  bool foundTckIMap    = true;

  if (not recoNonIsolecalcands) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol eg candidate";
    foundCand = false;
  }
  if (not electronNonIsoHandle) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol electron";
    foundHandle = false;
  }
  if (not HcalEleIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol Hcal electron";
    foundHcalIMap = false;
  }
  if (not L1NonIsoPixelSeedsMap) {
    foundPixelSCMap = false;
    EGerrMsg = EGerrMsg + "  HLTEgamma: No pixelSeed-SC for electron";
  }
  if (not TrackEleIsolMap) {
    EGerrMsg = EGerrMsg + "  HLTEgamma: No isol Track electron";
    foundTckIMap = false;
  }

  if (foundCand) {
    for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand = recoNonIsolecalcands->begin();
         recoecalcand!= recoNonIsolecalcands->end(); recoecalcand++) {
      //get the ref to the SC:
      reco::RecoEcalCandidateRef ref = reco::RecoEcalCandidateRef(recoNonIsolecalcands, distance(recoNonIsolecalcands->begin(), recoecalcand));
      reco::SuperClusterRef recrSC = ref->superCluster();
      //reco::SuperClusterRef recrSC = recoecalcand->superCluster();

      myHLTElectron ele;
      ele.hcalIsol   = -999; 
      ele.trackIsol  = -999;
      ele.L1Isolated = false; 
      ele.p          = -999;
      ele.pixelSeeds = -999;
      ele.newSC      = true;
      ele.Et         = recoecalcand->et();
      ele.eta        = recoecalcand->eta();
      ele.phi        = recoecalcand->phi();
      ele.E          = recrSC->energy();

      // fill the hcal Isolation
      if (foundHcalIMap) {
        //reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( reco::RecoEcalCandidateRef(recoNonIsolecalcands, distance(recoNonIsolecalcands->begin(), recoecalcand)) );
        reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( ref );
        if (mapi !=(*HcalEleIsolMap).end()) {ele.hcalIsol = mapi->val;}
      }
      // look if the SC has associated pixelSeeds
      int nmatch = 0;

      if (foundPixelSCMap) {
        for (reco::ElectronPixelSeedCollection::const_iterator it = L1NonIsoPixelSeedsMap->begin();
             it != L1NonIsoPixelSeedsMap->end(); it++) {
          const reco::SuperClusterRef & scRef = it->superCluster();
          if (&(*recrSC) ==  &(*scRef)) { nmatch++;}
        }
      }

      ele.pixelSeeds = nmatch;

      // look if the SC was promoted to an electron:
      if (foundHandle) {
        bool FirstElectron = true;
        reco::ElectronRef electronref;
        for (reco::ElectronCollection::const_iterator iElectron = electronNonIsoHandle->begin(); 
             iElectron != electronNonIsoHandle->end(); iElectron++) {
          // 1) find the SC from the electron
          electronref = reco::ElectronRef(electronNonIsoHandle, iElectron - electronNonIsoHandle->begin());
          const reco::SuperClusterRef theClus = electronref->superCluster(); //SC from the electron;
          if (&(*recrSC) ==  &(*theClus)) { // ref is the RecoEcalCandidateRef corresponding to the electron
            if (FirstElectron) { // the first electron is stored in ele, keeping the ele.newSC = true
              FirstElectron = false;
              ele.p = electronref->track()->momentum().R();
              // fill the track Isolation
              if (foundTckIMap) {
                reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find( electronref);
                if (mapTr !=(*TrackEleIsolMap).end()) { ele.trackIsol = mapTr->val;}
              }
            } else {
              // FirstElectron is false, i.e. the SC of this electron is common to another electron.
              // A new myHLTElectron is inserted in the theHLTElectrons vector setting newSC = false
              myHLTElectron ele2;
              ele2.hcalIsol   = ele.hcalIsol;
              ele2.trackIsol  = -999;
              ele2.Et         = ele.Et;
              ele2.eta        = ele.eta;
              ele2.phi        = ele.phi;
              ele2.E          = ele.E;
              ele2.L1Isolated = ele.L1Isolated;
              ele2.pixelSeeds = ele.pixelSeeds;
              ele2.newSC      = false;
              ele2.p          = electronref->track()->momentum().R();
              // fill the track Isolation
              if (foundTckIMap) {
                reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find( electronref);
                if (mapTr !=(*TrackEleIsolMap).end()) { ele2.trackIsol = mapTr->val;}
              }
              theHLTElectronsLargeWindows.push_back(ele2);
            }
          }
        } // end of loop over electrons
      } // end of if (foundHandle) {

      // store the electron into the vector
      theHLTElectronsLargeWindows.push_back(ele);
    } // end of loop over ecalCandidates
  } // end of if (foundCand) {
}

