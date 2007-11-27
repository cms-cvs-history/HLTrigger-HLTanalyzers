#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <stdlib.h>
#include <string.h>

#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "HLTrigger/HLTanalyzers/interface/HLTEgamma.h"


HLTEgamma::HLTEgamma() {
  evtCounter=0;

  //set parameter defaults 
  _Monte=false;
  _Debug=false;
}

/*  Setup the analysis to put the branch-variables into the tree. */
void HLTEgamma::setup(const edm::ParameterSet& pSet, TTree* HltTree) {

  CandIso_ = pSet.getParameter<std::string> ("CandIso");
  CandNonIso_ = pSet.getParameter<std::string> ("CandNonIso");
  EcalIso_ = pSet.getParameter<std::string> ("EcalIso");
  EcalNonIso_ = pSet.getParameter<std::string> ("EcalNonIso");
  HcalIsoPho_ = pSet.getParameter<std::string> ("HcalIsoPho");
  HcalNonIsoPho_ = pSet.getParameter<std::string> ("HcalNonIsoPho");
  IsoPhoTrackIsol_ = pSet.getParameter<std::string> ("IsoPhoTrackIsol");
  NonIsoPhoTrackIsol_ = pSet.getParameter<std::string> ("NonIsoPhoTrackIsol");
  IsoElectronTag_ = pSet.getParameter<std::string> ("IsoElectrons");
  NonIsoElectronTag_ = pSet.getParameter<std::string> ("NonIsoElectrons");
  IsoEleHcalTag_ = pSet.getParameter<std::string> ("HcalIsoEle");
  NonIsoEleHcalTag_ = pSet.getParameter<std::string> ("HcalNonIsoEle");
  IsoEleTrackIsolTag_ = pSet.getParameter<std::string> ("IsoEleTrackIsol");
  NonIsoEleTrackIsolTag_ = pSet.getParameter<std::string> ("NonIsoEleTrackIsol");

  edm::ParameterSet myEmParams = pSet.getParameter<edm::ParameterSet>("RunParameters") ;
  vector<std::string> parameterNames = myEmParams.getParameterNames() ;
  
  for ( vector<std::string>::iterator iParam = parameterNames.begin();
	iParam != parameterNames.end(); iParam++ ){
    if  ( (*iParam) == "Monte" ) _Monte =  myEmParams.getParameter<bool>( *iParam );
    else if ( (*iParam) == "Debug" ) _Debug =  myEmParams.getParameter<bool>( *iParam );
  }
  
  const int kMaxEl = 10000;
  elpt = new float[kMaxEl];
  elphi = new float[kMaxEl];
  eleta = new float[kMaxEl];
  elet = new float[kMaxEl];
  ele = new float[kMaxEl];
  const int kMaxPhot = 10000;
  photonpt = new float[kMaxPhot];
  photonphi = new float[kMaxPhot];
  photoneta = new float[kMaxPhot];
  photonet = new float[kMaxPhot];
  photone = new float[kMaxPhot];
  const int kMaxhPhot = 500;
  hphotet = new float[kMaxhPhot];
  hphoteta = new float[kMaxhPhot];
  hphotphi = new float[kMaxhPhot];
  hphoteiso = new float[kMaxhPhot];
  hphothiso = new float[kMaxhPhot];
  hphottiso = new float[kMaxhPhot];
  hphotl1iso = new int[kMaxhPhot];
  const int kMaxhEle = 500;
  heleet = new float[kMaxhEle];
  heleeta = new float[kMaxhEle];
  helephi = new float[kMaxhEle];
  heleE = new float[kMaxhEle];
  helep = new float[kMaxhEle];
  helehiso = new float[kMaxhEle];
  heletiso = new float[kMaxhEle];
  helel1iso = new int[kMaxhEle];

  // Egamma-specific branches of the tree 
  HltTree->Branch("NrecoElec",&nele,"NrecoElec/I");
  HltTree->Branch("recoElecPt",elpt,"recoElecPt[NrecoElec]/F");
  HltTree->Branch("recoElecPhi",elphi,"recoElecPhi[NrecoElec]/F");
  HltTree->Branch("recoElecEta",eleta,"recoElecEta[NrecoElec]/F");
  HltTree->Branch("recoElecEt",elet,"recoElecEt[NrecoElec]/F");
  HltTree->Branch("recoElecE",ele,"recoElecE[NrecoElec]/F");
  HltTree->Branch("NrecoPhot",&nphoton,"NrecoPhot/I");
  HltTree->Branch("recoPhotPt",photonpt,"recoPhotPt[NrecoPhot]/F");
  HltTree->Branch("recoPhotPhi",photonphi,"recoPhotPhi[NrecoPhot]/F");
  HltTree->Branch("recoPhotEta",photoneta,"recoPhotEta[NrecoPhot]/F");
  HltTree->Branch("recoPhotEt",photonet,"recoPhotEt[NrecoPhot]/F");
  HltTree->Branch("recoPhotE",photone,"recoPhotE[NrecoPhot]/F");
  HltTree->Branch("NohPhot",&nhltgam,"NohPhot/I");
  HltTree->Branch("ohPhotEt",hphotet,"ohPhotEt[NohPhot]/F");
  HltTree->Branch("ohPhotEta",hphoteta,"ohPhotEta[NohPhot]/F");
  HltTree->Branch("ohPhotPhi",hphotphi,"ohPhotPhi[NohPhot]/F");
  HltTree->Branch("ohPhotEiso",hphoteiso,"ohPhotEiso[NohPhot]/F");
  HltTree->Branch("ohPhotHiso",hphothiso,"ohPhotHiso[NohPhot]/F");
  HltTree->Branch("ohPhotTiso",hphottiso,"ohPhotTiso[NohPhot]/F");
  HltTree->Branch("ohPhotL1iso",hphotl1iso,"ohPhotL1iso[NohPhot]/I");
  HltTree->Branch("NohEle",&nhltele,"NohEle/I");
  HltTree->Branch("ohEleEt",heleet,"ohEleEt[NohEle]/F");
  HltTree->Branch("ohEleEta",heleeta,"ohEleEta[NohEle]/F");
  HltTree->Branch("ohElePhi",helephi,"ohElePhi[NohEle]/F");
  HltTree->Branch("ohEleE",heleE,"ohEleE[NohEle]/F");
  HltTree->Branch("ohEleP",helep,"ohEleP[NohEle]/F");
  HltTree->Branch("ohEleHiso",helehiso,"ohEleHiso[NohEle]/F");
  HltTree->Branch("ohEleTiso",heletiso,"ohEleTiso[NohEle]/F");
  HltTree->Branch("ohEleL1iso",helel1iso,"ohEleLiso[NohEle]/I");

}

/* **Analyze the event** */
void HLTEgamma::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup,
			const ElectronCollection& Electron,
			const PhotonCollection& Photon,
			const CaloGeometry& geom,
			TTree* HltTree) {

  //std::cout << " Beginning HLTEgamma " << std::endl;

  if (&Electron) {
    ElectronCollection myelectrons;
    myelectrons=Electron;
    nele = myelectrons.size();
    std::sort(myelectrons.begin(),myelectrons.end(),EtGreater());
    typedef ElectronCollection::const_iterator ceiter;
    int iel=0;
    for (ceiter i=myelectrons.begin(); i!=myelectrons.end(); i++) {
      elpt[iel] = i->pt();
      elphi[iel] = i->phi();
      eleta[iel] = i->eta();
      elet[iel] = i->et();
      ele[iel] = i->energy();
      iel++;
    }
  }
  else {nele = 0;}

  if (&Photon) {
    PhotonCollection myphotons;
    myphotons=Photon;
    nphoton = myphotons.size();
    std::sort(myphotons.begin(),myphotons.end(),EtGreater());
    typedef PhotonCollection::const_iterator phiter;
    int ipho=0;
    for (phiter i=myphotons.begin(); i!=myphotons.end(); i++) {
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
  MakeL1IsolatedPhotons(iEvent,iSetup);
  MakeL1NonIsolatedPhotons(iEvent,iSetup);
  nhltgam = theHLTPhotons.size();
  std::sort(theHLTPhotons.begin(),theHLTPhotons.end(),EtGreater());
  for(int u=0; u<nhltgam; u++){
    hphotet[u] = theHLTPhotons[u].Et;
    hphoteta[u] = theHLTPhotons[u].eta;
    hphotphi[u] = theHLTPhotons[u].phi;
    hphoteiso[u] = theHLTPhotons[u].ecalIsol;
    hphothiso[u] = theHLTPhotons[u].hcalIsol;
    hphottiso[u] = theHLTPhotons[u].trackIsol;
    hphotl1iso[u] = theHLTPhotons[u].L1Isolated;
  }

  theHLTElectrons.clear();
  MakeL1IsolatedElectrons(iEvent,iSetup);
  MakeL1NonIsolatedElectrons(iEvent,iSetup);
  nhltele = theHLTElectrons.size();
  std::sort(theHLTElectrons.begin(),theHLTElectrons.end(),EtGreater());
  for(int u=0; u<nhltele; u++){
    heleet[u] = theHLTElectrons[u].Et;
    heleeta[u] = theHLTElectrons[u].eta;  
    helephi[u] = theHLTElectrons[u].phi;
    heleE[u] = theHLTElectrons[u].E;
    helep[u] = theHLTElectrons[u].p;
    helehiso[u] = theHLTElectrons[u].hcalIsol;
    heletiso[u] = theHLTElectrons[u].trackIsol;
    helel1iso[u] = theHLTElectrons[u].L1Isolated;
  }


}




void HLTEgamma::MakeL1IsolatedPhotons(edm::Event const& iEvent, edm::EventSetup const& iSetup){

  string EGerrMsg("");
  bool foundCand=true;bool foundEcalIMap=true;bool foundHcalIMap=true;bool foundTckIMap=true;
  edm::Handle<reco::RecoEcalCandidateCollection> recoIsolecalcands;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> EcalIsolMap,HcalIsolMap,TrackIsolMap;
  try {iEvent.getByLabel(CandIso_,recoIsolecalcands);} catch (...) {
    EGerrMsg=EGerrMsg + "  HLTEgamma: No isol eg candidate";
    foundCand=false;
  }
  try {iEvent.getByLabel(EcalIso_,EcalIsolMap);} catch (...) {
    EGerrMsg=EGerrMsg + "  HLTEgamma: No Ecal isol map";
    foundEcalIMap=false;
  }
  try {iEvent.getByLabel(HcalIsoPho_,HcalIsolMap);} catch (...) {
    EGerrMsg=EGerrMsg + "  HLTEgamma: No Hcal isol photon map";
    foundHcalIMap=false;
  }
  try {iEvent.getByLabel(IsoPhoTrackIsol_,TrackIsolMap);} catch (...) {
    EGerrMsg=EGerrMsg + "  HLTEgamma: No Track isol photon map";
    foundTckIMap=false;
  }

  // Iterator to the isolation-map
  reco::RecoEcalCandidateIsolationMap::const_iterator mapi;

  if(foundCand){
    //Loop over SuperCluster and fill the HLTPhotons
    for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand= recoIsolecalcands->begin(); 
	 recoecalcand!=recoIsolecalcands->end(); recoecalcand++) {

      myHLTPhoton pho;
      pho.ecalIsol=-999;pho.hcalIsol=-999;pho.trackIsol=-999;
      pho.L1Isolated = true;

      pho.Et=recoecalcand->et();
      pho.eta=recoecalcand->eta();
      pho.phi=recoecalcand->phi();
   
      //Method to get the reference to the candidate
      reco::RecoEcalCandidateRef ref= reco::RecoEcalCandidateRef(recoIsolecalcands,distance(recoIsolecalcands->begin(),recoecalcand));

      // First/Second member of the Map: Ref-to-Candidate(mapi)/Isolation(->val)
      //Fill the ecal Isolation
      if (foundEcalIMap){
	mapi = (*EcalIsolMap).find( ref);
	if (mapi !=(*EcalIsolMap).end()) { pho.ecalIsol=mapi->val;}
      }
      //Fill the hcal Isolation
      if (foundHcalIMap){
	mapi = (*HcalIsolMap).find( ref);
	if (mapi !=(*HcalIsolMap).end()) { pho.hcalIsol=mapi->val;}
      }
      //Fill the track Isolation
      if (foundTckIMap){
	mapi = (*TrackIsolMap).find( ref);
	if (mapi !=(*TrackIsolMap).end()) { pho.trackIsol=mapi->val;}
      }

      //store the photon into the vector
      theHLTPhotons.push_back(pho);

    }

  }

}

void HLTEgamma::MakeL1NonIsolatedPhotons(edm::Event const& iEvent, edm::EventSetup const& iSetup){

  string EGerrMsg("");
  bool foundCand=true;bool foundEcalIMap=true;bool foundHcalIMap=true;bool foundTckIMap=true;
  edm::Handle<reco::RecoEcalCandidateCollection> recoNonIsolecalcands;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> EcalNonIsolMap,HcalNonIsolMap,TrackNonIsolMap;
  try {iEvent.getByLabel(CandNonIso_,recoNonIsolecalcands);} catch (...) { 
    EGerrMsg=EGerrMsg + "  HLTEgamma: No non-isol eg candidate";
    foundCand=false;
  }
  try {iEvent.getByLabel(EcalNonIso_,EcalNonIsolMap);} catch (...) { 
    EGerrMsg=EGerrMsg + "  HLTEgamma: No Ecal non-isol map";
    foundEcalIMap=false;
  }
  try {iEvent.getByLabel(HcalNonIsoPho_,HcalNonIsolMap);} catch (...) { 
    EGerrMsg=EGerrMsg + "  HLTEgamma: No Hcal non-isol photon map";
    foundHcalIMap=false;
  }
  try {iEvent.getByLabel(NonIsoPhoTrackIsol_,TrackNonIsolMap);} catch (...) { 
    EGerrMsg=EGerrMsg + "  HLTEgamma: No Track non-isol photon map";
    foundTckIMap=false;
  }

  reco::RecoEcalCandidateIsolationMap::const_iterator mapi;

  if(foundCand){
    for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand= recoNonIsolecalcands->begin(); 
	 recoecalcand!=recoNonIsolecalcands->end(); recoecalcand++) {//Loop over SuperCluster and fill the HLTPhotons

      myHLTPhoton pho;
      pho.ecalIsol=-999;pho.hcalIsol=-999;pho.trackIsol=-999;
      pho.L1Isolated = false;
   
      pho.Et=recoecalcand->et();
      pho.eta=recoecalcand->eta();
      pho.phi=recoecalcand->phi();

      reco::RecoEcalCandidateRef ref= reco::RecoEcalCandidateRef(recoNonIsolecalcands,distance(recoNonIsolecalcands->begin(),recoecalcand));
      
      //Fill the ecal Isolation
      if (foundEcalIMap){
	mapi = (*EcalNonIsolMap).find( ref);
	if (mapi !=(*EcalNonIsolMap).end()) { pho.ecalIsol=mapi->val;}
      }
      //Fill the hcal Isolation
      if (foundHcalIMap){
	mapi = (*HcalNonIsolMap).find( ref);
	if (mapi !=(*HcalNonIsolMap).end()) { pho.hcalIsol=mapi->val;}
      }
      //Fill the track Isolation
      if (foundTckIMap){
	mapi = (*TrackNonIsolMap).find( ref);
	if (mapi !=(*TrackNonIsolMap).end()) { pho.trackIsol=mapi->val;}
      }

      //store the photon into the vector
      theHLTPhotons.push_back(pho);
      
    }

  }

}

void HLTEgamma::MakeL1IsolatedElectrons(edm::Event const& iEvent, edm::EventSetup const& iSetup){
  // If there are electrons, then the isolation maps and the SC should be in the event; if not it is an error

  string EGerrMsg("");
  bool foundCand=true;bool foundHandle=true;bool foundHcalIMap=true;bool foundTckIMap=true;
  edm::Handle<reco::ElectronCollection> electronIsoHandle;
  edm::Handle<reco::RecoEcalCandidateCollection> recoIsolecalcands;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> HcalEleIsolMap;
  edm::Handle<reco::ElectronIsolationMap> TrackEleIsolMap;
  try {iEvent.getByLabel(CandIso_,recoIsolecalcands);} catch (...) {
    EGerrMsg=EGerrMsg + "  HLTEgamma: No isol eg candidate";
    foundCand=false;
  }
  try {iEvent.getByLabel(IsoElectronTag_,electronIsoHandle);} catch (...) { 
    EGerrMsg=EGerrMsg + "  HLTEgamma: No isol electron";
    foundHandle=false;
  }
  try {iEvent.getByLabel(IsoEleHcalTag_,HcalEleIsolMap);} catch (...) { 
    EGerrMsg=EGerrMsg + "  HLTEgamma: No isol Hcal electron";
    foundHcalIMap=false;
  }
  try {iEvent.getByLabel(IsoEleTrackIsolTag_,TrackEleIsolMap);} catch (...) { 
    EGerrMsg=EGerrMsg + "  HLTEgamma: No isol Track electron";
    foundTckIMap=false;
  }

  if (foundHandle){
    for(reco::ElectronCollection::const_iterator iElectron = electronIsoHandle->begin(); iElectron != 
	electronIsoHandle->end();iElectron++){

      // 1) find the SC from the electron
      reco::ElectronRef electronref(reco::ElectronRef(electronIsoHandle,iElectron - electronIsoHandle->begin()));
      const reco::SuperClusterRef theClus = electronref->superCluster(); //SC from the electron;

      // 2) find the RecoEcalCandidate (wrapper for SC) corresponding to the SC of the electron
      reco::RecoEcalCandidateRef ref;
      if (foundCand){
	for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand= recoIsolecalcands->begin();
	     recoecalcand!=recoIsolecalcands->end(); recoecalcand++) {
	  ref = reco::RecoEcalCandidateRef(recoIsolecalcands,distance(recoIsolecalcands->begin(),recoecalcand));
	  reco::SuperClusterRef recrSC = ref->superCluster();
	  if(&(*recrSC) ==  &(*theClus)) {// ref is the RecoEcalCandidateRef corresponding to the electron
	    break;
	  }
	}
	//now ref is the RecoEcalCandidateRef correspoding to the electron

	myHLTElectron ele;
	ele.hcalIsol=-999; ele.trackIsol=-999;
	ele.L1Isolated = true;

	ele.Et=ref->et();
	ele.eta=ref->eta();
	ele.phi=ref->phi();
	ele.E = theClus->energy();
	ele.p=electronref->track()->momentum().R();

	//Fill the hcal Isolation
	if (foundHcalIMap){
	  reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleIsolMap).find( ref);
	  if(mapi !=(*HcalEleIsolMap).end()) { ele.hcalIsol=mapi->val;}
	}
	//Fill the track Isolation
	// find the
	if(foundTckIMap){
	  reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleIsolMap).find( electronref);
	  if(mapTr !=(*TrackEleIsolMap).end()) { ele.trackIsol=mapTr->val;}
	}

	//store the electron into the vector
	theHLTElectrons.push_back(ele);

      }

    }
  }

}

void HLTEgamma::MakeL1NonIsolatedElectrons(edm::Event const& iEvent, edm::EventSetup const& iSetup){
  // If there are electrons, then the isolation maps and the SC should be in the event; if not it is an error

  string EGerrMsg("");
  bool foundCand=true;bool foundHandle=true;bool foundHcalIMap=true;bool foundTckIMap=true;
  edm::Handle<reco::ElectronCollection> electronNonIsoHandle;
  edm::Handle<reco::RecoEcalCandidateCollection> recoNonIsolecalcands;
  edm::Handle<reco::RecoEcalCandidateIsolationMap> HcalEleNonIsolMap;
  edm::Handle<reco::ElectronIsolationMap> TrackEleNonIsolMap;
  try {iEvent.getByLabel(CandNonIso_,recoNonIsolecalcands);} catch (...) {
    EGerrMsg=EGerrMsg + "  HLTEgamma: No non-isol eg candidate";
    foundCand=false;
  }
  try {iEvent.getByLabel(NonIsoElectronTag_,electronNonIsoHandle);} catch (...) { 
    EGerrMsg=EGerrMsg + "  HLTEgamma: No non-isol electron";
    foundHandle=false;
  }
  try {iEvent.getByLabel(NonIsoEleHcalTag_,HcalEleNonIsolMap);} catch (...) { 
    EGerrMsg=EGerrMsg + "  HLTEgamma: No non-isol Hcal electron";
    foundHcalIMap=false;
  }
  try {iEvent.getByLabel(NonIsoEleTrackIsolTag_,TrackEleNonIsolMap);} catch (...) { 
    EGerrMsg=EGerrMsg + "  HLTEgamma: No non-isol Track electron";
    foundTckIMap=false;
  }

  if (foundHandle){
    for(reco::ElectronCollection::const_iterator iElectron = electronNonIsoHandle->begin(); iElectron != 
	electronNonIsoHandle->end();iElectron++){

      // 1) find the SC from the electron
      reco::ElectronRef electronref(reco::ElectronRef(electronNonIsoHandle,iElectron - electronNonIsoHandle->begin()));
      const reco::SuperClusterRef theClus = electronref->superCluster(); //SC from the electron;

      //2) find the RecoEcalCandidate (wrapper for SC) corresponding to the SC of the electron
      reco::RecoEcalCandidateRef ref;
      if (foundCand){
	for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand= recoNonIsolecalcands->begin();
	     recoecalcand!=recoNonIsolecalcands->end(); recoecalcand++) {
	  ref = reco::RecoEcalCandidateRef(recoNonIsolecalcands,distance(recoNonIsolecalcands->begin(),recoecalcand));
	  reco::SuperClusterRef recrSC = ref->superCluster();
	  if(&(*recrSC) ==  &(*theClus)) {// ref is the RecoEcalCandidateRef correspoding to the electron
	    break;
	  }
	}
	//now ref is the RecoEcalCandidateRef correspoding to the electron

	myHLTElectron ele;
	ele.hcalIsol=-999; ele.trackIsol=-999;
	ele.L1Isolated = false;
	
	ele.Et=ref->et();
	ele.eta=ref->eta();
	ele.phi=ref->phi();
	ele.E = theClus->energy();
	ele.p=electronref->track()->momentum().R();

	//Fill the hcal Isolation
	if(foundHcalIMap){
	  reco::RecoEcalCandidateIsolationMap::const_iterator mapi = (*HcalEleNonIsolMap).find( ref);
	  if(mapi !=(*HcalEleNonIsolMap).end()) { ele.hcalIsol=mapi->val;}
	}
	//Fill the track Isolation
	if(foundTckIMap){
	  reco::ElectronIsolationMap::const_iterator mapTr = (*TrackEleNonIsolMap).find( electronref);
	  if(mapTr !=(*TrackEleNonIsolMap).end()) { ele.trackIsol=mapTr->val;}
	}

	//store the electron into the vector
	theHLTElectrons.push_back(ele);
      }

    }
  }

}

