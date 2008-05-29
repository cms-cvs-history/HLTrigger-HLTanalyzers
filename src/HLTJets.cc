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

#include "HLTrigger/HLTanalyzers/interface/HLTJets.h"

HLTJets::HLTJets() {
  evtCounter=0;

  //set parameter defaults 
  _Monte=false;
  _Debug=false;
  _CalJetMin=0.;
  _GenJetMin=0.;
}

/*  Setup the analysis to put the branch-variables into the tree. */
void HLTJets::setup(const edm::ParameterSet& pSet, TTree* HltTree) {

  edm::ParameterSet myJetParams = pSet.getParameter<edm::ParameterSet>("RunParameters") ;
  vector<std::string> parameterNames = myJetParams.getParameterNames() ;
  
  for ( vector<std::string>::iterator iParam = parameterNames.begin();
	iParam != parameterNames.end(); iParam++ ){
    if  ( (*iParam) == "Monte" ) _Monte =  myJetParams.getParameter<bool>( *iParam );
    else if ( (*iParam) == "Debug" ) _Debug =  myJetParams.getParameter<bool>( *iParam );
    else if ( (*iParam) == "CalJetMin" ) _CalJetMin =  myJetParams.getParameter<double>( *iParam );
    else if ( (*iParam) == "GenJetMin" ) _GenJetMin =  myJetParams.getParameter<double>( *iParam );
  }

  const int kMaxJetCal = 10000;
  jcalpt = new float[kMaxJetCal];
  jcalphi = new float[kMaxJetCal];
  jcaleta = new float[kMaxJetCal];
  jcalet = new float[kMaxJetCal];
  jcale = new float[kMaxJetCal];
  const int kMaxJetgen = 10000;
  jgenpt = new float[kMaxJetgen];
  jgenphi = new float[kMaxJetgen];
  jgeneta = new float[kMaxJetgen];
  jgenet = new float[kMaxJetgen];
  jgene = new float[kMaxJetgen];
  const int kMaxTower = 10000;
  towet = new float[kMaxTower];
  toweta = new float[kMaxTower];
  towphi = new float[kMaxTower];
  towen = new float[kMaxTower];
  towem = new float[kMaxTower];
  towhd = new float[kMaxTower];
  towoe = new float[kMaxTower];
  const int kMaxTau1 = 500;
  l2t1emiso = new float[kMaxTau1];
  l25t1Pt = new float[kMaxTau1];
  l25t1tckiso = new int[kMaxTau1];
  l3t1Pt = new float[kMaxTau1];
  l3t1tckiso = new int[kMaxTau1];
  t1Eta = new float[kMaxTau1];
  t1Pt = new float[kMaxTau1];
  t1Phi = new float[kMaxTau1];

  // Jet- MEt-specific branches of the tree 
  HltTree->Branch("NrecoJetCal",&njetcal,"NrecoJetCal/I");
  HltTree->Branch("NrecoJetGen",&njetgen,"NrecoJetGen/I");
  HltTree->Branch("NrecoTowCal",&ntowcal,"NrecoTowCal/I");
  HltTree->Branch("recoJetCalPt",jcalpt,"recoJetCalPt[NrecoJetCal]/F");
  HltTree->Branch("recoJetCalPhi",jcalphi,"recoJetCalPhi[NrecoJetCal]/F");
  HltTree->Branch("recoJetCalEta",jcaleta,"recoJetCalEta[NrecoJetCal]/F");
  HltTree->Branch("recoJetCalEt",jcalet,"recoJetCalEt[NrecoJetCal]/F");
  HltTree->Branch("recoJetCalE",jcale,"recoJetCalE[NrecoJetCal]/F");
  HltTree->Branch("recoJetGenPt",jgenpt,"recoJetGenPt[NrecoJetGen]/F");
  HltTree->Branch("recoJetGenPhi",jgenphi,"recoJetGenPhi[NrecoJetGen]/F");
  HltTree->Branch("recoJetGenEta",jgeneta,"recoJetGenEta[NrecoJetGen]/F");
  HltTree->Branch("recoJetGenEt",jgenet,"recoJetGenEt[NrecoJetGen]/F");
  HltTree->Branch("recoJetGenE",jgene,"recoJetGenE[NrecoJetGen]/F");
  HltTree->Branch("recoTowEt",towet,"recoTowEt[NrecoTowCal]/F");
  HltTree->Branch("recoTowEta",toweta,"recoTowEta[NrecoTowCal]/F");
  HltTree->Branch("recoTowPhi",towphi,"recoTowPhi[NrecoTowCal]/F");
  HltTree->Branch("recoTowE",towen,"recoTowE[NrecoTowCal]/F");
  HltTree->Branch("recoTowEm",towem,"recoTowEm[NrecoTowCal]/F");
  HltTree->Branch("recoTowHad",towhd,"recoTowHad[NrecoTowCal]/F");
  HltTree->Branch("recoTowOE",towoe,"recoTowOE[NrecoTowCal]/F");
  HltTree->Branch("recoMetCal",&mcalmet,"recoMetCal/F");
  HltTree->Branch("recoMetCalPhi",&mcalphi,"recoMetCalPhi/F");
  HltTree->Branch("recoMetCalSum",&mcalsum,"recoMetCalSum/F");
  HltTree->Branch("recoMetGen",&mgenmet,"recoMetGen/F");
  HltTree->Branch("recoMetGenPhi",&mgenphi,"recoMetGenPhi/F");
  HltTree->Branch("recoMetGenSum",&mgensum,"recoMetGenSum/F");
  HltTree->Branch("recoHTCal",&htcalet,"recoHTCal/F");
  HltTree->Branch("recoHTCalPhi",&htcalphi,"recoHTCalPhi/F");
  HltTree->Branch("recoHTCalSum",&htcalsum,"recoHTCalSum/F");

  // Taus
  HltTree->Branch("NohTau1",&nohtau1,"NohTau1/I");
  HltTree->Branch("ohTau1Eta",t1Eta,"ohTau1Eta[NohTau1]/F");
  HltTree->Branch("ohTau1Phi",t1Phi,"ohTau1Phi[NohTau1]/F");
  HltTree->Branch("ohTau1Pt",t1Pt,"ohTau1Pt[NohTau1]/F");
  HltTree->Branch("ohTau1Eiso",l2t1emiso,"ohTau1Eiso[NohTau1]/F");
  HltTree->Branch("ohTau1L25Tpt",l25t1Pt,"ohTau1L25Tpt[NohTau1]/F");
  HltTree->Branch("ohTau1L25Tiso",l25t1tckiso,"ohTau1L25Tiso[NohTau1]/I");
  HltTree->Branch("ohTau1L3Tpt",l3t1Pt,"ohTau1L3Tpt[NohTau1]/F");
  HltTree->Branch("ohTau1L3Tiso",l3t1tckiso,"ohTau1L3Tiso[NohTau1]/I");

  //for(int ieta=0;ieta<NETA;ieta++){cout << " ieta " << ieta << " eta min " << CaloTowerEtaBoundries[ieta] <<endl;}

}

/* **Analyze the event** */
void HLTJets::analyze(const CaloJetCollection& calojets,
		      const GenJetCollection& genjets,
		      const CaloMETCollection& recmets,
		      const GenMETCollection& genmets,
		      const METCollection& ht,
		      const reco::HLTTauCollection& myHLT1Tau,
// 		      const reco::HLTTauCollection& myHLT2Tau,
		      const CaloTowerCollection& caloTowers,
		      TTree* HltTree) {

  if (_Debug) std::cout << " Beginning HLTJets " << std::endl;

  //initialize branch variables
  njetcal=0; njetgen=0;ntowcal=0;
  mcalmet=0.; mcalphi=0.;
  mgenmet=0.; mgenphi=0.;
  htcalet=0.,htcalphi=0.,htcalsum=0.;

  if (&calojets) {
    CaloJetCollection mycalojets;
    mycalojets=calojets;
    std::sort(mycalojets.begin(),mycalojets.end(),PtGreater());
//     njetcal = mycalojets.size();
    typedef CaloJetCollection::const_iterator cjiter;
    int jcal=0;
    for ( cjiter i=mycalojets.begin(); i!=mycalojets.end(); i++) {

      if (i->pt()>_CalJetMin){
	jcalpt[jcal] = i->pt();
	jcalphi[jcal] = i->phi();
	jcaleta[jcal] = i->eta();
	jcalet[jcal] = i->et();
	jcale[jcal] = i->energy();
	jcal++;
      }

    }
    njetcal = jcal;
  }
  else {njetcal = 0;}

  if (&caloTowers){
    ntowcal = caloTowers.size();
    int jtow = 0;
    for ( CaloTowerCollection::const_iterator tower=caloTowers.begin(); tower!=caloTowers.end(); tower++) {
      towet[jtow] = tower->et();
      toweta[jtow] = tower->eta();
      towphi[jtow] = tower->phi();
      towen[jtow] = tower->energy();
      towem[jtow] = tower->emEnergy();
      towhd[jtow] = tower->hadEnergy();
      towoe[jtow] = tower->outerEnergy();
      jtow++;
    }
  }
  else {ntowcal = 0;}

  if (&recmets) {
    typedef CaloMETCollection::const_iterator cmiter;
    for ( cmiter i=recmets.begin(); i!=recmets.end(); i++) {
      mcalmet = i->pt();
      mcalphi = i->phi();
      mcalsum = i->sumEt();
    }
  }

  if (&ht) {
    typedef METCollection::const_iterator iter;
    for ( iter i=ht.begin(); i!=ht.end(); i++) {
      htcalet = i->pt();
      htcalphi = i->phi();
      htcalsum = i->sumEt();
    }
  }

  if (_Monte){

    if (&genjets) {
      GenJetCollection mygenjets;
      mygenjets=genjets;
      std::sort(mygenjets.begin(),mygenjets.end(),PtGreater());
//       njetgen = mygenjets.size();
      typedef GenJetCollection::const_iterator gjiter;
      int jgen=0;
      for ( gjiter i=mygenjets.begin(); i!=mygenjets.end(); i++) {

	if (i->pt()>_GenJetMin){
	  jgenpt[jgen] = i->pt();
	  jgenphi[jgen] = i->phi();
	  jgeneta[jgen] = i->eta();
	  jgenet[jgen] = i->et();
	  jgene[jgen] = i->energy();
	  jgen++;
	}

      }
       njetgen = jgen;
    }
    else {njetgen = 0;}

    if (&genmets) {
      typedef GenMETCollection::const_iterator gmiter;
      for ( gmiter i=genmets.begin(); i!=genmets.end(); i++) {
	mgenmet = i->pt();
	mgenphi = i->phi();
	mgensum = i->sumEt();
      }
    }

  }


/////////////////////////////// Open-HLT Taus ///////////////////////////////


    if (&myHLT1Tau){
      
      nohtau1 = myHLT1Tau.size();
      HLTTauCollection mytau1jets;
      mytau1jets=myHLT1Tau;
      std::sort(mytau1jets.begin(),mytau1jets.end(),GetPtGreater());
      typedef HLTTauCollection::const_iterator t1it;
      int it1=0;
      for(t1it i=mytau1jets.begin(); i!=mytau1jets.end(); i++){
	//Ask for Eta,Phi and Et of the tau:
	t1Eta[it1] = i->getEta();
	t1Phi[it1] = i->getPhi();
	t1Pt[it1] = i->getPt();
	//Ask for L2 EMIsolation cut: Nominal cut : < 5
	l2t1emiso[it1] = i->getEMIsolationValue();
	//Get L25 LeadTrackPt : Nominal cut : > 20 GeV
	l25t1Pt[it1] = i->getL25LeadTrackPtValue();
	//Get TrackIsolation response (returns 0 = failed or 1= passed)
	l25t1tckiso[it1] = i->getL25TrackIsolationResponse();
	//Get L3 LeadTrackPt : Nominal cut : > 20 GeV
	l3t1Pt[it1] = i->getL3LeadTrackPtValue();
	//Get TrackIsolation response (returns 0 = failed or 1= passed)
	l3t1tckiso[it1] = i->getL3TrackIsolationResponse();
	//MET : > 65
	it1++;
      }
      
    }
    else {nohtau1 = 0;}

}
