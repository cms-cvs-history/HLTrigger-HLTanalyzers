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
    std::vector<std::string> parameterNames = myJetParams.getParameterNames() ;
    
    for ( std::vector<std::string>::iterator iParam = parameterNames.begin();
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
    jcale = new float[kMaxJetCal];
    jcalemf = new float[kMaxJetCal]; 
    jcaln90 = new float[kMaxJetCal]; 
    
    jcorcalpt = new float[kMaxJetCal]; 
    jcorcalphi = new float[kMaxJetCal]; 
    jcorcaleta = new float[kMaxJetCal]; 
    jcorcale = new float[kMaxJetCal]; 
    jcorcalemf = new float[kMaxJetCal]; 
    jcorcaln90 = new float[kMaxJetCal]; 
    
    const int kMaxJetgen = 10000;
    jgenpt = new float[kMaxJetgen];
    jgenphi = new float[kMaxJetgen];
    jgeneta = new float[kMaxJetgen];
    jgene = new float[kMaxJetgen];
    const int kMaxTower = 10000;
    towet = new float[kMaxTower];
    toweta = new float[kMaxTower];
    towphi = new float[kMaxTower];
    towen = new float[kMaxTower];
    towem = new float[kMaxTower];
    towhd = new float[kMaxTower];
    towoe = new float[kMaxTower];
    const int kMaxTau = 500;
    l2tauemiso = new float[kMaxTau];
    l25tauPt = new float[kMaxTau];
    l3tautckiso = new int[kMaxTau];
    tauEta = new float[kMaxTau];
    tauPt = new float[kMaxTau];
    tauPhi = new float[kMaxTau];
    
    const int kMaxPFTau = 500;
    pfTauEta         =  new float[kMaxPFTau];
    pfTauPhi         =  new float[kMaxPFTau];
    pfTauPt          =  new float[kMaxPFTau];
    pfTauJetPt       =  new float[kMaxPFTau];
    pfTauLeadTrackPt =  new float[kMaxPFTau];
    pfTauLeadPionPt  =  new float[kMaxPFTau];
    pfTauTrkIso         =  new int[kMaxPFTau];
    pfTauGammaIso         =  new int[kMaxPFTau];
    pfMHT   = -100;
    
    
    // Jet- MEt-specific branches of the tree 
    HltTree->Branch("NrecoJetCal",&njetcal,"NrecoJetCal/I");
    HltTree->Branch("NrecoJetGen",&njetgen,"NrecoJetGen/I");
    HltTree->Branch("NrecoTowCal",&ntowcal,"NrecoTowCal/I");
    HltTree->Branch("recoJetCalPt",jcalpt,"recoJetCalPt[NrecoJetCal]/F");
    HltTree->Branch("recoJetCalPhi",jcalphi,"recoJetCalPhi[NrecoJetCal]/F");
    HltTree->Branch("recoJetCalEta",jcaleta,"recoJetCalEta[NrecoJetCal]/F");
    HltTree->Branch("recoJetCalE",jcale,"recoJetCalE[NrecoJetCal]/F");
    HltTree->Branch("recoJetCalEMF",jcalemf,"recoJetCalEMF[NrecoJetCal]/F");
    HltTree->Branch("recoJetCalN90",jcaln90,"recoJetCalN90[NrecoJetCal]/F");
    
    HltTree->Branch("recoJetGenPt",jgenpt,"recoJetGenPt[NrecoJetGen]/F");
    HltTree->Branch("recoJetGenPhi",jgenphi,"recoJetGenPhi[NrecoJetGen]/F");
    HltTree->Branch("recoJetGenEta",jgeneta,"recoJetGenEta[NrecoJetGen]/F");
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
    //for(int ieta=0;ieta<NETA;ieta++){std::cout << " ieta " << ieta << " eta min " << CaloTowerEtaBoundries[ieta] <<std::endl;}
    
    HltTree->Branch("NrecoJetCorCal",&ncorjetcal,"NrecoJetCorCal/I"); 
    HltTree->Branch("recoJetCorCalPt",jcorcalpt,"recoJetCorCalPt[NrecoJetCorCal]/F"); 
    HltTree->Branch("recoJetCorCalPhi",jcorcalphi,"recoJetCorCalPhi[NrecoJetCorCal]/F"); 
    HltTree->Branch("recoJetCorCalEta",jcorcaleta,"recoJetCorCalEta[NrecoJetCorCal]/F"); 
    HltTree->Branch("recoJetCorCalE",jcorcale,"recoJetCorCalE[NrecoJetCorCal]/F"); 
    HltTree->Branch("recoJetCorCalEMF",jcorcalemf,"recoJetCorCalEMF[NrecoJetCorCal]/F");
    HltTree->Branch("recoJetCorCalN90",jcorcaln90,"recoJetCorCalN90[NrecoJetCorCal]/F");
    
    // Taus
    HltTree->Branch("NohTau",&nohtau,"NohTau/I");
    HltTree->Branch("ohTauEta",tauEta,"ohTauEta[NohTau]/F");
    HltTree->Branch("ohTauPhi",tauPhi,"ohTauPhi[NohTau]/F");
    HltTree->Branch("ohTauPt",tauPt,"ohTauPt[NohTau]/F");
    HltTree->Branch("ohTauEiso",l2tauemiso,"ohTauEiso[NohTau]/F");
    HltTree->Branch("ohTauL25Tpt",l25tauPt,"ohTauL25Tpt[NohTau]/F");
    HltTree->Branch("ohTauL3Tiso",l3tautckiso,"ohTauL3Tiso[NohTau]/I");
    
    HltTree->Branch("NohPFTau",&nohPFTau,"NohPFTau/I");
    HltTree->Branch("pfTauPt",pfTauPt,"pfTauPt[NohPFTau]/F");
    HltTree->Branch("pfTauEta",pfTauEta,"pfTauEta[NohPFTau]/F");
    HltTree->Branch("pfTauPhi",pfTauPhi,"pfTauPhi[NohPFTau]/F");
    HltTree->Branch("pfTauLeadTrackPt",pfTauLeadTrackPt,"pfTauLeadTrackPt[NohPFTau]/F");
    HltTree->Branch("pfTauLeadPionPt",pfTauLeadPionPt,"pfTauLeadPionPt[NohPFTau]/F");
    HltTree->Branch("pfTauTrkIso",pfTauTrkIso,"pfTauTrkIso[NohPFTau]/I");
    HltTree->Branch("pfTauGammaIso",pfTauGammaIso,"pfTauGammaIso[NohPFTau]/I");
    HltTree->Branch("pfTauJetPt",pfTauJetPt,"pfTauJetPt[NohPFTau]/F");    
    HltTree->Branch("pfMHT",&pfMHT,"pfMHT/F");
    
    
}

/* **Analyze the event** */
void HLTJets::analyze(const edm::Handle<reco::CaloJetCollection>      & calojets,
                      const edm::Handle<reco::CaloJetCollection>      & calocorjets,
                      const edm::Handle<reco::GenJetCollection>       & genjets,
                      const edm::Handle<reco::CaloMETCollection>      & recmets,
                      const edm::Handle<reco::GenMETCollection>       & genmets,
                      const edm::Handle<reco::METCollection>          & ht,
                      const edm::Handle<reco::HLTTauCollection>       & taujets,
                      const edm::Handle<reco::PFTauCollection>        & pfTaus,
                      const edm::Handle<CaloTowerCollection>          & caloTowers,
                      double thresholdForSavingTowers, 
                       double                minPtCH,
                       double               minPtGamma,
                      TTree * HltTree) {
    
    if (_Debug) std::cout << " Beginning HLTJets " << std::endl;
    
    //initialize branch variables
    njetcal=0; ncorjetcal=0; njetgen=0;ntowcal=0;
    mcalmet=0.; mcalphi=0.;
    mgenmet=0.; mgenphi=0.;
    htcalet=0.,htcalphi=0.,htcalsum=0.;
    
    if (calojets.isValid()) {
        reco::CaloJetCollection mycalojets;
        mycalojets=*calojets;
        std::sort(mycalojets.begin(),mycalojets.end(),PtGreater());
        typedef reco::CaloJetCollection::const_iterator cjiter;
        int jcal=0;
        for ( cjiter i=mycalojets.begin(); i!=mycalojets.end(); i++) {
            
            if (i->pt()>_CalJetMin){
                jcalpt[jcal] = i->pt();
                jcalphi[jcal] = i->phi();
                jcaleta[jcal] = i->eta();
                jcale[jcal] = i->energy();
                jcalemf[jcal] = i->emEnergyFraction();
                jcaln90[jcal] = i->n90();
                jcal++;
            }
            
        }
        njetcal = jcal;
    }
    else {njetcal = 0;}
    
    if (calocorjets.isValid()) {
        reco::CaloJetCollection mycalocorjets;
        mycalocorjets=*calocorjets;
        std::sort(mycalocorjets.begin(),mycalocorjets.end(),PtGreater());
        typedef reco::CaloJetCollection::const_iterator ccorjiter;
        int jcorcal=0;
        for ( ccorjiter i=mycalocorjets.begin(); i!=mycalocorjets.end(); i++) {
            
            if (i->pt()>_CalJetMin){
                jcorcalpt[jcorcal] = i->pt();
                jcorcalphi[jcorcal] = i->phi();
                jcorcaleta[jcorcal] = i->eta();
                jcorcale[jcorcal] = i->energy();
                jcorcalemf[jcorcal] = i->emEnergyFraction();
                jcorcaln90[jcorcal] = i->n90();
                jcorcal++;
            }
            
        }
        ncorjetcal = jcorcal;
    }
    else {ncorjetcal = 0;}
    
    if (caloTowers.isValid()) {
        //    ntowcal = caloTowers->size();
        int jtow = 0;
        for ( CaloTowerCollection::const_iterator tower=caloTowers->begin(); tower!=caloTowers->end(); tower++) {
            if(tower->energy() > thresholdForSavingTowers)
            {
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
        ntowcal = jtow;
    }
    else {ntowcal = 0;}
    
    if (recmets.isValid()) {
        typedef reco::CaloMETCollection::const_iterator cmiter;
        for ( cmiter i=recmets->begin(); i!=recmets->end(); i++) {
            mcalmet = i->pt();
            mcalphi = i->phi();
            mcalsum = i->sumEt();
        }
    }
    
    if (ht.isValid()) {
        typedef reco::METCollection::const_iterator iter;
        for ( iter i=ht->begin(); i!=ht->end(); i++) {
            htcalet = i->pt();
            htcalphi = i->phi();
            htcalsum = i->sumEt();
        }
    }
    
    if (_Monte){
        
        if (genjets.isValid()) {
            reco::GenJetCollection mygenjets;
            mygenjets=*genjets;
            std::sort(mygenjets.begin(),mygenjets.end(),PtGreater());
            typedef reco::GenJetCollection::const_iterator gjiter;
            int jgen=0;
            for ( gjiter i=mygenjets.begin(); i!=mygenjets.end(); i++) {
                
                if (i->pt()>_GenJetMin){
                    jgenpt[jgen] = i->pt();
                    jgenphi[jgen] = i->phi();
                    jgeneta[jgen] = i->eta();
                    jgene[jgen] = i->energy();
                    jgen++;
                }
                
            }
            njetgen = jgen;
        }
        else {njetgen = 0;}
        
        if (genmets.isValid()) {
            typedef reco::GenMETCollection::const_iterator gmiter;
            for ( gmiter i=genmets->begin(); i!=genmets->end(); i++) {
                mgenmet = i->pt();
                mgenphi = i->phi();
                mgensum = i->sumEt();
            }
        }
        
    }
    
    
    /////////////////////////////// Open-HLT Taus ///////////////////////////////
    
    if (taujets.isValid()) {      
        nohtau = taujets->size();
        reco::HLTTauCollection mytaujets;
        mytaujets=*taujets;
        std::sort(mytaujets.begin(),mytaujets.end(),GetPtGreater());
        typedef reco::HLTTauCollection::const_iterator tauit;
        int itau=0;
        for(tauit i=mytaujets.begin(); i!=mytaujets.end(); i++){
            //Ask for Eta,Phi and Et of the tau:
            tauEta[itau] = i->getEta();
            tauPhi[itau] = i->getPhi();
            tauPt[itau] = i->getPt();
            //Ask for L2 EMIsolation cut: Nominal cut : < 5
            l2tauemiso[itau] = i->getEMIsolationValue();
            //Get L25 LeadTrackPt : Nominal cut : > 20 GeV
            l25tauPt[itau] = i->getL25LeadTrackPtValue();
            //Get TrackIsolation response (returns 0 = failed or 1= passed)
            l3tautckiso[itau] = i->getL3TrackIsolationResponse();
            //MET : > 65
            itau++;
        }      
    }
    else {nohtau = 0;}

    
    ////////////////Particle Flow Taus ////////////////////////////////////
    if(pfTaus.isValid()) {
        float minTrkPt = minPtCH;
        float minGammaPt = minPtGamma;
        nohPFTau  = pfTaus->size();
        reco::PFTauCollection taus = *pfTaus;
        std::sort(taus.begin(),taus.end(),GetPFPtGreater());
        typedef reco::PFTauCollection::const_iterator pftauit;
        int ipftau=0;
        float pfMHTx = 0;
        float pfMHTy = 0;
        for(pftauit i=taus.begin(); i!=taus.end(); i++){
            //Ask for Eta,Phi and Et of the tau:
            pfTauEta[ipftau] = i->eta();
            pfTauPhi[ipftau] = i->phi();
            pfTauPt[ipftau] = i->pt();
            pfTauJetPt[ipftau] = i->pfTauTagInfoRef()->pfjetRef()->pt();

            pfMHTx = pfMHTx + i->pfTauTagInfoRef()->pfjetRef()->px();
            pfMHTy = pfMHTy + i->pfTauTagInfoRef()->pfjetRef()->py();
            
  /*
            if( (i->leadPFCand()).isNonnull())
                pfTauLeadPionPt[ipftau] = i->leadPFCand()->pt();            
*/
            if( (i->leadPFNeutralCand()).isNonnull())
                pfTauLeadPionPt[ipftau] = i->leadPFNeutralCand()->pt();        
            if((i->leadPFChargedHadrCand()).isNonnull())
                pfTauLeadTrackPt[ipftau] = i->leadPFChargedHadrCand()->pt();
            int myTrks=0;
            for (unsigned int iTrk = 0; iTrk < i->isolationPFChargedHadrCands().size(); iTrk++)
            {
                if(i->isolationPFChargedHadrCands()[iTrk]->pt() > minTrkPt) myTrks++;
            }
                
            pfTauTrkIso[ipftau] = myTrks;
            int myGammas=0;
            for (unsigned int iGamma = 0; iGamma < i->isolationPFGammaCands().size(); iGamma++)
            {
                if(i->isolationPFGammaCands()[iGamma]->pt() > minGammaPt) myGammas++;
            }                        
            pfTauGammaIso[ipftau] = myGammas;
        } 
        pfMHT = sqrt(pfMHTx*pfMHTx + pfMHTy*pfMHTy);
        
    }
    
 
        
}
