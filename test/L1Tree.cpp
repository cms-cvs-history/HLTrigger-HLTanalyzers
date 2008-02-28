#define L1Tree_cxx
#include "L1Tree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <TFormula.h>

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

void L1Tree::Loop(vector<int> * iCount
		  ,vector<int> * sPureCount, vector<int> * pureCount
		  ,vector< vector<int> > * overlapCount
		  ,vector<string> trignames,vector<int> prescales,int NEntries
		  ,bool doMuonCut,bool doElecCut) {
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = (Long64_t)NEntries; 
   if (NEntries <= 0)
     nentries = fChain->GetEntries();

   vector<int> iCountNoPrescale;
   for (int it = 0; it < Ntrig; it++){
     iCountNoPrescale.push_back(0);
   }
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%100000 == 0) cout<<"Processing entry "<<jentry<<"/"<<nentries<<"\r"<<flush<<endl;

      // Cut on muon quality
      // init
      for (int i=0;i<10;i++) {
	NL1GoodSingleMu = 0;
	L1GoodSingleMuPt[i] = -999.;
	L1GoodSingleMuE[i] = -999.;
	L1GoodSingleMuEta[i] = -999.;
	L1GoodSingleMuPhi[i] = -999.;
	L1GoodSingleMuIsol[i] = -999;
	L1GoodSingleMuMip[i] = -999;
	L1GoodSingleMuFor[i] = -999;
	L1GoodSingleMuRPC[i] = -999;
	L1GoodSingleMuQal[i] = -999;     
      }
      // Cut on muon quality
      for (int i=0;i<NL1Mu;i++) {
	if ( L1MuQal[i]==4 || L1MuQal[i]==5 || L1MuQal[i]==6 || L1MuQal[i]==7 ) {
	  L1GoodSingleMuPt[NL1GoodSingleMu] = L1MuPt[i];
	  L1GoodSingleMuE[NL1GoodSingleMu] = L1MuE[i];
	  L1GoodSingleMuEta[NL1GoodSingleMu] = L1MuEta[i];
	  L1GoodSingleMuPhi[NL1GoodSingleMu] = L1MuPhi[i];
	  L1GoodSingleMuIsol[NL1GoodSingleMu] = L1MuIsol[i];
	  L1GoodSingleMuMip[NL1GoodSingleMu] = L1MuMip[i];
	  L1GoodSingleMuFor[NL1GoodSingleMu] = L1MuFor[i];
	  L1GoodSingleMuRPC[NL1GoodSingleMu] = L1MuRPC[i];
	  L1GoodSingleMuQal[NL1GoodSingleMu] = L1MuQal[i];
	  NL1GoodSingleMu++;
	}
      }

      // init
      for (int i=0;i<10;i++) {
	NL1GoodDoubleMu = 0;
	L1GoodDoubleMuPt[i] = -999.;
	L1GoodDoubleMuE[i] = -999.;
	L1GoodDoubleMuEta[i] = -999.;
	L1GoodDoubleMuPhi[i] = -999.;
	L1GoodDoubleMuIsol[i] = -999;
	L1GoodDoubleMuMip[i] = -999;
	L1GoodDoubleMuFor[i] = -999;
	L1GoodDoubleMuRPC[i] = -999;
	L1GoodDoubleMuQal[i] = -999;     
      }
      // Cut on muon quality
      for (int i=0;i<NL1Mu;i++) {
	if ( L1MuQal[i]==3 || L1MuQal[i]==5 || L1MuQal[i]==6 || L1MuQal[i]==7 ) {
	  L1GoodDoubleMuPt[NL1GoodDoubleMu] = L1MuPt[i];
	  L1GoodDoubleMuE[NL1GoodDoubleMu] = L1MuE[i];
	  L1GoodDoubleMuEta[NL1GoodDoubleMu] = L1MuEta[i];
	  L1GoodDoubleMuPhi[NL1GoodDoubleMu] = L1MuPhi[i];
	  L1GoodDoubleMuIsol[NL1GoodDoubleMu] = L1MuIsol[i];
	  L1GoodDoubleMuMip[NL1GoodDoubleMu] = L1MuMip[i];
	  L1GoodDoubleMuFor[NL1GoodDoubleMu] = L1MuFor[i];
	  L1GoodDoubleMuRPC[NL1GoodDoubleMu] = L1MuRPC[i];
	  L1GoodDoubleMuQal[NL1GoodDoubleMu] = L1MuQal[i];
	  NL1GoodDoubleMu++;
	}
      }
      
      /*
      for (int i=0;i<NMCpart;i++) {
	if (MCpid[i]==11 || MCpid[i]==-11 || MCpid[i]==22) {
	  if (MCpt[i]>5.) {
	    //cout<<jentry<<": "<<MCpt[i]<<endl;
	    NMCmu3Andel1Events++;
	    break;
	  }
	}
      }
      for (int i=0;i<NMCpart;i++) {
	if (MCpid[i]==11 || MCpid[i]==-11) {
	  if (MCpt[i]>5.) {
	    //cout<<jentry<<": "<<MCpt[i]<<endl;
	    NMCel1Events++;
	    break;
	  }
	}
      }
      for (int i=0;i<NMCpart;i++) {
	if (MCpid[i]==13 || MCpid[i]==-13) {
	  if (MCpt[i]>5.) {
	    //cout<<jentry<<": "<<MCpt[i]<<endl;
	    NMCmu3Events++;
	    break;
	  }
	}
      }
      */
      //if (MCmu3>0) cout<<"***** MCMu3: "<<jentry<<endl;
      //if (MCel1>0) cout<<"----- MCEl1: "<<jentry<<endl;
      //if (MCel1>0 && MCmu3>0) cout<<"===== MCEl1 & MCMu3: "<<jentry<<endl;
      /*
      if (MCmu3>0) NMCmu3Events++;
      if (MCel1>0) NMCel1Events++;
      if (MCel1>0 && MCmu3>0) NMCmu3Andel1Events++;
      if (jentry%10000==0) {
	cout<<setprecision(8)<<"NMCmu3Events: "<<NMCmu3Events<<" / "<<jentry<<" = "<<(double)NMCmu3Events/(double)jentry<<endl;
	cout<<setprecision(8)<<"NMCel1Events: "<<NMCel1Events<<" / "<<jentry<<" = "<<(double)NMCel1Events/(double)jentry<<endl;
	cout<<setprecision(8)<<"NMCmu3Andel1Events: "<<NMCmu3Andel1Events<<" / "<<jentry<<" = "<<(double)NMCmu3Andel1Events/(double)jentry<<endl;
      }
      
      */
      /*
      if (jentry%10000==0) {
	cout<<setprecision(8)<<"NMCmu5Events: "<<NMCmu3Events<<" / "<<jentry<<" = "<<(double)NMCmu3Events/(double)jentry<<endl;
	cout<<setprecision(8)<<"NMCel5Events: "<<NMCel1Events<<" / "<<jentry<<" = "<<(double)NMCel1Events/(double)jentry<<endl;
	cout<<setprecision(8)<<"NMCeg5Events: "<<NMCmu3Andel1Events<<" / "<<jentry<<" = "<<(double)NMCmu3Andel1Events/(double)jentry<<endl;
      }
      */
      
      // 1. Loop to check which Bit fired
      // Triggernames are assigned to trigger cuts in unambigous way!
      // If you define a new trigger also define a new unambigous name!
      for (int it = 0; it < Ntrig; it++){
	triggerBit[it] = false;
	triggerBitNoPrescale[it] = false;
	previousBitsFired[it] = false;
	allOtherBitsFired[it] = false;

	/* This is the trigger part*/
	if ( ((!doMuonCut) || (MCmu3==0))
	     && ((!doElecCut) || (MCel1==0)) ) { // Avoids muon/elec double counting

	  // To be extended!
	  if (trignames[it].compare("OrAllMu") == 0) {
	    if ( (L1GoodSingleMuPt[0]>=0.) ||
		 ( NL1GoodDoubleMu>=2 && L1GoodDoubleMuPt[1]>=3.) ) {
	      triggerBitNoPrescale[it] = true;
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }	  
	    }
	  }



	  
	  if (trignames[it].compare("SingleMu0") == 0) {
	    if ( (NL1GoodSingleMu>=1 && L1GoodSingleMuPt[0]>=0.) ) {
	      triggerBitNoPrescale[it] = true;
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }	  
	    }
	  }
	  else if (trignames[it].compare("SingleMu3") == 0) { 
	    if ( (NL1GoodSingleMu>=1 && L1GoodSingleMuPt[0]>=3.) ) {
	      triggerBitNoPrescale[it] = true;
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleMu5") == 0) { 
	    //if (jentry%10000==0) cout<<setprecision(8)<<"SingleMu5: "<<iCount->at(it)<<" / "<<jentry<<" = "<<(double)iCount->at(it)/(double)jentry<<endl;
	    if ( (NL1GoodSingleMu>=1 && L1GoodSingleMuPt[0]>=5.) ) {
	      triggerBitNoPrescale[it] = true;
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleMu7") == 0) { 
	    if ( (NL1GoodSingleMu>=1 && L1GoodSingleMuPt[0]>=7.) ) {
	      triggerBitNoPrescale[it] = true;
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleMu3") == 0) { 
	    if ( (L1_SingleMu3==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleMu5") == 0) { 
	    if ( (L1_SingleMu5==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleMu7") == 0) { 
	    //if (jentry%10000==0) cout<<setprecision(8)<<"L1_SingleMu7: "<<iCount->at(it)<<" / "<<jentry<<" = "<<(double)iCount->at(it)/(double)jentry<<endl;
	    if ( (L1_SingleMu7==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleMu10") == 0) { 
	    if ( (L1_SingleMu10==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleMu14") == 0) { 
	    if ( (L1_SingleMu14==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleMu20") == 0) { 
	    if ( (L1_SingleMu20==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleMu25") == 0) { 
	    if ( (L1_SingleMu25==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_DoubleMu3") == 0) { 
	    if ( (L1_DoubleMu3==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_TripleMu3") == 0) { 
	    if ( (L1_TripleMu3==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_Mu3_Jet15") == 0) { 
	    if ( (L1_Mu3_Jet15==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_Mu5_Jet15") == 0) { 
	    if ( (L1_Mu5_Jet15==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_Mu3_Jet70") == 0) { 
	    if ( (L1_Mu3_Jet70==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_Mu5_Jet20") == 0) { 
	    if ( (L1_Mu5_Jet20==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_Mu3_IsoEG5") == 0) { 
	    if ( (L1_Mu3_IsoEG5==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_Mu5_IsoEG10") == 0) { 
	    if ( (L1_Mu5_IsoEG10==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_Mu3_EG12") == 0) { 
	    if ( (L1_Mu3_EG12==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleEG5") == 0) { 
	    //if (jentry%10000==0) cout<<setprecision(8)<<"SingleEG5: "<<iCount->at(it)<<" / "<<jentry<<" = "<<(double)iCount->at(it)/(double)jentry<<endl;
	    if ( (L1IsolEmEt[0]>=5. || L1NIsolEmEt[0]>=5.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleEG6") == 0) { 
	    if ( (L1IsolEmEt[0]>=6. || L1NIsolEmEt[0]>=6.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleEG7") == 0) { 
	    if ( (L1IsolEmEt[0]>=7. || L1NIsolEmEt[0]>=7.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleEG8") == 0) { 
	    if ( (L1IsolEmEt[0]>=8. || L1NIsolEmEt[0]>=8.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleEG9") == 0) { 
	    if ( (L1IsolEmEt[0]>=9. || L1NIsolEmEt[0]>=9.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleEG10") == 0) { 
	    if ( (L1IsolEmEt[0]>=10. || L1NIsolEmEt[0]>=10.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleEG12") == 0) { 
	    if ( (L1IsolEmEt[0]>=12. || L1NIsolEmEt[0]>=12.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleEG15") == 0) { 
	    if ( (L1IsolEmEt[0]>=15. || L1NIsolEmEt[0]>=15.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleIsoEG12") == 0) { 
	    //if (jentry%10000==0) cout<<setprecision(8)<<"SingleIsoEG12: "<<iCount->at(it)<<" / "<<jentry<<" = "<<(double)iCount->at(it)/(double)jentry<<endl;
	    if ( (L1IsolEmEt[0]>=12.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleIsoEG15") == 0) { 
	    if ( (L1IsolEmEt[0]>=15.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleIsoEG8") == 0) { 
	    if ( (L1_SingleIsoEG8==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleIsoEG10") == 0) { 
	    if ( (L1_SingleIsoEG10==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleIsoEG12") == 0) { 
	    if ( (L1_SingleIsoEG12==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleIsoEG15") == 0) { 
	    //if (jentry%10000==0) cout<<setprecision(8)<<"L1_SingleIsoEG15: "<<iCount->at(it)<<" / "<<jentry<<" = "<<(double)iCount->at(it)/(double)jentry<<endl;
	    if ( (L1_SingleIsoEG15==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleIsoEG20") == 0) { 
	    if ( (L1_SingleIsoEG20==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleIsoEG25") == 0) { 
	    if ( (L1_SingleIsoEG25==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleEG5") == 0) { 
	    if ( (L1_SingleEG5==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleEG8") == 0) { 
	    if ( (L1_SingleEG8==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleEG10") == 0) { 
	    if ( (L1_SingleEG10==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleEG12") == 0) { 
	    if ( (L1_SingleEG12==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleEG15") == 0) { 
	    if ( (L1_SingleEG15==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleEG20") == 0) { 
	    if ( (L1_SingleEG20==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleEG25") == 0) { 
	    if ( (L1_SingleEG25==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleJet10") == 0) { 
	    if ( (L1TauEt[0]>=10. || L1CenJetEt[0]>=10. || L1ForJetEt[0]>=10.) ) {
	      triggerBitNoPrescale[it] = true;
	      if ( (iCountNoPrescale[it] % prescales[it] == 0) ) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleJet15") == 0) { 
	    if ( (L1TauEt[0]>=15. || L1CenJetEt[0]>=15. || L1ForJetEt[0]>=15.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleJet20") == 0) { 
	    if ( (L1TauEt[0]>=20. || L1CenJetEt[0]>=20. || L1ForJetEt[0]>=20.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleJet25") == 0) { 
	    if ( (L1TauEt[0]>=25. || L1CenJetEt[0]>=25. || L1ForJetEt[0]>=25.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleJet30") == 0) { 
	    if ( (L1TauEt[0]>=30. || L1CenJetEt[0]>=30. || L1ForJetEt[0]>=30.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleJet35") == 0) { 
	    if ( (L1TauEt[0]>=35. || L1CenJetEt[0]>=35. || L1ForJetEt[0]>=35.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleJet40") == 0) { 
	    if ( (L1TauEt[0]>=40. || L1CenJetEt[0]>=40. || L1ForJetEt[0]>=40.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleJet50") == 0) { 
	    if ( (L1TauEt[0]>=50. || L1CenJetEt[0]>=50. || L1ForJetEt[0]>=50.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleJet70") == 0) { 
	    if ( (L1TauEt[0]>=70. || L1CenJetEt[0]>=70. || L1ForJetEt[0]>=70.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleJet100") == 0) {
	    if ( (L1TauEt[0]>=100. || L1CenJetEt[0]>=100. || L1ForJetEt[0]>=100.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleJet150") == 0) { 
	    if ( (L1TauEt[0]>=150. || L1CenJetEt[0]>=150. || L1ForJetEt[0]>=150.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleJet15") == 0) { 
	    if ( (L1_SingleJet15==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleJet20") == 0) { 
	    if ( (L1_SingleJet20==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleJet30") == 0) { 
	    if ( (L1_SingleJet30==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleJet50") == 0) { 
	    if ( (L1_SingleJet50==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleJet70") == 0) { 
	    if ( (L1_SingleJet70==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleJet100") == 0) { 
	    if ( (L1_SingleJet100==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleJet150") == 0) { 
	    if ( (L1_SingleJet150==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleJet200") == 0) { 
	    if ( (L1_SingleJet200==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleTauJet10") == 0) { 
	    if ( (L1TauEt[0]>=10.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleTauJet20") == 0) { 
	    if ( (L1TauEt[0]>=20. ) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleTauJet25") == 0) { 
	    if ( (L1TauEt[0]>=25. ) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleTauJet30") == 0) { 
	    if ( (L1TauEt[0]>=30. ) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleTauJet35") == 0) { 
	    if ( (L1TauEt[0]>=35. ) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleTauJet40") == 0) { 
	    if ( (L1TauEt[0]>=40. ) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleTauJet50") == 0) { 
	    if ( (L1TauEt[0]>=50. ) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleTauJet60") == 0) { 
	    if ( (L1TauEt[0]>=60. ) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleTauJet70") == 0) { 
	    if ( (L1TauEt[0]>=70. ) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleTauJet80") == 0) { 
	    if ( (L1TauEt[0]>=80. ) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("SingleTauJet100") == 0) { 
	    if ( (L1TauEt[0]>=100. ) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleTauJet10") == 0) { 
	    if ( (L1_SingleTauJet10==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleTauJet20") == 0) { 
	    if ( (L1_SingleTauJet20==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleTauJet30") == 0) { 
	    if ( (L1_SingleTauJet30==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleTauJet40") == 0) { 
	    if ( (L1_SingleTauJet40==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleTauJet60") == 0) { 
	    if ( (L1_SingleTauJet60==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleTauJet80") == 0) { 
	    if ( (L1_SingleTauJet80==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_SingleTauJet100") == 0) { 
	    if ( (L1_SingleTauJet100==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_HTT100") == 0) { 
	    if ( (L1_HTT100==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_HTT200") == 0) { 
	    if ( (L1_HTT200==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_HTT250") == 0) { 
	    if ( (L1_HTT250==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_HTT300") == 0) { 
	    if ( (L1_HTT300==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_HTT400") == 0) { 
	    if ( (L1_HTT400==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_HTT500") == 0) { 
	    if ( (L1_HTT500==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM40_Jet30") == 0) { 
	    if ( (L1Met>=40.) && (L1TauEt[0]>=30. || L1CenJetEt[0]>=30. || L1ForJetEt[0]>=30.) ) {
	      triggerBitNoPrescale[it] = true;
	      if ( (iCountNoPrescale[it] % prescales[it] == 0) ) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM45_Jet30") == 0) { 
	    if ( (L1Met>=45.) && (L1TauEt[0]>=30. || L1CenJetEt[0]>=30. || L1ForJetEt[0]>=30.) ) {
	      triggerBitNoPrescale[it] = true;
	      if ( (iCountNoPrescale[it] % prescales[it] == 0) ) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM50_Jet40") == 0) { 
	    if ( (L1Met>=50.) && (L1TauEt[0]>=40. || L1CenJetEt[0]>=40. || L1ForJetEt[0]>=40.) ) {
	      triggerBitNoPrescale[it] = true;
	      if ( (iCountNoPrescale[it] % prescales[it] == 0) ) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM10") == 0) { 
	    if ( (L1Met>=10.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM15") == 0) { 
	    if ( (L1Met>=15.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM20") == 0) { 
	    if ( (L1Met>=20.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM25") == 0) { 
	    if ( (L1Met>=25.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM30") == 0) { 
	    if ( (L1Met>=30.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM35") == 0) { 
	    if ( (L1Met>=35.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM40") == 0) { 
	    if ( (L1Met>=40.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM45") == 0) { 
	    if ( (L1Met>=45.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM50") == 0) { 
	    if ( (L1Met>=50.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("ETM60") == 0) { 
	    if ( (L1Met>=60.) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_ETM20") == 0) { 
	    if ( (L1_ETM20==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_ETM30") == 0) { 
	    if ( (L1_ETM30==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_ETM40") == 0) { 
	    if ( (L1_ETM40==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_ETM50") == 0) { 
	    if ( (L1_ETM50==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_ETM60") == 0) { 
	    if ( (L1_ETM60==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_DoubleIsoEG8") == 0) { 
	    if ( (L1_DoubleIsoEG8==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_DoubleIsoEG10") == 0) { 
	    if ( (L1_DoubleIsoEG10==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_DoubleEG5") == 0) { 
	    if ( (L1_DoubleEG5==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_DoubleEG10") == 0) { 
	    if ( (L1_DoubleEG10==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_DoubleEG15") == 0) { 
	    if ( (L1_DoubleEG15==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_DoubleJet70") == 0) { 
	    if ( (L1_DoubleJet70==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_DoubleJet100") == 0) { 
	    if ( (L1_DoubleJet100==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_DoubleTauJet40") == 0) { 
	    if ( (L1_DoubleTauJet40==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_IsoEG10_Jet15") == 0) { 
	    if ( (L1_IsoEG10_Jet15==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_IsoEG10_Jet20") == 0) { 
	    if ( (L1_IsoEG10_Jet20==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_IsoEG10_Jet30") == 0) { 
	    if ( (L1_IsoEG10_Jet30==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_IsoEG10_Jet70") == 0) { 
	    if ( (L1_IsoEG10_Jet70==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_IsoEG10_TauJet20") == 0) { 
	    if ( (L1_IsoEG10_TauJet20==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_IsoEG10_TauJet30") == 0) { 
	    if ( (L1_IsoEG10_TauJet30==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_TauJet20_ETM20") == 0) { 
	    if ( (L1_TauJet20_ETM20==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_TauJet30_ETM30") == 0) { 
	    if ( (L1_TauJet30_ETM30==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_TauJet30_ETM40") == 0) { 
	    if ( (L1_TauJet30_ETM40==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_HTT100_ETM30") == 0) { 
	    if ( (L1_HTT100_ETM30==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_TripleJet50") == 0) { 
	    if ( (L1_TripleJet50==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_QuadJet30") == 0) { 
	    if ( (L1_QuadJet30==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("QuadJet30") == 0) { 
	    int nCount=0;
	    for (int i=0;i<4;i++) {
	      if (L1CenJetEt[i]>=30.) nCount++;
	      if (L1ForJetEt[i]>=30.) nCount++;
	      if (L1TauEt[i]>=30.) nCount++;
	      if (nCount>=4) break;
	    }
	    if ( nCount>=4 ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("QuadJet40") == 0) { 
	    int nCount=0;
	    for (int i=0;i<4;i++) {
	      if (L1CenJetEt[i]>=40.) nCount++;
	      if (L1ForJetEt[i]>=40.) nCount++;
	      if (L1TauEt[i]>=40.) nCount++;
	      if (nCount>=4) break;
	    }
	    if ( nCount>=4 ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("QuadJet50") == 0) { 
            int nCount=0;
            for (int i=0;i<4;i++) {
              if (L1CenJetEt[i]>=50.) nCount++;
              if (L1ForJetEt[i]>=50.) nCount++;
              if (L1TauEt[i]>=50.) nCount++;
              if (nCount>=4) break;
            }
            if ( nCount>=4 ) {
              triggerBitNoPrescale[it] = true;
              if (iCountNoPrescale[it] % prescales[it] == 0) {
                triggerBit[it] = true;
              }
            }
	  }
	  else if (trignames[it].compare("L1_ExclusiveDoubleIsoEG6") == 0) { 
	    if ( (L1_ExclusiveDoubleIsoEG6==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
	      triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_ExclusiveDoubleJet60") == 0) { 
	    if ( (L1_ExclusiveDoubleJet60==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_ExclusiveJet25_Gap_Jet25") == 0) { 
	    if ( (L1_ExclusiveJet25_Gap_Jet25==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_IsoEG10_Jet20_ForJet10") == 0) { 
	    if ( (L1_IsoEG10_Jet20_ForJet10==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }	  
	  }
	  else if (trignames[it].compare("L1_MinBias_HTT10") == 0) { 
	    if ( (L1_MinBias_HTT10==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("L1_ZeroBias") == 0) { 
	    if ( (L1_ZeroBias==1) ) { 
	      triggerBitNoPrescale[it] = true; 
	      if (iCountNoPrescale[it] % prescales[it] == 0) {
		triggerBit[it] = true;
	      }
	    }
	  }
	  else if (trignames[it].compare("MinBias_SingleHF1") == 0) { 
	    for (int i=0;i<NrecoTowCal;i++) {
	      if (recoTowEta[i]>3.0 || recoTowEta[i]<-3.0 ) {
		if (recoTowEt[i]>=1.0) {
		  triggerBitNoPrescale[it] = true;
		  if ( iCountNoPrescale[it] % prescales[it] == 0) {
		    triggerBit[it] = true;
		    break;
		  }
		}
	      }
	    }
	  }
	  else if (trignames[it].compare("MinBias_DoubleHF1") == 0) { 
	    bool forwSideFired = false;
	    bool backSideFired = false;
	    for (int i=0;i<NrecoTowCal;i++) {
	      if (recoTowEta[i]>3.0 && recoTowEt[i]>=1.0) forwSideFired = true;
	      if (recoTowEta[i]<-3.0 && recoTowEt[i]>=1.0) backSideFired = true;
	      if (forwSideFired && backSideFired) {
		triggerBitNoPrescale[it] = true;
		if ( iCountNoPrescale[it] % prescales[it] == 0 ) {
		  triggerBit[it] = true;
		  break;
		}
	      }
	    }
	  }
	  else if (trignames[it].compare("MinBias_SingleHF2") == 0) { 
	    for (int i=0;i<NrecoTowCal;i++) {
	      if (recoTowEta[i]>3.0 || recoTowEta[i]<-3.0 ) {
		if (recoTowEt[i]>=2.0) {
		  triggerBitNoPrescale[it] = true;
		  if ( iCountNoPrescale[it] % prescales[it] == 0) {
		    triggerBit[it] = true;
		    break;
		  }
		}
	      }
	    }
	  }
	  else if (trignames[it].compare("MinBias_DoubleHF2") == 0) { 
	    bool forwSideFired = false;
	    bool backSideFired = false;
	    for (int i=0;i<NrecoTowCal;i++) {
	      if (recoTowEta[i]>3.0 && recoTowEt[i]>=2.0) forwSideFired = true;
	      if (recoTowEta[i]<-3.0 && recoTowEt[i]>=2.0) backSideFired = true;
	      if (forwSideFired && backSideFired) {
		triggerBitNoPrescale[it] = true;
		if ( iCountNoPrescale[it] % prescales[it] == 0 ) {
		  triggerBit[it] = true;
		  break;
		}
	      }
	    }
	  }


	  else if (trignames[it].compare("MinBias_SingleHF3") == 0) { 
	    for (int i=0;i<NrecoTowCal;i++) {
	      if (recoTowEta[i]>3.0 || recoTowEta[i]<-3.0 ) {
		if (recoTowEt[i]>=3.0) {
		  triggerBitNoPrescale[it] = true;
		  if ( iCountNoPrescale[it] % prescales[it] == 0) {
		    triggerBit[it] = true;
		    break;
		  }
		}
	      }
	    }
	  }
	  else if (trignames[it].compare("MinBias_DoubleHF3") == 0) { 
	    bool forwSideFired = false;
	    bool backSideFired = false;
	    for (int i=0;i<NrecoTowCal;i++) {
	      if (recoTowEta[i]>3.0 && recoTowEt[i]>=3.0) forwSideFired = true;
	      if (recoTowEta[i]<-3.0 && recoTowEt[i]>=3.0) backSideFired = true;
	      if (forwSideFired && backSideFired) {
		triggerBitNoPrescale[it] = true;
		if ( iCountNoPrescale[it] % prescales[it] == 0 ) {
		  triggerBit[it] = true;
		  break;
		}
	      }
	    }
	  }

	  
	}
      }
      
      // 2. Loop to check overlaps
      for (int it = 0; it < Ntrig; it++){
	if (triggerBitNoPrescale[it]) {
	  (iCountNoPrescale[it])++;
	}
	if (triggerBit[it]) {
	  (iCount->at(it))++;
	  for (int it2 = 0; it2 < Ntrig; it2++){
	    if ( (it2<it) && triggerBit[it2] )
	      previousBitsFired[it] = true;
	    if ( (it2!=it) && triggerBit[it2] )
	      allOtherBitsFired[it] = true;    
	    if (triggerBit[it2])
	      (overlapCount->at(it))[it2] += 1;
	  }
	  if (!(previousBitsFired[it]))
	    (sPureCount->at(it))++;
	  if (!(allOtherBitsFired[it]))
	    (pureCount->at(it))++; 
	}
      }
   }
}
