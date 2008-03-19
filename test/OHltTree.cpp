#define OHltTree_cxx
#include "OHltTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <TFormula.h>

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

void OHltTree::Loop(vector<int> * iCount, vector<int> * sPureCount, vector<int> * pureCount
		   ,vector< vector<int> > * overlapCount
		   ,vector<TString> trignames
			 ,map<TString,int> map_TrigPrescls,int NEntries
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

	// Loop over events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    if (jentry%1000 == 0) cout<<"Processing entry "<<jentry<<"/"<<nentries<<"\r"<<flush<<endl;

		// Set reference to TBranch that contains bit of corresponding path
		SetMapBitOfStandardHLTPath();

    // Loop to check which bit fired
    for (int it = 0; it < Ntrig; it++){

      triggerBit[it] = false;
      triggerBitNoPrescale[it] = false;
      previousBitsFired[it] = false;
      allOtherBitsFired[it] = false;

			// Do not double overlap muon triggered events 
			if ( (doMuonCut && MCmu3!=0) && ( trignames[it].Contains("Muon") || trignames[it].Contains("Mu") )  ) continue;

			// Do not double overlap electron triggered events 
			if( (doElecCut && MCel1!=0) && ( trignames[it].Contains("Electron") || trignames[it].Contains("EM") )  ) continue;

	  	if ( (map_BitOfStandardHLTPath.find(trignames[it])->second==1) ) { 

	    	triggerBitNoPrescale[it] = true;

	    	if ((iCountNoPrescale[it]) % map_TrigPrescls.find(trignames[it])->second == 0) { 

					triggerBit[it] = true; 

	    	} 
	  	}	

    }  // end for it

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
	if (not previousBitsFired[it])
	  (sPureCount->at(it))++;
	if (not allOtherBitsFired[it])
	  (pureCount->at(it))++; 
      }
    }
  }

}

