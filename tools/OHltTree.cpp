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
			 ,map<TString,int> map_pathHLTPrescl
			 ,int NEntries
		   ,bool doMuonCut,bool doElecCut
			 ,double muonPt, double muonDr) {
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

	//int tempFlag;
	//TBranch *tempBranch;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;


    if (jentry%1000 == 0) cout<<"Processing entry "<<jentry<<"/"<<nentries<<"\r"<<flush<<endl;

    // 1. Loop to check which Bit fired
    // Triggernames are assigned to trigger cuts in unambigous way!
    // If you define a new trigger also define a new unambigous name!
		//SetStandardHLTPath();
		SetMapBitOfStandardHLTPath();

		bool passed_1jet_L1SingleJet15=false;


		///////////////////////////////////////////////////////////////////
		// Loop over HLT paths and do rate counting
		///////////////////////////////////////////////////////////////////
    for (int it = 0; it < Ntrig; it++){

      triggerBit[it] = false;
      triggerBitNoPrescale[it] = false;
      previousBitsFired[it] = false;
      allOtherBitsFired[it] = false;
			/*
			*/
			if ( (doMuonCut && MCmu3!=0) && ( trignames[it].Contains("Muon") || trignames[it].Contains("Mu") )  ) continue;
			if( (doElecCut && MCel1!=0) && ( trignames[it].Contains("Electron") || trignames[it].Contains("EM") )  ) continue;
			if( !passed_1jet_L1SingleJet15 && trignames[it].Contains("HLT2jetAve30") ) continue;
			//if ( trignames[it].Contains("Cand")   ) continue;

	  	if ( (map_BitOfStandardHLTPath.find(trignames[it])->second==1) ) { 

	    	triggerBitNoPrescale[it] = true;

	    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 

					if (trignames[it].CompareTo("HLT1jetPE7") == 0) passed_1jet_L1SingleJet15=true;

					triggerBit[it] = true; 

	    	} 
	  	}	
			else if (trignames[it].CompareTo("OpenHLT1jet") == 0) {

			  if (  (L1_SingleJet150==1)  &&  (recoJetCalPt[0]>200.)   ) {

			   	 triggerBitNoPrescale[it] = true;

			   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
			   	   triggerBit[it] = true;
			   	 }

				} // end if L1_SinleJet150
			}	  
			else if (trignames[it].CompareTo("OpenHLT2jet") == 0) {
			  if ( (recoJetCalPt[0]>150. && recoJetCalPt[1]>150.) && (L1_SingleJet150==1 || L1_DoubleJet70==1)  ) {
			    triggerBitNoPrescale[it] = true;
			    if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
			      triggerBit[it] = true;
			    }
			  }
			}	  
			else if (trignames[it].CompareTo("OpenHLT1MuonIso") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu7==1 ) {
				
		
					if(AmuonPassedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
					else if (HLT1MuonIso == 1) {
		
						cout << " Warning OpenHLT1MuonIso NOT passed, but HLT1MuonIso passed !!!! " << endl;
						PrintOhltVariables(2,muon);
						PrintOhltVariables(3,muon);
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("OpenHLT2MuonNonIso") == 0) { 
		
				// L1 Seed
			  if ( L1_DoubleMu3==1 ) {
				
		
					if(DoubleMuonPassedRelaxedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
					//else if( (HLT2MuonNonIso==1) ) {
								//cout << "******** WARNING:  OpenHLT2MuonNonIso not passed,  HLT2MuonNonIso passed *******************" << endl;
								//cout << "run " << run << ", event " << event << endl;
								//PrintOhltVariables(2,muon);
								//PrintOhltVariables(3,muon);
								//cout << "**********************************************************************************" << endl;
					//}
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("OpenHLT1Electron") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG12==1 ) {
				
					if(AelectronPassedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
					//else if( (HLT1Electron==1) ) {
								//cout << "******** WARNING:  OpenHLT1Electron not passed,  HLT1Electron passed *******************" << endl;
								//cout << "run " << run << ", event " << event << endl;
								//PrintOhltVariables(3,electron);
								//cout << "**********************************************************************************" << endl;
					//}
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("OpenHLT2ElectronRelaxed") == 0) { 
		
				// L1 Seed
			  if ( L1_DoubleEG10==1 ) {
				
					if(DoubleElectronPassedRelaxedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
					//else if( (HLT2ElectronRelaxed==1) ) {
								//cout << "******** WARNING:  OpenHLT2ElectronRelaxed not passed,  HLT2ElectronRelaxed passed *******************" << endl;
								//cout << "run " << run << ", event " << event << endl;
								//PrintOhltVariables(3,electron);
								//cout << "**********************************************************************************" << endl;
					//}
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("OpenHLT1Photon") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG12==1 ) {
				
		
					if(AphotonPassedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
					//else if( (HLT1Photon==1) ) {
								//cout << "******** WARNING:  OpenHLT1Photon not passed,  HLT1Photon passed *******************" << endl;
								//cout << "run " << run << ", event " << event << endl;
								//PrintOhltVariables(3,photon);
								//cout << "**********************************************************************************" << endl;
					//}
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("OpenHLT1PhotonRelaxed") == 0) { 
		
		
				// L1 Seed
			  if ( L1_SingleEG15==1 ) {
				
		
					if(AphotonPassedRelaxedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
					//else if( (HLT1PhotonRelaxed==1) ) {
								//cout << "******** WARNING:  OpenHLT1PhotonRelaxed not passed,  HLT1PhotonRelaxed passed *******************" << endl;
								//cout << "run " << run << ", event " << event << endl;
								//PrintOhltVariables(3,photon);
								//cout << "**********************************************************************************" << endl;
					//}
				} // end if L1 seed
				} // end else if
				else if (trignames[it].CompareTo("AptHLT1MET20") == 0) {
		
					  if (  (L1_ETM20==1)  &&  (recoMetCal>20)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						}
				}	  
				else if (trignames[it].CompareTo("AptHLT1MET35") == 0) {
		
					  if (  (L1_ETM30==1)  &&  (recoMetCal>35)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						}
				}	  
				else if (trignames[it].CompareTo("AptHLT1MET50") == 0) {
		
					  if (  (L1_ETM40==1)  &&  (recoMetCal>50)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						}
				}	  
				else if (trignames[it].CompareTo("AptHLT1MET65") == 0) {
		
					  if (  (L1_ETM40==1)  &&  (recoMetCal>65)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						}
				}	  
				else if (trignames[it].CompareTo("AptHLT1MET75") == 0) {
		
					  if (  (L1_ETM40==1)  &&  (recoMetCal>75)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} 
				}	  
				else if (trignames[it].CompareTo("AptHLT1Level1jet10") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >=10 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >=10 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >=10 ) rc = true;
						}
		
					  if (  (rc==true) /* &&  (recoJetCalPt[0]>????) */ ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1Level1jet15") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >=15 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >=15 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >=15 ) rc = true;
						}
		
					  if (  (rc==true) /* &&  (recoJetCalPt[0]>????) */ ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1Level1jet20") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >= 20 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >= 20 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >= 20 ) rc = true;
						}
		
					  if (  (rc==true) /* &&  (recoJetCalPt[0]>????) */ ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1Level1jet30") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >= 30 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >= 30 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >= 30 ) rc = true;
						}
		
					  if (  (rc==true) /* &&  (recoJetCalPt[0]>????) */ ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1Level1jet50") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >=50 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >=50 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >=50 ) rc = true;
						}
		
					  if (  (rc==true) /* &&  (recoJetCalPt[0]>????) */ ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1Level1jet70") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >=70 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >=70 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >=70 ) rc = true;
						}
		
					  if (  (rc==true) /* &&  (recoJetCalPt[0]>????) */ ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1Level1jet100") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >=100 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >=100 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >=100 ) rc = true;
						}
		
					  if (  (rc==true) /* &&  (recoJetCalPt[0]>????) */ ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1Level1jet150") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >=150 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >=150 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >=150 ) rc = true;
						}
		
					  if (  (rc==true) /* &&  (recoJetCalPt[0]>????) */ ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1Level1jet200") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >=200 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >=200 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >=200 ) rc = true;
						}
		
					  if (  (rc==true) /* &&  (recoJetCalPt[0]>????) */ ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1jet20") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >=10 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >=10 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >=10 ) rc = true;
						}
		
					  if (  (rc==true)  &&  (NrecoJetCal>0 && recoJetCalPt[0]>=20) ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1jet30") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >=15 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >=15 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >=15 ) rc = true;
						}
		
					  if (  (rc==true)  &&  (NrecoJetCal>0 && recoJetCalPt[0]>=30)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1jet40") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >= 15 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >= 15 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >= 15 ) rc = true;
						}
		
					  if (  (rc==true)  &&  (NrecoJetCal>0 && recoJetCalPt[0]>=40)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1jet50") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >= 30 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >= 30 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >= 30 ) rc = true;
						}
		
					  if (  (rc==true)  &&  (NrecoJetCal>0 && recoJetCalPt[0]>=50)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1jet60") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >= 30 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >= 30 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >= 30 ) rc = true;
						}
		
					  if (  (rc==true)  &&  (NrecoJetCal>0 && recoJetCalPt[0]>=60)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1jet80") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >=50 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >=50 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >=50 ) rc = true;
						}
		
					  if (  (rc==true)  &&  (NrecoJetCal>0 && recoJetCalPt[0]>=80)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1jet110") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >=70 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >=70 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >=70 ) rc = true;
						}
		
					  if (  (rc==true)  &&  (NrecoJetCal>0 && recoJetCalPt[0]>=110)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1jet150") == 0) {

						bool rc = false;

						if(NL1CenJet > 0) {

							if(L1CenJetEt[0] >=100 ) rc = true;
						}
						if(NL1ForJet > 0) {

							if(L1ForJetEt[0] >=100 ) rc = true;
						}
						if(NL1Tau> 0) {

							if(L1TauEt[0] >=100 ) rc = true;
						}
		
					  if (  (rc==true)  &&  (NrecoJetCal>0 && recoJetCalPt[0]>=150)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1jet180") == 0) {

					  if (  L1_SingleJet70  &&  (NrecoJetCal>0 && recoJetCalPt[0]>=180)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT1jet250") == 0) {

					  if (  L1_SingleJet70  &&  (NrecoJetCal>0 && recoJetCalPt[0]>=250)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT2jetAve15") == 0) {

					  if (  L1_SingleJet15  &&  (NrecoJetCal>1 && (recoJetCalPt[0]+recoJetCalPt[1])/2 >=15)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT2jetAve30") == 0) {

					  if (  L1_SingleJet30  &&  (NrecoJetCal>1 && (recoJetCalPt[0]+recoJetCalPt[1])/2 >=30)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT2jetAve50") == 0) {

					  if (  L1_SingleJet50  &&  (NrecoJetCal>1 && (recoJetCalPt[0]+recoJetCalPt[1])/2 >=50)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT2jetAve70") == 0) {

					  if (  L1_SingleJet70  &&  (NrecoJetCal>1 && (recoJetCalPt[0]+recoJetCalPt[1])/2 >=70)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT2jetAve130") == 0) {

					  if (  L1_SingleJet70  &&  (NrecoJetCal>1 && (recoJetCalPt[0]+recoJetCalPt[1])/2 >=130)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT2jetAve220") == 0) {

					  if (  L1_SingleJet70  &&  (NrecoJetCal>1 && (recoJetCalPt[0]+recoJetCalPt[1])/2 >=220)  ) {
		
					   	 triggerBitNoPrescale[it] = true;
		
					   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
					   	   triggerBit[it] = true;
					   	 }
		
						} // end if 
				}	  
				else if (trignames[it].CompareTo("AptHLT2jet100") == 0) {
	  			if ( (recoJetCalPt[0]>100. && recoJetCalPt[1]>100.) && (L1_SingleJet100==1 || L1_DoubleJet70==1)  ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLT4jet20") == 0) {
	  			if ( (recoJetCalPt[0]>20. && recoJetCalPt[1]>20. && recoJetCalPt[2]>20. && recoJetCalPt[3]>20.) && (HasL1jet(10) )  ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLT4jet30") == 0) {
	  			if ( (recoJetCalPt[0]>30. && recoJetCalPt[1]>30. && recoJetCalPt[2]>30. && recoJetCalPt[3]>30.) && (HasL1jet(15) )  ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLT5jet20") == 0) {
	  			if ( (recoJetCalPt[0]>20. && recoJetCalPt[1]>20. && recoJetCalPt[2]>20. && recoJetCalPt[3]>20. && recoJetCalPt[4]>20.) && (HasL1jet(10) )  ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLT5jet30") == 0) {
	  			if ( (recoJetCalPt[0]>30. && recoJetCalPt[1]>30. && recoJetCalPt[2]>30. && recoJetCalPt[3]>30. && recoJetCalPt[4]>30.) && (HasL1jet(15) )  ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLT6jet20") == 0) {
	  			if ( (recoJetCalPt[0]>20. && recoJetCalPt[1]>20. && recoJetCalPt[2]>20. && recoJetCalPt[3]>20. && recoJetCalPt[4]>20. && recoJetCalPt[5]>20.) && (HasL1jet(10) )  ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLT6jet30") == 0) {
	  			if ( (recoJetCalPt[0]>30. && recoJetCalPt[1]>30. && recoJetCalPt[2]>30. && recoJetCalPt[3]>30. && recoJetCalPt[4]>30. && recoJetCalPt[5]>30.) && (HasL1jet(15) )  ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLT1RelaxedLepton3jet20") == 0) {
	  			if ( (recoJetCalPt[0]>20. && recoJetCalPt[1]>20. && recoJetCalPt[2]>20. &&  HasL1jet(10) )  && (AmuonPassedAptOhltCuts(3,0.02,0) || AelectronPassedAptOhltCuts(5,0,0.06,3,1.5,2.45) )   ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLT1RelaxedLepton3jet30") == 0) {
	  			if ( (recoJetCalPt[0]>30. && recoJetCalPt[1]>30. && recoJetCalPt[2]>30.  && HasL1jet(15) ) && (AmuonPassedAptOhltCuts(3,0.02,0) || AelectronPassedAptOhltCuts(5,0,0.06,3,1.5,2.45) )  ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLTXRelMu3_3jet30") == 0) {
	  			if ( (recoJetCalPt[0]>30. && recoJetCalPt[1]>30. && recoJetCalPt[2]>30.  && HasL1jet(15) ) && (AmuonPassedAptOhltCuts(3,0.02,0)  )  ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLTXRelEl5_3jet30") == 0) {
	  			if ( (recoJetCalPt[0]>30. && recoJetCalPt[1]>30. && recoJetCalPt[2]>30.  && HasL1jet(15) ) && ( AelectronPassedAptOhltCuts(5,0,0.06,3,1.5,2.45) )  ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLT1RelaxedLepton4jet20") == 0) {
	  			if ( (recoJetCalPt[0]>20. && recoJetCalPt[1]>20. && recoJetCalPt[2]>20. && recoJetCalPt[3]>20. &&  HasL1jet(10) )  && (AmuonPassedAptOhltCuts(3,0.02,0) || AelectronPassedAptOhltCuts(5,0,0.06,3,1.5,2.45) )   ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLT1RelaxedLepton4jet30") == 0) {
	  			if ( (recoJetCalPt[0]>30. && recoJetCalPt[1]>30. && recoJetCalPt[2]>30. && recoJetCalPt[3]>30. && HasL1jet(15) ) && (AmuonPassedAptOhltCuts(3,0.02,0) || AelectronPassedAptOhltCuts(5,0,0.06,3,1.5,2.45) )  ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLT1RelaxedLepton5jet20") == 0) {
	  			if ( (recoJetCalPt[0]>20. && recoJetCalPt[1]>20. && recoJetCalPt[2]>20. && recoJetCalPt[3]>20.  && recoJetCalPt[4]>20. &&  HasL1jet(10) )  && (AmuonPassedAptOhltCuts(3,0.02,0) || AelectronPassedAptOhltCuts(5,0,0.06,3,1.5,2.45) )   ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				else if (trignames[it].CompareTo("AptHLT1RelaxedLepton5jet30") == 0) {
	  			if ( (recoJetCalPt[0]>30. && recoJetCalPt[1]>30. && recoJetCalPt[2]>30. && recoJetCalPt[3]>30. && recoJetCalPt[4]>30. && HasL1jet(15) ) && (AmuonPassedAptOhltCuts(3,0.02,0) || AelectronPassedAptOhltCuts(5,0,0.06,3,1.5,2.45) )  ) {
	    			triggerBitNoPrescale[it] = true;
	    			if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
	      			triggerBit[it] = true;
	    			}
	  			}
				}	  
				/*
				*/
				else if (trignames[it].CompareTo("AptHLT1MuonLevel1Open") == 0) { 
		
				// L1 Seed
				  if ( CandHLT1MuonLevel1==1 ) {
				
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
		    		} 
		
					} // end if L1 seed
				} // end else if
				else if (trignames[it].CompareTo("AptHLT1MuonLevel1") == 0) { 
		
				// L1 Seed
				  if ( L1_SingleMu7==1 || L1_DoubleMu3==1 ) {
				
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
		    		} 
		
					} // end if L1 seed
				} // end else if
				else if (trignames[it].CompareTo("AptHLT2Muon3") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(DoubleMuonPassedRelaxedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1MuonIso3_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(3,0.02,1)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1MuonIso5_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(5,0.02,1)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1MuonIso7_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(7,0.02,1)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1MuonIso9_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(9,0.02,1)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1MuonIso11_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(11,0.02,1)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1MuonIso3") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(3,999999,1)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1MuonIso5") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(5,999999,1)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1MuonIso7") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(7,999999,1)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1MuonIso9") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(9,999999,1)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1MuonIso11") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(11,999999,1)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
				else if (trignames[it].CompareTo("AptHLT1Muon3_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(3,0.02)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Muon5_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(5,0.02)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Muon7_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(7,0.02)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Muon9_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(9,0.02)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Muon11_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(11,0.02)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Muon13_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(13,0.02)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Muon15_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(15,0.02)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
				else if (trignames[it].CompareTo("AptHLT1Muon3") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(3,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Muon5") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(5,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Muon7") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu3==1 ) {
				
					if(AmuonPassedAptOhltCuts(7,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Muon9") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu7==1 ) {
				
					if(AmuonPassedAptOhltCuts(9,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Muon11") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu7==1 ) {
				
					if(AmuonPassedAptOhltCuts(11,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Muon13") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu7==1 ) {
				
					if(AmuonPassedAptOhltCuts(13,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Muon15") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu7==1 ) {
				
					if(AmuonPassedAptOhltCuts(15,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1L2Muon11") == 0) { 
		
				double muonPt = 11;
				// L1 Seed
			  if ( L1_SingleMu7==1 ) {
				
					bool rcL2 = false;

					// Loop over all L2 muons
					for (int i=0;i<NohMuL2;i++) {


							// L2 condition
						if ( TMath::Abs(ohMuL2Eta[i])<2.5 && ohMuL2Pt[i]+3.9*ohMuL2PtErr[i]*ohMuL2Pt[i]>muonPt )  {

							rcL2 = true;

						} // end if L2
					} //end for all muons in L2

					if(rcL2) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1L2Muon16") == 0) { 
		
				double muonPt = 16;
				// L1 Seed
			  if ( L1_SingleMu7==1 ) {
				
					bool rcL2 = false;

					// Loop over all L2 muons
					for (int i=0;i<NohMuL2;i++) {


							// L2 condition
						if ( TMath::Abs(ohMuL2Eta[i])<2.5 && ohMuL2Pt[i]+3.9*ohMuL2PtErr[i]*ohMuL2Pt[i]>muonPt )  {

							rcL2 = true;

						} // end if L2
					} //end for all muons in L2

					if(rcL2) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron5_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(5)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron7_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(7)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron8_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(8)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron8_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(8,0,0.06,3,999,999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron10_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG8==1 ) {
				
					if(AelectronPassedAptOhltCuts(10)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron9_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(9)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron11_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(11)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron13_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(13)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron15_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(15)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron20_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(20)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron25_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(25)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron30_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(30)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1ElectronIso8_10_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG8==1 ) {
				
					if(AelectronPassedAptOhltCuts(10,1,0.06,3,999,999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1ElectronIso8_12_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG8==1 ) {
				
					if(AelectronPassedAptOhltCuts(12,1,0.06,3,999,999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1ElectronIso8_15_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG8==1 ) {
				
					if(AelectronPassedAptOhltCuts(15,1,0.06,3,999,999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1ElectronIso8_10_DblThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG8==1 ) {
				
					if(AelectronPassedAptOhltCuts(10,1,0.12,6,999,999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1ElectronIso8_12_DblThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG8==1 ) {
				
					if(AelectronPassedAptOhltCuts(12,1,0.12,6,999,999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1ElectronIso8_15_DblThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG8==1 ) {
				
					if(AelectronPassedAptOhltCuts(15,1,0.12,6,999,999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron12_15_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG12==1 ) {
				
					if(AelectronPassedAptOhltCuts(15,0,0.06,3,999,999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron12_17_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG12==1 ) {
				
					if(AelectronPassedAptOhltCuts(17,0,0.06,3,999,999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron12_15_DblThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG12==1 ) {
				
					if(AelectronPassedAptOhltCuts(15,0,0.12,6,999,999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron12_17_DblThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG12==1 ) {
				
					if(AelectronPassedAptOhltCuts(17,0,0.12,6,999,999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT2ElectronIso8_10_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_DoubleIsoEG8==1 ) {
				
					if(AelectronPassedAptOhltCuts(10)>1) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT2Electron10_12_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_DoubleEG10==1 ) {
				
					if(AelectronPassedAptOhltCuts(12)>1) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT2Electron5_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(5) > 1) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT2Electron10_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(10) > 1) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Electron12_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AelectronPassedAptOhltCuts(12)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon5_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AphotonPassedAptOhltCuts(5,0,0)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon10_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AphotonPassedAptOhltCuts(10,0,0)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon12_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG8==1 ) {
				
					if(AphotonPassedAptOhltCuts(12,0,0)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon15_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG15==1 ) {
				
					if(AphotonPassedAptOhltCuts(15,0,0)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon20_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AphotonPassedAptOhltCuts(20,0,0)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon25_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AphotonPassedAptOhltCuts(25,0,0)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon30_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG20==1 ) {
				
					if(AphotonPassedAptOhltCuts(30,0,0)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon5NoTrkIso") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AphotonPassedAptOhltCuts(5,0,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon10_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AphotonPassedAptOhltCuts(10,0,0,1.5,6,4)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon10NoTrkIso") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG8==1 ) {
				
					if(AphotonPassedAptOhltCuts(10,0,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon12NoTrkIso") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG5==1 ) {
				
					if(AphotonPassedAptOhltCuts(12,0,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon15NoTrkIso") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG10==1 ) {
				
					if(AphotonPassedAptOhltCuts(15,0,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon20NoTrkIso") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG10==1 ) {
				
					if(AphotonPassedAptOhltCuts(20,0,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon25NoTrkIso") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG15==1 ) {
				
					if(AphotonPassedAptOhltCuts(25,0,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon30NoTrkIso") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG20==1 ) {
				
					if(AphotonPassedAptOhltCuts(30,0,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon35NoTrkIso") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG25==1 ) {
				
					if(AphotonPassedAptOhltCuts(35,0,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon40NoTrkIso") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG25==1 ) {
				
					if(AphotonPassedAptOhltCuts(40,0,999999)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1PhotonIso8_15_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG8==1 ) {
				
					if(AphotonPassedAptOhltCuts(15,1,0,1.5,6,4)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1PhotonIso8_17_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG8==1 ) {
				
					if(AphotonPassedAptOhltCuts(17,1,0,1.5,6,4)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1PhotonIso8_20_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG8==1 ) {
				
					if(AphotonPassedAptOhltCuts(20,1,0,1.5,6,4)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1PhotonIso8_20_DblThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG8==1 ) {
				
					if(AphotonPassedAptOhltCuts(20,1,0,3.0,12,8)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1PhotonIso8_25_DblThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG8==1 ) {
				
					if(AphotonPassedAptOhltCuts(25,1,0,3.0,12,8)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1PhotonIso8_30_DblThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG8==1 ) {
				
					if(AphotonPassedAptOhltCuts(30,1,0,3.0,12,8)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon12_20_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG12==1 ) {
				
					if(AphotonPassedAptOhltCuts(20,0,0,1.5,6,4)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon12_25_DefThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG12==1 ) {
				
					if(AphotonPassedAptOhltCuts(25,0,0,1.5,6,4)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon12_15_DblThrr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG12==1 ) {
				
					if(AphotonPassedAptOhltCuts(15,0,0,3.0,12,8)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon12_25_DblThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG12==1 ) {
				
					if(AphotonPassedAptOhltCuts(25,0,0,3.0,12,8)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon12_30_DblThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG12==1 ) {
				
					if(AphotonPassedAptOhltCuts(30,0,0,3.0,12,8)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon12_35_DblThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG12==1 ) {
				
					if(AphotonPassedAptOhltCuts(35,0,0,3.0,12,8)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT1Photon12_40_DblThr") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleEG12==1 ) {
				
					if(AphotonPassedAptOhltCuts(40,0,0,3.0,12,8)) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT2PhotonIso8_15_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_DoubleIsoEG8==1 ) {
				
					if(AphotonPassedAptOhltCuts(15,1,0)>1) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT2Photon10NoTrkIso") == 0) { 
		
				// L1 Seed
			  if ( L1_DoubleEG5==1 ) {
				
					if(AphotonPassedAptOhltCuts(10,0,99999)>1) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT2Photon10_20_Open") == 0) { 
		
				// L1 Seed
			  if ( L1_DoubleEG10==1 ) {
				
					if(AphotonPassedAptOhltCuts(20,0,0)>1) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT2Photon10_20NoTrkIso") == 0) { 
		
				// L1 Seed
			  if ( L1_DoubleEG10==1 ) {
				
					if(AphotonPassedAptOhltCuts(20,0,999)>1) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("AptHLT2Photon10_30NoTrkIso") == 0) { 
		
				// L1 Seed
			  if ( L1_DoubleEG10==1 ) {
				
					if(AphotonPassedAptOhltCuts(30,0,999)>1) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
		
				} // end if L1 seed
			} // end else if

			/*
			////////////////////////////////////////////////////
			// Section for Validating OpenHLT
			////////////////////////////////////////////////////
			else if (trignames[it].CompareTo("OpenHLT1MuonIso") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleMu7==1 ) {
				
		
					if(AmuonPassedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
					//else if( (HLT1MuonIso==1) ) {
								//cout << "******** WARNING:  OpenHLT1MuonIso not passed,  HLT1MuonIso passed *******************" << endl;
								//cout << "run " << run << ", event " << event << endl;
								//PrintOhltVariables(2,muon);
								//PrintOhltVariables(3,muon);
								//cout << "**********************************************************************************" << endl;
					//}
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("OpenHLT2MuonNonIso") == 0) { 
		
				// L1 Seed
			  if ( L1_DoubleMu3==1 ) {
				
		
					if(DoubleMuonPassedRelaxedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
					//else if( (HLT2MuonNonIso==1) ) {
								//cout << "******** WARNING:  OpenHLT2MuonNonIso not passed,  HLT2MuonNonIso passed *******************" << endl;
								//cout << "run " << run << ", event " << event << endl;
								//PrintOhltVariables(2,muon);
								//PrintOhltVariables(3,muon);
								//cout << "**********************************************************************************" << endl;
					//}
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("OpenHLT1jet") == 0) {
			  if (  (L1_SingleJet150==1)  ) {
			  	if ( (recoJetCalPt[0]>200.)   ) {
			   	 triggerBitNoPrescale[it] = true;
			   	 if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
			   	   triggerBit[it] = true;
			   	 }
			  	}
					//else if( (HLT1jet==1) ) {
								//cout << "******** WARNING:  OpenHLT1jet not passed,  HLT1jet passed *******************" << endl;
								//cout << "run " << run << ", event " << event << endl;
								//PrintOhltVariables(3,jet);
								//cout << "**********************************************************************************" << endl;
					//}
				} // end if L1_SinleJet150
			}	  
			else if (trignames[it].CompareTo("OpenHLT2jet") == 0) {
			  if ( (recoJetCalPt[0]>150. && recoJetCalPt[1]>150.) && (L1_SingleJet150==1 || L1_DoubleJet70==1)  ) {
			    triggerBitNoPrescale[it] = true;
			    if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) {
			      triggerBit[it] = true;
			    }
			  }
			}	  
			else if (trignames[it].CompareTo("OpenHLT1Electron") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG12==1 ) {
				
					if(AelectronPassedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
					//else if( (HLT1Electron==1) ) {
								//cout << "******** WARNING:  OpenHLT1Electron not passed,  HLT1Electron passed *******************" << endl;
								//cout << "run " << run << ", event " << event << endl;
								//PrintOhltVariables(3,electron);
								//cout << "**********************************************************************************" << endl;
					//}
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("OpenHLT2ElectronRelaxed") == 0) { 
		
				// L1 Seed
			  if ( L1_DoubleEG10==1 ) {
				
					if(DoubleElectronPassedRelaxedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
					//else if( (HLT2ElectronRelaxed==1) ) {
								//cout << "******** WARNING:  OpenHLT2ElectronRelaxed not passed,  HLT2ElectronRelaxed passed *******************" << endl;
								//cout << "run " << run << ", event " << event << endl;
								//PrintOhltVariables(3,electron);
								//cout << "**********************************************************************************" << endl;
					//}
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("OpenHLT1Photon") == 0) { 
		
				// L1 Seed
			  if ( L1_SingleIsoEG12==1 ) {
				
		
					if(AphotonPassedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
					//else if( (HLT1Photon==1) ) {
								//cout << "******** WARNING:  OpenHLT1Photon not passed,  HLT1Photon passed *******************" << endl;
								//cout << "run " << run << ", event " << event << endl;
								//PrintOhltVariables(3,photon);
								//cout << "**********************************************************************************" << endl;
					//}
				} // end if L1 seed
			} // end else if
			else if (trignames[it].CompareTo("OpenHLT1PhotonRelaxed") == 0) { 
		
		
				// L1 Seed
			  if ( L1_SingleEG15==1 ) {
				
		
					if(AphotonPassedRelaxedOhltCuts()) {
		
			    	triggerBitNoPrescale[it] = true;
		
			    	if ((iCountNoPrescale[it]) % map_pathHLTPrescl.find(trignames[it])->second == 0) { 
		
			      	triggerBit[it] = true; 
		
			    	} 
		
					}
					//else if( (HLT1PhotonRelaxed==1) ) {
								//cout << "******** WARNING:  OpenHLT1PhotonRelaxed not passed,  HLT1PhotonRelaxed passed *******************" << endl;
								//cout << "run " << run << ", event " << event << endl;
								//PrintOhltVariables(3,photon);
								//cout << "**********************************************************************************" << endl;
					//}
				} // end if L1 seed
			} // end else if
			*/

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

void OHltTree::PrintOhltVariables(int level, int type)
{
	switch(type) 
	{	
		case muon:

				if(level == 3) {

					cout << "Level 3: number of muons = " << NohMuL3 << endl;

					for (int i=0;i<NohMuL3;i++) {

						cout << "ohMuL3Pt["<<i<<"] = " << ohMuL3Pt[i] << endl;
						cout << "ohMuL3PtErr["<<i<<"] = " << ohMuL3PtErr[i] << endl;
						cout << "ohMuL3Pt+Err["<<i<<"] = " << ohMuL3Pt[i]+2.2*ohMuL3PtErr[i]*ohMuL3Pt[i] << endl;
			// swapped Eta and Phi in HltMuon.cc
						cout << "ohMuL3Phi["<<i<<"] = " << ohMuL3Phi[i] << endl;
						cout << "ohMuL3Eta["<<i<<"] = " << ohMuL3Eta[i] << endl;
			// remember to change it back when analyzing outputs of next OpenHLT tag
						cout << "ohMuL3Chg["<<i<<"] = " << ohMuL3Chg[i] << endl;
						cout << "ohMuL3Iso["<<i<<"] = " << ohMuL3Iso[i] << endl;
						cout << "ohMuL3Dr["<<i<<"] = " << ohMuL3Dr[i] << endl;
						cout << "ohMuL3Dz["<<i<<"] = " << ohMuL3Dz[i] << endl;

					}
				}
				else if(level == 2) {

						cout << "Level 2: number of muons = " << NohMuL2 << endl;
					for (int i=0;i<NohMuL2;i++) {

						cout << "ohMuL2Pt["<<i<<"] = " << ohMuL2Pt[i] << endl;
						cout << "ohMuL2PtErr["<<i<<"] = " << ohMuL2PtErr[i] << endl;
						cout << "ohMuL2Pt+Err["<<i<<"] = " << ohMuL2Pt[i]+3.9*ohMuL2PtErr[i]*ohMuL2Pt[i] << endl;
						cout << "ohMuL2Phi["<<i<<"] = " << ohMuL2Phi[i] << endl;
						cout << "ohMuL2Eta["<<i<<"] = " << ohMuL2Eta[i] << endl;
						cout << "ohMuL2Chg["<<i<<"] = " << ohMuL2Chg[i] << endl;
						cout << "ohMuL2Iso["<<i<<"] = " << ohMuL2Iso[i] << endl;
						cout << "ohMuL2Dr["<<i<<"] = " << ohMuL2Dr[i] << endl;
						cout << "ohMuL2Dz["<<i<<"] = " << ohMuL2Dz[i] << endl;

					}
				}
				else {
					cout << "PrintOhltVariables: Ohlt has Muon variables only for L2 and 3. Must provide one." << endl;
				}
				break;

	case electron:

		cout << "oh: number of electrons = " << NohEle << endl;
		for (int i=0;i<NohEle;i++) {

			cout << "ohEleEt["<<i<<"] = " << ohEleEt[i] << endl;
			cout << "ohElePhi["<<i<<"] = " << ohElePhi[i] << endl;
			cout << "ohEleEta["<<i<<"] = " << ohEleEta[i] << endl;
			cout << "ohEleE["<<i<<"] = " << ohEleE[i] << endl;
			cout << "ohEleP["<<i<<"] = " << ohEleP[i] << endl;
		 	cout << "ohElePt["<<i<<"] =" <<  ohEleP[i] * TMath::Sin(2*TMath::ATan(TMath::Exp(-1*ohEleEta[i]))) << endl;
			cout << "ohEleHiso["<<i<<"] = " << ohEleHiso[i] << endl;
			cout << "ohEleTiso["<<i<<"] = " << ohEleTiso[i] << endl;
			cout << "ohEleL1iso["<<i<<"] = " << ohEleL1iso[i] << endl;
			cout << "recoElecE["<<i<<"] = " << recoElecE[i] << endl;
			cout << "recoElecEt["<<i<<"] = " << recoElecEt[i] << endl;
			cout << "recoElecPt["<<i<<"] = " << recoElecPt[i] << endl;
			cout << "recoElecPhi["<<i<<"] = " << recoElecPhi[i] << endl;
			cout << "recoElecEta["<<i<<"] = " << recoElecEta[i] << endl;

		}
			break;

	case photon:

		cout << "oh: number of photons = " << NohPhot << endl;
		for (int i=0;i<NohPhot;i++) {

			cout << "ohPhotEt["<<i<<"] = " << ohPhotEt[i] << endl;
			cout << "ohPhotPhi["<<i<<"] = " << ohPhotPhi[i] << endl;
			cout << "ohPhotEta["<<i<<"] = " << ohPhotEta[i] << endl;
			cout << "ohPhotEiso["<<i<<"] = " << ohPhotEiso[i] << endl;
			cout << "ohPhotHiso["<<i<<"] = " << ohPhotHiso[i] << endl;
			cout << "ohPhotTiso["<<i<<"] = " << ohPhotTiso[i] << endl;
			cout << "ohPhotL1iso["<<i<<"] = " << ohPhotL1iso[i] << endl;
			cout << "recoPhotE["<<i<<"] = " << recoPhotE[i] << endl;
			cout << "recoPhotEt["<<i<<"] = " << recoPhotEt[i] << endl;
			cout << "recoPhotPt["<<i<<"] = " << recoPhotPt[i] << endl;
			cout << "recoPhotPhi["<<i<<"] = " << recoPhotPhi[i] << endl;
			cout << "recoPhotEta["<<i<<"] = " << recoPhotEta[i] << endl;

		}
			break;

	case jet:

		cout << "oh: number of recoJetCal = " << NrecoJetCal << endl;
		for (int i=0;i<NrecoJetCal;i++) {

			cout << "recoJetCalE["<<i<<"] = " << recoJetCalE[i] << endl;
			cout << "recoJetCalEt["<<i<<"] = " << recoJetCalEt[i] << endl;
			cout << "recoJetCalPt["<<i<<"] = " << recoJetCalPt[i] << endl;
			cout << "recoJetCalPhi["<<i<<"] = " << recoJetCalPhi[i] << endl;
			cout << "recoJetCalEta["<<i<<"] = " << recoJetCalEta[i] << endl;
		}
			break;

	/*
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
	*/


	default:

			cout << "PrintOhltVariables: You did not provide correct object type." <<endl;
			break;

	} // end switch
}

bool OHltTree::AmuonPassedOhltCuts()
{

	bool rcL2 = false;

	// Loop over all L2 muons
	for (int i=0;i<NohMuL2;i++) {


			// L2 condition
		if ( TMath::Abs(ohMuL2Eta[i])<2.5 && ohMuL2Pt[i]+3.9*ohMuL2PtErr[i]*ohMuL2Pt[i]>11. && ohMuL2Iso[i]==1)  {

			rcL2 = true;

		} // end if L2
	} //end for all muons in L2

	bool rcL3 = false;

	for (int i=0;i<NohMuL3;i++) {
				// L3 condition
		if( TMath::Abs(ohMuL3Eta[i])<2.5 && ohMuL3Pt[i]+2.2*ohMuL3PtErr[i]*ohMuL3Pt[i]>11. && ohMuL3Dr[i] < 0.02 && ohMuL3Iso[i]==1 )   { 

			rcL3 = true;

	 	}	// end if L3
	} // end for all muons in L3


	return (rcL2 && rcL3);


}

bool OHltTree::AmuonPassedAptOhltCuts(double muonPt, double muonDr, double muonIso)
{

	bool rcL2 = false;

	// Loop over all L2 muons
	for (int i=0;i<NohMuL2;i++) {


			// L2 condition
		if ( TMath::Abs(ohMuL2Eta[i])<2.5 && ohMuL2Pt[i]+3.9*ohMuL2PtErr[i]*ohMuL2Pt[i]>muonPt )  {

			rcL2 = true;

		} // end if L2
	} //end for all muons in L2

	bool rcL3 = false;

	for (int i=0;i<NohMuL3;i++) {
				// L3 condition
		if( TMath::Abs(ohMuL3Eta[i])<2.5 && ohMuL3Pt[i]+2.2*ohMuL3PtErr[i]*ohMuL3Pt[i]>muonPt && ohMuL3Dr[i] < muonDr && ohMuL3Iso[i] >= muonIso) { 

			rcL3 = true;

	 	}	// end if L3
	} // end for all muons in L3


	return (rcL2 && rcL3);


}

bool OHltTree::DoubleMuonPassedRelaxedOhltCuts()
{

	bool rcL2 = false;

	// at least 2 muons at L2
	if (NohMuL2<2) return rcL2;

	int nMuPassL2 = 0;
	// Loop over all L2 muons
	for (int i=0;i<NohMuL2;i++) {


			// L2 condition
		if ( TMath::Abs(ohMuL2Eta[i])<2.5 && ohMuL2Pt[i]+3.9*ohMuL2PtErr[i]*ohMuL2Pt[i]>3. )  {

			nMuPassL2++;

		} // end if L2
	} //end for all muons in L2

	if(nMuPassL2 >= 2) rcL2=true;


	bool rcL3 = false;

	// at least 2 muons at L3
	if (NohMuL3<2) return rcL3;

	int nMuPassL3 = 0;

	for (int i=0;i<NohMuL3;i++) {
				// L3 condition
				// here, Phi is really Eta, it was a bug in HLTMuon.cc
		if( TMath::Abs(ohMuL3Eta[i])<2.5 && ohMuL3Pt[i]+2.2*ohMuL3PtErr[i]*ohMuL3Pt[i]>3. /*&& ohMuL3Dr[i] < 0.02*/  )   { 

			nMuPassL3++;

	 	}	// end if L3
	} // end for all muons in L3

	if(nMuPassL3 >= 2) rcL3=true;

	return (rcL2 && rcL3);


}

bool OHltTree::AphotonPassedOhltCuts()
{

	bool rc = false;

	// Loop over all oh photons
	for (int i=0;i<NohPhot;i++) {

		if ( ohPhotEt[i] > 30 && ohPhotEiso[i] <1.5 && ohPhotTiso[i]==0 && ohPhotL1iso[i]==1 /*&& isEgammaL1_PhotonSuperclusterMatched(i)*/ ) {

			if( (TMath::Abs(ohPhotEta[i]) < 1.5 && ohPhotHiso[i]<6)  || 
					(1.5 < TMath::Abs(ohPhotEta[i]) && TMath::Abs(ohPhotEta[i]) < 2.5 && ohPhotHiso[i]<4) ) {

					rc =true;
		
			}
		}

	} // end for

	return rc;

}

int  OHltTree::AphotonPassedAptOhltCuts(double Et, double L1iso, double Tiso, double Eiso, double HisoBR, double HisoEC)
{

	int rc = 0;

	// Loop over all oh photons
	for (int i=0;i<NohPhot;i++) {

		if ( ohPhotEt[i] > Et)
			if ( ohPhotL1iso[i] >= L1iso )
				if( ohPhotTiso[i]<=Tiso )
					if( ohPhotEiso[i] < Eiso /* && isEgammaL1_PhotonSuperclusterMatched(i)*/ ) {

						if( (TMath::Abs(ohPhotEta[i]) < 1.5 && ohPhotHiso[i] < HisoEC )  || 
						(1.5 < TMath::Abs(ohPhotEta[i]) && TMath::Abs(ohPhotEta[i]) < 2.5 && ohPhotHiso[i] < HisoEC ) ) 

						rc++;
		
					}

	} // end for

	return rc;

}

bool OHltTree::AphotonPassedRelaxedOhltCuts()
{

	bool rc = false;

	// Loop over all oh photons
	for (int i=0;i<NohPhot;i++) {


		if ( ohPhotEt[i] > 40 && ohPhotEiso[i] <1.5 && ohPhotTiso[i]==0 /*&& isEgammaL1_PhotonSuperclusterMatched(i)*/ ) {

			if( (TMath::Abs(ohPhotEta[i]) < 1.5 && ohPhotHiso[i]<6)  || 
					(1.5 < TMath::Abs(ohPhotEta[i]) && TMath::Abs(ohPhotEta[i]) < 2.5 && ohPhotHiso[i]<4) ) {

					rc =true;
		
			}
		}

	} // end for

	return rc;

}

bool OHltTree::AelectronPassedOhltCuts()
{

	bool rc = false;
	//const int size = NohEle;
	//float ohElePt[size];

	double eoverpbarrelcut = 1.5;
	double eoverpendcapcut = 2.45;

	// Loop over all oh electrons
	for (int i=0;i<NohEle;i++) {

		//do not have pT available, must calculate it form P and eta
		//ohElePt[i] = ohEleP[i] * TMath::Sin(2*TMath::ATan(TMath::Exp(-1*ohEleEta[i])));

		
		if ( ohEleEt[i] > 15 && ohEleHiso[i] <3 && ohEleTiso[i] !=-999. && ohEleTiso[i]< 0.06 && ohEleL1iso[i]==1 /*&& isEgammaL1_ElectronSuperclusterMatched(i)*/) {

			if( (TMath::Abs(ohEleEta[i]) < 1.5 && ohEleE[i]/ohEleP[i]<eoverpbarrelcut)  || 
					(1.5 < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i]) < 2.5 && ohEleE[i]/ohEleP[i]<eoverpendcapcut) ) {

					rc =true;

			}
		
		} // end if

	} // end for

	return rc;

}

int OHltTree::AelectronPassedAptOhltCuts(double Et, double L1iso, double Tiso, double Hiso, double eoverpBR, double eoverpEC)
{

	int rc = 0;;
	//const int size = NohEle;
	//float ohElePt[size];

	//double eoverpbarrelcut = 1.5;
	//double eoverpendcapcut = 2.45;

	// Loop over all oh electrons
	for (int i=0;i<NohEle;i++) {

		//do not have pT available, must calculate it form P and eta
		//ohElePt[i] = ohEleP[i] * TMath::Sin(2*TMath::ATan(TMath::Exp(-1*ohEleEta[i])));

		
		if ( ohEleEt[i] >= Et) {

			if ( ohEleL1iso[i] >= L1iso )   // L1iso is 0 or 1

				if ( ohEleTiso[i] <= Tiso && ohEleTiso[i] != -999   && ohEleHiso[i] <= Hiso /*&& isEgammaL1_ElectronSuperclusterMatched(i)*/ ) 

					if( (TMath::Abs(ohEleEta[i]) <= 1.5 && ohEleE[i]/ohEleP[i] <= eoverpBR)  || 
							(1.5 < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i]) <= 2.5 && ohEleE[i]/ohEleP[i] <= eoverpEC) ) 

						rc++;

		} // end if

	} // end for

	return rc;

}

bool OHltTree::AelectronPassedRelaxedOhltCuts()
{

	bool rc = false;
	//const int size = NohEle;
	//float ohElePt[size];

	// Loop over all oh electrons
	for (int i=0;i<NohEle;i++) {

		//ohElePt[i] = ohEleP[i] * TMath::Sin(2*TMath::ATan(TMath::Exp(-1*ohEleEta[i])));

		if ( ohEleEt[i] > 17 && ohEleHiso[i] <3 && ohEleTiso[i]< 0.06) {

					rc =true;
		
		}

	} // end for

	return rc;

}

bool OHltTree::DoubleElectronPassedRelaxedOhltCuts()
{

	bool rc = false;
	//const int size = NohEle;
	//float ohElePt[size];

	double eoverpbarrelcut = 15000;  // Changed from 1.5 in CMSSW_16X
	double eoverpendcapcut = 24500;  // Changed from 2.45 in CMSSW_16X

	int nElePass = 0;
	// Loop over all oh electrons
	for (int i=0;i<NohEle;i++) {

		//do not have pT available, must calculate it form P and eta
		//ohElePt[i] = ohEleP[i] * TMath::Sin(2*TMath::ATan(TMath::Exp(-1*ohEleEta[i])));

		
		if ( ohEleEt[i] > 12.0 && ohEleHiso[i] <9.0 && ohEleTiso[i] !=-999. && ohEleTiso[i]< 0.4 /*&& isEgammaL1_ElectronSuperclusterMatched(i)*/) {

			if( (TMath::Abs(ohEleEta[i]) < 1.5 && ohEleE[i]/ohEleP[i]<eoverpbarrelcut)  || 
					(1.5 < TMath::Abs(ohEleEta[i]) && TMath::Abs(ohEleEta[i]) < 2.5 && ohEleE[i]/ohEleP[i]<eoverpendcapcut) ) {

					nElePass++;

			}
		
		} // end if

	} // end for

	if(nElePass>=2) rc =true;
	return rc;

}


bool OHltTree::isEgammaL1_ElectronSuperclusterMatched(int i) 
{

	bool rc = false;
	
	float deltaPhi = 1.044/2.;
	float deltaEta = -999.;
	
	//barrel
	if(TMath::Abs(ohEleEta[i]) < 1.4791) deltaEta = 0.522/2.;

	//end cap
	if(TMath::Abs(ohEleEta[i])> 1.4791 && TMath::Abs(ohEleEta[i]) < 2.5) deltaEta = 0.87/2.;

	// isolated case
	if(ohEleL1iso[i]==1) {

		for (int j=0;j<NL1IsolEm;j++) {

			if(L1IsolEmEt[j] > 5.)
				if(TMath::Abs(ohEleEta[i] - L1IsolEmEta[j]) < deltaEta  && TMath::Abs(ohElePhi[i] -L1IsolEmPhi[j]) < deltaPhi) rc = true;

		}

	}
	// nonisolated case
	else {

		for (int j=0;j<NL1IsolEm;j++) {

			if(TMath::Abs(ohEleEta[i] - L1NIsolEmEta[j]) < deltaEta  && TMath::Abs(ohElePhi[i] -L1NIsolEmPhi[j]) < deltaPhi) rc = true;

		}

	} // end else

	return rc;

}


bool OHltTree::isEgammaL1_PhotonSuperclusterMatched(int i) 
{

	bool rc = false;
	
	float deltaPhi = 1.044/2.;
	float deltaEta = -999.;
	
	//barrel
	if(TMath::Abs(ohPhotEta[i]) < 1.4791) deltaEta = 0.522/2.;

	//end cap
	if(TMath::Abs(ohPhotEta[i])> 1.4791 && TMath::Abs(ohPhotEta[i]) < 2.5) deltaEta = 0.87/2.;

	// isolated case
	if(ohPhotL1iso[i]==1) {

		for (int j=0;j<NL1IsolEm;j++) {

			if(TMath::Abs(ohPhotEta[i] - L1IsolEmEta[j]) < deltaEta  && TMath::Abs(ohPhotPhi[i] -L1IsolEmPhi[j]) < deltaPhi) rc = true;

		}

	}
	// nonisolated case
	else {

		for (int j=0;j<NL1IsolEm;j++) {

			if(TMath::Abs(ohPhotEta[i] - L1NIsolEmEta[j]) < deltaEta  && TMath::Abs(ohPhotPhi[i] -L1NIsolEmPhi[j]) < deltaPhi) rc = true;

		}

	} // end else

	return rc;

}

bool OHltTree::HasL1jet(double Pt) {

	bool rc = false;

	if(NL1CenJet > 0) {

		if(L1CenJetEt[0] >=Pt ) rc = true;
	}
	if(NL1ForJet > 0) {

		if(L1ForJetEt[0] >=Pt ) rc = true;
	}
	if(NL1Tau> 0) {

		if(L1TauEt[0] >=Pt ) rc = true;
	}

	return rc;

}
