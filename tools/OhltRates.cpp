/////////////////////////////////////////////////////////////////////////////////////////////////
//
//        Program to calculate rates of trigger paths using variables of OHltTree class,
//
//				Note: OHltTree class needs to be updated if any new variables become available 
//				in OpenHLT (HLTAnalyzer).
//				
//        Author:  Vladimir Rekovic,     Date: 2007/12/10
//					
//         
//
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include "OHltTree.h"
#include "OHltMenu.h"

#include "TH1.h"
#include "TChain.h"
#include "TCut.h"

#include <map>

using namespace std;

void BookOHltMenu_2e31_v02(OHltMenu *menu, double &iLumi, double &nBunches);
void BookOHltMenu_2e30_v02(OHltMenu *menu, double &iLumi, double &nBunches);
void BookOHltMenu_2e31_v03(OHltMenu *menu, double &iLumi, double &nBunches);
void BookOHltMenu_2e30_v03(OHltMenu *menu, double &iLumi, double &nBunches);
void BookOHltMenu_Validate(OHltMenu *menu, double &iLumi, double &nBunches);

Double_t eff(Int_t a, Int_t b){ 
  if (b==0.){return -1.;}
  Double_t af = Double_t(a);
  Double_t bf = Double_t(b);   
  Double_t effi = af/bf;
  return effi;
}
Double_t seff(Int_t a, Int_t b){
  if (b==0.){return -1.;}
  Double_t af = Double_t(a);
  Double_t bf = Double_t(b);   
  Double_t r = af/bf;
  Double_t unc = sqrt(af + (r*r*bf) )/bf;
  return unc;
}

Double_t eff(Double_t a, Double_t b){ 
  if (b==0.){return -1.;}
  Double_t af = Double_t(a);
  Double_t bf = Double_t(b);   
  Double_t effi = af/bf;
  return effi;
}

Double_t seff(Double_t a, Double_t b){
  if (b==0.){return -1.;}
  Double_t af = Double_t(a);
  Double_t bf = Double_t(b);   
  Double_t r = af/bf;
  Double_t unc = sqrt(af + (r*r*bf) )/bf;
  return unc;
}

void ShowUsage() {

		cout << "usage:  ./OhltRates nEvents menu <conditions> " << endl;

}

int main(int argc, char *argv[]){
  
	if(argc<3) {

		ShowUsage();
		return 0;

	}

  int NEntries = -1;
  if (argc>1) {
    NEntries = atoi(argv[1]);
  }
  TString sMenu = "";
  if (argc>2) {
    sMenu = TString(argv[2]);
  }
  TString sConditions = "";
  if (argc>3) {
    sConditions = TString(argv[3]);
  }
  int Version = 0;
	char sVersion[100];
  if (argc>4) {
    Version = atoi(argv[4]);
  }
  sprintf(sVersion,"_v%02d",Version);


  ////////////////////////////////////////////////////////////
  // Instanteneous Luminosity [cm^-2 s^-1]
  
  // const Double_t ILumi = 8.E32;
  // const Double_t ILumi = 2.E33;
  //const Double_t ILumi = 1.E32;
  
  // For accurate rate calculation
  const double bunchCrossingTime = 25.0E-09;  // 25 ns
  const double maxFilledBunches = 3557;
  //const double nFilledBunches = 1000;
  
  
  /**** Different Beam conditions: ****/
 
  double ILumi = 1.E27;
  double nFilledBunches = 1;
  
  //double ILumi = 2E31;
  //double nFilledBunches = 156;
  
  //const double ILumi = 1.7E30;
  //const double nFilledBunches = 43;
  
  //const double ILumi = 6.1E30;
  //const double nFilledBunches = 43;
  
  //const double ILumi = 1.1E31;
  //const double nFilledBunches = 156;
  
  //const double ILumi = 5.6E31;
  //const double nFilledBunches = 156;
  
  //const double ILumi = 1.1E32;
  //const double nFilledBunches = 156;
  
  //const double ILumi = 1.0E32;
  //const double nFilledBunches = 1000;
  
  /****   ***   ****/
	
	OHltMenu* menu = new OHltMenu();

	if(sMenu.CompareTo("2e31_v02") == 0) BookOHltMenu_2e31_v02(menu,ILumi,nFilledBunches);
	else if(sMenu.CompareTo("2e31_v03") == 0) BookOHltMenu_2e31_v03(menu,ILumi,nFilledBunches);
	else if(sMenu.CompareTo("2e30_v02") == 0) BookOHltMenu_2e30_v02(menu,ILumi,nFilledBunches);
	else if(sMenu.CompareTo("2e30_v03") == 0) BookOHltMenu_2e30_v03(menu,ILumi,nFilledBunches);
	else if(sMenu.CompareTo("Validate") == 0) BookOHltMenu_Validate(menu,ILumi,nFilledBunches);
	else {

		cout << "No valid menu specified.  Either creata a new menu or use existing one. Exiting!" << endl;
		ShowUsage();
		return 0;

	}


  vector<TString> trignames = menu->GetHlts(); 
	int Ntrig = (int) trignames.size();

  map<TString,int> map_TrigPrescls = menu->GetHltPrescaleMap(); 

  double collisionRate = (nFilledBunches / maxFilledBunches) / bunchCrossingTime ;  // Hz
  
  ////////////////////////////////////////////////////////////
	cout << endl << endl << endl;
	cout << "--------------------------------------------------------------------------" << endl;
	cout << "NEntries = " << NEntries << endl;
	cout << "Menu = " << sMenu << endl;
	if(sConditions.CompareTo("") ==0) sConditions=TString("Ideal");
	cout << "Conditions = " << sConditions << endl;
	cout << "Version = " << sVersion << endl;
	cout << "Inst Luminosity = " << ILumi << ",  Bunches = " << nFilledBunches <<  endl;
	cout << "--------------------------------------------------------------------------" << endl;
	cout << endl << endl << endl;

	// Wait 2 sec for user to read announcement
	sleep(2);




	/*
  ////////////////////////////////////////////////////////////
  // Cross-sections [pb]
      
	vector<Double_t> xsec;
	vector<Double_t> skmeff; // Skim efficiencies

	//xsec.push_back(1.01E5); // PYTHIA cross-section for QCD in 170 < ^pt < 230
	//skmeff.push_back(1.);

	xsec.push_back(7.923E10); // PYTHIA cross-section for MinBias (pb)
	skmeff.push_back(0.105);  // skmeff for sbQCD 4e29
	//skmeff.push_back(0.00208);  // for sbQCD sample 1e32

	xsec.push_back(6.338E7); // PYTHIA cross-section times filter pp->muX (7.923E10*0.0008)
	skmeff.push_back(1.); // for ppMuX 4e29
	//skmeff.push_back(0.268); // for ppMuX sample 1e32

	xsec.push_back(7.685E8); // PYTHIA cross-section times filter pp->eleX (7.923E10*0.0097)
	skmeff.push_back(0.558); //for ppEleX 4e29
	//skmeff.push_back(0.0413); //for ppEleX sample 1e32
	*/

	//xsec.push_back(7.923E10); // PYTHIA cross-section for MinBias (pb)
	//skmeff.push_back(0.00375);
	//xsec.push_back(1.797E3); // PYTHIA cross-section for Z -> 2mu
	//skmeff.push_back(1.);
	//xsec.push_back(1.01E5); // PYTHIA cross-section for QCD in 170 < ^pt < 230
	//skmeff.push_back(1.);
	//xsec.push_back(1.712E4); // PYTHIA cross-section for W -> e nu
	//skmeff.push_back(1.0);
	//xsec.push_back(1.787E3); // PYTHIA cross-section for Z -> 2e
	//skmeff.push_back(1.);


	/*
	xsec.push_back(5.52E10); // PYTHIA cross-section for QCD in 0 < ^pt < 15
	skmeff.push_back(1.);
	//skmeff.push_back(0.00062); // for reutuned 1e32
	xsec.push_back(1.46E9); // PYTHIA cross-section for QCD in 15 < ^pt < 20
	skmeff.push_back(1.);
	//skmeff.push_back(0.0232); // for reutuned 1e32
	*/
	/*
	xsec.push_back(6.32E8); // PYTHIA cross-section for QCD in 20 < ^pt < 30
	//skmeff.push_back(1.);
	skmeff.push_back(0.0687); // for reutuned 1e32
	xsec.push_back(1.63E8); // PYTHIA cross-section for QCD in 30 < ^pt < 50
	//skmeff.push_back(1.);
	skmeff.push_back(0.260); // for reutuned 1e32
	xsec.push_back(2.16E7); // PYTHIA cross-section for QCD in 50 < ^pt < 80
	//skmeff.push_back(1.);
	skmeff.push_back(0.701); // for reutuned 1e32
	xsec.push_back(3.08E6); // PYTHIA cross-section for QCD in 80 < ^pt < 120
	skmeff.push_back(1.);
	xsec.push_back(4.94E5); // PYTHIA cross-section for QCD in 120 < ^pt < 170
	skmeff.push_back(1.);
	xsec.push_back(1.01E5); // PYTHIA cross-section for QCD in 170 < ^pt < 230
	skmeff.push_back(1.);
	xsec.push_back(2.45E4); // PYTHIA cross-section for QCD in 230 < ^pt < 300
	skmeff.push_back(1.);
	//xsec.push_back(6.24E3); // PYTHIA cross-section for QCD in 300 < ^pt < 380
	//skmeff.push_back(1.);
	//xsec.push_back(1.78E3); // PYTHIA cross-section for QCD in 380 < ^pt < 470
	//skmeff.push_back(1.);
	//xsec.push_back(6.83E2); // PYTHIA cross-section for QCD in 470 < ^pt < 600
	//skmeff.push_back(1.);
	//xsec.push_back(2.04E2); // PYTHIA cross-section for QCD in 600 < ^pt < 800
	//skmeff.push_back(1.);
	//xsec.push_back(3.51E1); // PYTHIA cross-section for QCD in 800 < ^pt < 1000
	//skmeff.push_back(1.);

	//xsec.push_back(1.702E8); // PYTHIA cross-section for Gamma+Jet 0_15
	//skmeff.push_back(1.);
	xsec.push_back(2.570E5); // PYTHIA cross-section for Gamma+Jet 15_20
	skmeff.push_back(1.);
	xsec.push_back(1.319E5); // PYTHIA cross-section for Gamma+Jet 20_30
	skmeff.push_back(1.);
	xsec.push_back(4.114E4); // PYTHIA cross-section for Gamma+Jet 30_50
	skmeff.push_back(1.);
	xsec.push_back(7.210E3); // PYTHIA cross-section for Gamma+Jet 50_80
	skmeff.push_back(1.);
	xsec.push_back(1.307E3); // PYTHIA cross-section for Gamma+Jet 80_120
	skmeff.push_back(1.);
	xsec.push_back(2.758E2); // PYTHIA cross-section for Gamma+Jet 120_170
	skmeff.push_back(1.);

	//xsec.push_back(8.709E1); // PYTHIA cross-section for Gamma+Jet 170_300
	//skmeff.push_back(1.);
	//xsec.push_back(8.285E0); // PYTHIA cross-section for Gamma+Jet 300_500
	//skmeff.push_back(1.);
	//xsec.push_back(8.778E-1); // PYTHIA cross-section for Gamma+Jet 500_7000
	//skmeff.push_back(1.);


	xsec.push_back(1.712E4); // PYTHIA cross-section for W -> e nu
	skmeff.push_back(1.0);
	xsec.push_back(1.717E4); // PYTHIA cross-section for W -> mu nu
	skmeff.push_back(1.);
	//skmeff.push_back(0.6418);
	xsec.push_back(1.787E3); // PYTHIA cross-section for Z -> 2e
	skmeff.push_back(1.);
	xsec.push_back(1.797E3); // PYTHIA cross-section for Z -> 2mu
	//skmeff.push_back(0.4614);
	skmeff.push_back(1.);

	//xsec.push_back(1.712E4); // PYTHIA cross-section for W -> e nu
	//skmeff.push_back(1.0);
	//xsec.push_back(1.787E3); // PYTHIA cross-section for Z -> 2e
	//skmeff.push_back(1.);
	*/

	/*
	//From Silvia's email
	CSA07-Zee:
	Xsection 1787pb
	Filter on |eta|<2.5, filter eff  .648

	RelVal-Zee: same xsection as previous, no cuts on eta (filter eff 1)

	CSA07-Zmumu:
	Xsection 1797pb
	Filter on |eta|<2.5, filter eff  .4614

	CSA07-Wmunu:
	Xsection 17170.
	Filter on |eta|<2.5, filter eff  0.6418

	RelVal-Wenu:
	Xsection 17120.
	no cuts on eta (filter eff 1)
	*/

	/*
	// From Pedram's file ~bargassa/Alternate/CMSSW_1_3_1_HLT5/src/Trigger_CorrelRate.C
	sec.push_back(7.9E3); // PYTHIA cross-section for W -> e nu
	skmeff.push_back(0.941858);
	xsec.push_back(9.8E3); // PYTHIA cross-section for W -> mu nu
	skmeff.push_back(0.83052);
	xsec.push_back(8.2E2); // PYTHIA cross-section for Z -> 2e
	skmeff.push_back(0.995076);
	xsec.push_back(7.9E2); // PYTHIA cross-section for Z -> 2mu
	skmeff.push_back(0.987032);
	xsec.push_back(2.4E7); // PYTHIA cross-section for pp->muX
	skmeff.push_back(0.398551);

// Add below the x-section of any other "signal" file to be read

  
  for (unsigned int ip = 0; ip < skmeff.size(); ip++) {
    xsec[ip] *= skmeff[ip];
  }
  
  // Convert cross-sections to cm^2
  for (unsigned int i = 0; i < skmeff.size(); i++){xsec[i] *= 1.E-36;}
	*/

	// Cross sections [pb]
	/////////////////////////////////////////
	vector<Double_t> xsec;
	vector<Double_t> skmeff; // Skim efficiencies
  
  // Input Files that have HLTTree in them
	// (outputs of OpenHLT (HLTAnalyzer)
	/////////////////////////////////////////
  vector<TChain*> TabChain;
  vector<bool> doMuonCut; vector<bool> doElecCut;
  vector<TString> ProcFil;


	TString DIR_4e29 = "/afs/hep.wisc.edu/cms/rekovic/3.8e29/";

	if( sConditions.CompareTo("Validate") == 0) {

		// Ideal Conditions
		/////////////////////////////////////////////////////////////////////////////


		// Zmumu
		ProcFil.clear();
		ProcFil.push_back("/afs/hep.wisc.edu/cms/rekovic/HLTAnalyzer_From174MC/Analyzer_Zmumu.10k.root");
	  TabChain.push_back(new TChain("HltTree"));
	  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
	   	TabChain.back()->Add(ProcFil[ipfile]);
	  }
		cout << "Requires muon/electron cut: 0th sample" << endl;
	  doMuonCut.push_back(false); doElecCut.push_back(false);

		xsec.push_back(7.923E10); // PYTHIA cross-section for MinBias (pb)
		skmeff.push_back(1.0);  // skmeff for sbQCD 4e29
		
		// Zmumu
		ProcFil.clear();
		ProcFil.push_back("/afs/hep.wisc.edu/cms/rekovic/HLTAnalyzer_From174MC/Analyzer_Zee.10k.root");
	  TabChain.push_back(new TChain("HltTree"));
	  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
	   	TabChain.back()->Add(ProcFil[ipfile]);
	  }
		cout << "Requires muon/electron cut: 0th sample" << endl;
	  doMuonCut.push_back(false); doElecCut.push_back(false);

		xsec.push_back(7.923E10); // PYTHIA cross-section for MinBias (pb)
		skmeff.push_back(1.0);  // skmeff for sbQCD 4e29
		

		/*
		// ppMuX
		ProcFil.clear();
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppMuX_4e29_12kHz.root");
	  TabChain.push_back(new TChain("HltTree"));
	  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
	    TabChain.back()->Add(ProcFil[ipfile]);
	  }
	  doMuonCut.push_back(false); doElecCut.push_back(false);

		xsec.push_back(6.338E7); // PYTHIA cross-section times filter pp->muX (7.923E10*0.0008)
		skmeff.push_back(1.); // for ppMuX 4e29
	

		// ppEleX
		ProcFil.clear();
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppEleX_L1Skim_4e29_12kHz.root");
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppEleX_L1Skim_4e29_12kHz_1.root");
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppEleX_L1Skim_4e29_12kHz_2.root");
	  TabChain.push_back(new TChain("HltTree"));
	  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
	    TabChain.back()->Add(ProcFil[ipfile]);
	  }
	  doMuonCut.push_back(false); doElecCut.push_back(false);
	
		xsec.push_back(7.685E8); // PYTHIA cross-section times filter pp->eleX (7.923E10*0.0097)
		skmeff.push_back(0.558); //for ppEleX 4e29
		*/

	} // endif
	else if( sConditions.CompareTo("Ideal") == 0) {

		// Ideal Conditions
		/////////////////////////////////////////////////////////////////////////////


		// sbQCD
		ProcFil.clear();
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_QCDSingleBin_L1Skim_4e29_12kHz.root");
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_QCDSingleBin_L1Skim_4e29_12kHz_1.root");
	  TabChain.push_back(new TChain("HltTree"));
	  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
	   	TabChain.back()->Add(ProcFil[ipfile]);
	  }
		cout << "Requires muon/electron cut: 0th sample" << endl;
	  doMuonCut.push_back(true); doElecCut.push_back(true);

		xsec.push_back(7.923E10); // PYTHIA cross-section for MinBias (pb)
		skmeff.push_back(0.105);  // skmeff for sbQCD 4e29
		

		// ppMuX
		ProcFil.clear();
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppMuX_4e29_12kHz.root");
	  TabChain.push_back(new TChain("HltTree"));
	  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
	    TabChain.back()->Add(ProcFil[ipfile]);
	  }
	  doMuonCut.push_back(false); doElecCut.push_back(false);

		xsec.push_back(6.338E7); // PYTHIA cross-section times filter pp->muX (7.923E10*0.0008)
		skmeff.push_back(1.); // for ppMuX 4e29
	

		// ppEleX
		ProcFil.clear();
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppEleX_L1Skim_4e29_12kHz.root");
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppEleX_L1Skim_4e29_12kHz_1.root");
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppEleX_L1Skim_4e29_12kHz_2.root");
	  TabChain.push_back(new TChain("HltTree"));
	  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
	    TabChain.back()->Add(ProcFil[ipfile]);
	  }
	  doMuonCut.push_back(false); doElecCut.push_back(false);
	
		xsec.push_back(7.685E8); // PYTHIA cross-section times filter pp->eleX (7.923E10*0.0097)
		skmeff.push_back(0.558); //for ppEleX 4e29

	} // endif
	else if(sConditions.CompareTo("10pb") == 0) {

		// 10 pb Conditions
		////////////////////////////////////////////////////////////////////////////////


		// sbQCD
		ProcFil.clear();
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_QCDSingleBin_L1Skim_4e29_12kHz-Frontier10pbConditionOpenHLT_JetRecoRerun.root");
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_QCDSingleBin_L1Skim_4e29_12kHz-Frontier10pbConditionOpenHLT_JetRecoRerun_1.root");
	  TabChain.push_back(new TChain("HltTree"));
	  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
	    TabChain.back()->Add(ProcFil[ipfile]);
	  }
		cout << "Requires muon/electron cut: 0th sample" << endl;
	  doMuonCut.push_back(true); doElecCut.push_back(true);

		xsec.push_back(7.923E10); // PYTHIA cross-section for MinBias (pb)
		skmeff.push_back(0.105);  // skmeff for sbQCD 4e29
		
	
		// ppMuX
		ProcFil.clear();
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppMuX_4e29_12kHz-Frontier10pbConditionOpenHLT_NoEgammaCut_LowSeeds_WithJetRecoRerun.root");
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppMuX_4e29_12kHz-Frontier10pbConditionOpenHLT_NoEgammaCut_LowSeeds_WithJetRecoRerun_1.root");
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppMuX_4e29_12kHz-Frontier10pbConditionOpenHLT_NoEgammaCut_LowSeeds_WithJetRecoRerun_2.root");
	  TabChain.push_back(new TChain("HltTree"));
	  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
	    TabChain.back()->Add(ProcFil[ipfile]);
	  }
	  doMuonCut.push_back(false); doElecCut.push_back(false);

		xsec.push_back(6.338E7); // PYTHIA cross-section times filter pp->muX (7.923E10*0.0008)
		skmeff.push_back(1.); // for ppMuX 4e29
	
	
		// ppEleX
		ProcFil.clear();
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppEleX_L1Skim_4e29_12kHz-Frontier10pbConditionOpenHLT_NoEgammaCut_LowSeeds_WithJetRecoRerun.root");
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppEleX_L1Skim_4e29_12kHz-Frontier10pbConditionOpenHLT_NoEgammaCut_LowSeeds_WithJetRecoRerun_1.root");
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppEleX_L1Skim_4e29_12kHz-Frontier10pbConditionOpenHLT_NoEgammaCut_LowSeeds_WithJetRecoRerun_2.root");
		ProcFil.push_back(DIR_4e29+"HLTAnalyzer_ppEleX_L1Skim_4e29_12kHz-Frontier10pbConditionOpenHLT_NoEgammaCut_LowSeeds_WithJetRecoRerun_3.root");
	  TabChain.push_back(new TChain("HltTree"));
	  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
	    TabChain.back()->Add(ProcFil[ipfile]);
	  }
	  doMuonCut.push_back(false); doElecCut.push_back(false);
	
		xsec.push_back(7.685E8); // PYTHIA cross-section times filter pp->eleX (7.923E10*0.0097)
		skmeff.push_back(0.558); //for ppEleX 4e29


	}
	else {

		cout << "No valid Conditions specified.  Either create input files for new Conditions or use existing ones. Exiting!" << endl;
		ShowUsage();
		return 0;

	}


	/*
	ProcFil.push_back("/afs/hep.wisc.edu/cms/rekovic/3.8e29/HLTAnalyzer_PhotonJets_pt0-15_4e29.root");
	ProcFil.push_back("/afs/hep.wisc.edu/cms/rekovic/3.8e29/HLTAnalyzer_PhotonJets_pt15-20_4e29.root");
	ProcFil.push_back("/afs/hep.wisc.edu/cms/rekovic/3.8e29/HLTAnalyzer_PhotonJets_pt20-30_4e29.root");
	ProcFil.push_back("/afs/hep.wisc.edu/cms/rekovic/3.8e29/HLTAnalyzer_PhotonJets_pt30-50_4e29.root");
	ProcFil.push_back("/afs/hep.wisc.edu/cms/rekovic/3.8e29/HLTAnalyzer_PhotonJets_pt50-80_4e29.root");
	*/






	// multiply xsec by skim Eff
	cout << "Number of files (datasets) to process " << TabChain.size() << endl;

  for (unsigned int ip = 0; ip < skmeff.size(); ip++) {
    xsec[ip] *= skmeff[ip];
  }
  
  // Convert cross-sections to cm^2
  for (unsigned int i = 0; i < skmeff.size(); i++){xsec[i] *= 1.E-36;}

  ////////////////////////////////////////////////////////////

  vector<int> * iCount = new vector<int>();
  vector<int> * sPureCount = new vector<int>();
  vector<int> * pureCount = new vector<int>();
  vector<int> otmp;
  vector< vector<int> > * overlapCount = new vector< vector<int> >();
  for (int it = 0; it < Ntrig; it++){
    iCount->push_back(0);
    sPureCount->push_back(0);
    pureCount->push_back(0);
    otmp.push_back(0);
  }
  for (int it = 0; it < Ntrig; it++){
    overlapCount->push_back(otmp);
  }

  vector<Double_t> Rat,sRat,seqpRat,sseqpRat,pRat,spRat,cRat;
  vector<Double_t> Odenp;
  vector< vector<Double_t> > Onum;
  for (int it = 0; it < Ntrig; it++){
    Rat.push_back(0.);
    sRat.push_back(0.);
    seqpRat.push_back(0.);
    sseqpRat.push_back(0.);
    pRat.push_back(0.);
    spRat.push_back(0.);
    cRat.push_back(0.);
    Odenp.push_back(0.);
  }
  for (int it = 0; it < Ntrig; it++){
    Onum.push_back(Odenp);
  }
  
  for (unsigned int ip = 0; ip < TabChain.size(); ip++){
    cout<<"Available sample "<<ip<<", file " << TabChain[ip] <<endl;
    cout<<" xsec = "  << scientific << xsec[ip]/skmeff[ip]/1.E-36 << fixed << ",  skmeff = "<< skmeff[ip] <<", doMuonCut = " << doMuonCut[ip] << ", doElecCut = " << doElecCut[ip] << endl;
	}
  vector<OHltTree*> hltt;
  for (unsigned int ip = 0; ip < TabChain.size(); ip++){
    for (int it = 0; it < Ntrig; it++){
      iCount->at(it) = 0;
      sPureCount->at(it) = 0;
      pureCount->at(it) = 0;
    }
    // For binwise analysis
    vector<Double_t> Rat_bin,sRat_bin,seqpRat_bin,sseqpRat_bin,pRat_bin,spRat_bin,cRat_bin;
    for (int it = 0; it < Ntrig; it++){
      Rat_bin.push_back(0.);
      sRat_bin.push_back(0.);
      seqpRat_bin.push_back(0.);
      sseqpRat_bin.push_back(0.);
      pRat_bin.push_back(0.);
      spRat_bin.push_back(0.);
      cRat_bin.push_back(0.);
    }

    hltt.push_back(new OHltTree((TTree*)TabChain[ip],Ntrig));

    int deno = NEntries; 
		int chainEntries = (int)hltt[ip]->fChain->GetEntries(); 
    if (NEntries <= 0 || NEntries > chainEntries) {
      deno = chainEntries;
    }
		cout<<"---------------------------------------------------------------" << endl;
    cout<<"Processing bin "<<ip<<" ( "<< deno <<" events ) "<<", file " << TabChain[ip] <<" (has "<<hltt[ip]->fChain->GetEntries()<<" events ) "<<endl;
		cout<<scientific;
		cout.precision(5);
    cout<<" xsec = "  << xsec[ip]/skmeff[ip]/1.E-36 << fixed << ",  skmeff = "<< skmeff[ip] <<", doMuonCut = " << doMuonCut[ip] << ", doElecCut = " << doElecCut[ip] << endl;
		cout<<"---------------------------------------------------------------" << endl;

    hltt[ip]->Loop(iCount,sPureCount,pureCount,overlapCount,trignames,map_TrigPrescls,deno,doMuonCut[ip],doElecCut[ip]);
    double mu = bunchCrossingTime * xsec[ip] * ILumi * maxFilledBunches / nFilledBunches;
    for (int it = 0; it < Ntrig; it++){
      // Get global overlaps
      for (int jt = 0; jt != Ntrig; ++jt){
	if (jt==it){
	  (Onum.at(it))[jt] = (((double)(iCount->at(it)) * xsec[ip])); 
	} else {
	  (Onum.at(it))[jt] += ( (double)(overlapCount->at(it).at(jt)) * xsec[ip]);     
	}
      }
      Odenp[it] += ((double)(iCount->at(it)) * xsec[ip]); // ovelap denominator
      
      Rat[it] += collisionRate*(1. - exp(- mu * eff(iCount->at(it),deno)));  // Single rates
      sRat[it] += pow(collisionRate*mu * seff(iCount->at(it),deno),2.);    //
      seqpRat[it] += collisionRate*(1. - exp(- mu * eff(sPureCount->at(it),deno)));  // Single rates
			if(it<4) {
				
				cout << "i=" << it << " Rate=" << Rat[it] << " +/- " << sqrt(sRat[it]) << ", passed evts=" << iCount->at(it) << ", total evts= " << deno << endl;

			}
      sseqpRat[it] += pow(collisionRate*mu * seff(sPureCount->at(it),deno),2.);    //
      pRat[it] += collisionRate*(1. - exp(- mu * eff(pureCount->at(it),deno)));  // Single rates
      spRat[it] += pow(collisionRate*mu * seff(pureCount->at(it),deno),2.);    //

      // Binwise
      Rat_bin[it] += collisionRate*(1. - exp(- mu * eff(iCount->at(it),deno)));  // Single rates
      sRat_bin[it] += pow(collisionRate*mu * seff(iCount->at(it),deno),2.);    //
      seqpRat_bin[it] += collisionRate*(1. - exp(- mu * eff(sPureCount->at(it),deno)));  // Single rates
      sseqpRat_bin[it] += pow(collisionRate*mu * seff(sPureCount->at(it),deno),2.);    //
      pRat_bin[it] += collisionRate*(1. - exp(- mu * eff(pureCount->at(it),deno)));  // Single rates
      spRat_bin[it] += pow(collisionRate*mu * seff(pureCount->at(it),deno),2.);    //

    }

    // Print binwise rates:
    // Loop over triggers
    Double_t RTOT_bin = 0.;
    Double_t sRTOT_bin = 0.;
    Double_t curat_bin = 0.;
    for (int it = 0; it < Ntrig; it++){
      curat_bin += seqpRat_bin[it];
      cRat_bin[it] = curat_bin;
      RTOT_bin += seqpRat_bin[it];                                            // Total Rate
      sRTOT_bin += sseqpRat_bin[it];
    }
    sRTOT_bin = sqrt(sRTOT_bin);

    // Print binwise
    cout.setf(ios::floatfield,ios::fixed);
    cout<<setprecision(3);
    for (int it=0; it < Ntrig; it++){
      cout  << setw(3) << it << ")" << setw(30) << trignames[it]  << " (" << setw(8) << map_TrigPrescls.find(trignames[it])->second << ")"
	    << " :   Indiv.: " << setw(8) << Rat_bin[it] << " +/- " << setw(5) << sqrt(sRat_bin[it]) 
	    << "   seqPure: " << setw(8) << seqpRat_bin[it]
	    << "   Pure: " << setw(8) << pRat_bin[it] 
	    << "   Cumul: " << setw(8) << cRat_bin[it] << "\n"<<flush;
    }
    cout << "\n"<<flush;
    cout << setw(60) << "TOTAL RATE : " << setw(5) << RTOT_bin << " +- " << sRTOT_bin << " Hz" << "\n";
    cout << "\n"<<flush;
    
  } // end for TabChain.size()

  //exit(0);


  // Loop over triggers
  Double_t RTOT = 0.;
  Double_t sRTOT = 0.;
  Double_t curat = 0.;
  for (int it = 0; it < Ntrig; it++){
    curat += seqpRat[it];
    cRat[it] = curat;
    RTOT += seqpRat[it];                                            // Total Rate
    sRTOT += sseqpRat[it];
  }


  sRTOT = sqrt(sRTOT);
    
  ////////////////////////////////////////////////////////////
  // Results
  cout<<setprecision(3);

  cout << endl;
  cout << "Trigger global overlaps : " << endl;
  for (int it = 0; it != Ntrig; ++it){
    for (int jt = 0; jt != Ntrig; ++jt){
      if (jt>=it) {
	// Overlap O(ij) = T(i) x T(j) / T(j)
	//cout << "i=" << it << " j=" << jt << "     " << eff((Onum.at(it))[jt],Odenp[jt]) << endl;   
      }
    }
  }

  cout.setf(ios::floatfield,ios::fixed);
  cout<<setprecision(1);

  cout << "\n";
  cout << "Trigger Rates [Hz] : " << "\n";
  cout << "------------------------------------------------------------------------------------------------------------------\n";
  // This is with the accurate formula: 
  for (int it=0; it < Ntrig; it++){
    cout  << setw(3) << it << ")" << setw(30) << trignames[it]  << " (" << setw(5) << map_TrigPrescls.find(trignames[it])->second << ")" 
	  << " :   Indiv.: " << setw(8) << Rat[it] << " +/- " << setw(8) << sqrt(sRat[it]) 
	  << "   seqPure: " << setw(8) << seqpRat[it]
	  << "   Pure: " << setw(8) << pRat[it] 
	  << "   Cumul: " << setw(8) << cRat[it] << "\n"<<flush;
  }
  cout << "\n"<<flush;
  cout << setw(60) << "TOTAL RATE : " << setw(5) << RTOT << " +- " << sRTOT << " Hz" << "\n";
  cout << "------------------------------------------------------------------------------------------------------------------\n"<<flush;

	char sLumi[10];
	sprintf(sLumi,"%1.1e",ILumi);
	TString hltTableFileName= TString("hltTable_") + TString(sLumi) + TString("_") + sConditions + TString("Conditions") + sVersion;
	TString texFile = hltTableFileName + TString(".tex");
	TString dviFile = hltTableFileName + TString(".dvi");
	TString psFile  = hltTableFileName + TString(".ps");
  ofstream outFile(texFile.Data());
  if (!outFile){cout<<"Error opening output file"<< endl;}
  outFile <<setprecision(1);
  outFile.setf(ios::floatfield,ios::fixed);
  outFile << "\\documentclass[amsmath,amssymb]{revtex4}" << endl;
  outFile << "\\usepackage{longtable}" << endl;
  outFile << "\\usepackage{color}" << endl;
  outFile << "\\begin{document}" << endl;
  outFile << "\\newcommand{\\met}{\\ensuremath{E\\kern-0.6em\\lower-.1ex\\hbox{\\/}\\_T}}" << endl;


    
		outFile << "\\begin{footnotesize}" << endl;
    outFile << "\\begin{longtable}{|c|l|c|c|c|c|c|}" << endl;
		outFile << "\\caption[Cuts]{New paths are introduced in addition to standard '1e32' paths.  Detailed description of the newly introduced paths is given at the end of the table.  Prescale listed is the total prescale (L1xHLT).  CSA07 data are used and HLTrigger table of CMSSW 16X.  Available HLT bandwith is 150 Hz = ((1 GB/s / fact. 3) - 100 MB/s for ALCA) / 1.5 MB/event.  L1 is fixed to 12 kHz. } \\label{CUTS} \\\\ " << endl;


    outFile << "\\hline \\multicolumn{7}{|c|}{\\bf \\boldmath HLT for L = "<< sLumi  << "}\\\\  \\hline" << endl;
    outFile << "{\\bf Status} & " << endl;
    outFile << "{\\bf Path Name} & " << endl;
    outFile << "{\\bf L1 condtition} & " << endl;
    outFile << "\\begin{tabular}{c} {\\bf Thresholds} \\\\ {\\bf $[$GeV$]$} \\end{tabular} & " << endl;
    outFile << "\\begin{tabular}{c} {\\bf Total} \\\\ {\\bf Prescale} \\end{tabular} & " << endl;
    outFile << "\\begin{tabular}{c} {\\bf HLT Rate} \\\\ {\\bf $[$Hz$]$} \\end{tabular} &" << endl;
    outFile << "\\begin{tabular}{c} {\\bf Total Rate} \\\\ {\\bf $[$Hz$]$} \\end{tabular} \\\\ \\hline" << endl;
		outFile << "\\endfirsthead " << endl;

    outFile << "\\multicolumn{7}{r}{\\bf \\bfseries --continued from previous page (L = " << sLumi << ")"  << "}\\\\ \\hline " << endl;
    outFile << "{\\bf Status} & " << endl;
    outFile << "{\\bf Path Name} & " << endl;
    outFile << "{\\bf L1 condtition} & " << endl;
    outFile << "\\begin{tabular}{c} {\\bf Thresholds} \\\\ {\\bf $[$GeV$]$} \\end{tabular} & " << endl;
    outFile << "\\begin{tabular}{c} {\\bf Total} \\\\ {\\bf Prescale} \\end{tabular} & " << endl;
    outFile << "\\begin{tabular}{c} {\\bf HLT Rate} \\\\ {\\bf $[$Hz$]$} \\end{tabular} &" << endl;
    outFile << "\\begin{tabular}{c} {\\bf Total Rate} \\\\ {\\bf $[$Hz$]$} \\end{tabular} \\\\ \\hline" << endl;
		outFile << "\\endhead " << endl;

		outFile << "\\hline \\multicolumn{6}{|r|}{{Continued on next page}} \\\\ \\hline " << endl;
		outFile << "\\endfoot " << endl;

		outFile << "\\hline " << endl;
		outFile << "\\endlastfoot " << endl;

  	for (int it=0; it < Ntrig; it++){

			TString tempTrigName = trignames[it];
			TString tempL1BitName = menu->GetHltL1BitMap().find(trignames[it])->second;
			TString tempThreshold = menu->GetHltThresholdMap().find(trignames[it])->second;

    	if(tempTrigName.Contains("Apt")) {

				outFile << "new & " ;
				tempTrigName.ReplaceAll("AptHLT","");

			}
			else {

				outFile << " 1e32 & " ;

			}

			tempTrigName.ReplaceAll("_","\\_");
			tempL1BitName.ReplaceAll("_","\\_");
			//tempTrigName.ReplaceAll("HLT","");
			outFile << "\\color{blue}"  << tempTrigName << " & " << "${\\it " << tempL1BitName << "}$ "<< " & " << tempThreshold << " & " <<  map_TrigPrescls.find(trignames[it])->second  << " & " << Rat[it] << " {$\\pm$ " << sqrt(sRat[it]) << "} & " << cRat[it] << "\\\\" << endl;
  	}

    outFile << "\\hline \\multicolumn{6}{|c|}{\\bf \\boldmath Total HLT rate (Hz) } & "<<  RTOT << " {$\\pm$ " << sRTOT << "} \\\\  \\hline" << endl;
    outFile << "\\hline " << endl;
		outFile << "\\multicolumn{7}{|l|}{ ($\\star$): {\\it Mu3\\_IsoEG5, Mu5\\_IsoEG10, Mu3\\_IsoEG12 }} \\\\ \\hline " << endl;
		outFile << "\\multicolumn{7}{|l|}{ ($\\dagger$): {\\it SingleJet150, DoubleJet70, TripleJet50 }} \\\\ \\hline " << endl;
		outFile << "\\multicolumn{7}{|l|}{ ($\\ddagger$): {\\it SingleJet150, DoubleJet70, TripleJet50, QuadJet30 }} \\\\ \\hline " << endl;
		outFile << "\\multicolumn{7}{|l|}{ ($\\Diamond$): {\\it SingleJet150, DoubleJet70, TripleJet50, QuadJet30, HTT300 }} \\\\ \\hline " << endl;
    outFile << "\\hline " << endl;
    outFile << "\\hline " << endl;
		outFile << "\\multicolumn{7}{|c|}{ {\\it \\bf HLT Requirements for new introduced paths (in GeV) }  } \\\\ \\hline " << endl;
		outFile << "\\multicolumn{7}{|l|}{ 1MuonA: $|\\eta|<2.5$, $L2Pt+3.9Err<$A, $L3Pt+2.2Err<$A  } \\\\ \\hline " << endl;
		outFile << "\\multicolumn{7}{|l|}{ 1jetA: $recoJetCalPt<$A  } \\\\ \\hline " << endl;
		outFile << "\\multicolumn{7}{|l|}{ 1ElectronA\\_DefThr: $Et<$A, $HCAL<3$, $E/p<1.5(2.45)$, $TrkIso<0.06$} \\\\ " << endl;
		outFile << "\\multicolumn{7}{|l|}{ 1ElectronA\\_DblThr: $Et<$A, $HCAL<6$, $E/p<3.0(4.90)$, $TrkIso<0.12$} \\\\  " << endl;
		outFile << "\\multicolumn{7}{|l|}{ 1ElectronA\\_Open: $Et<$A} \\\\ \\hline " << endl;
		outFile << "\\multicolumn{7}{|l|}{ 1PhotonA\\_DefThr: $Et<$A, $ECAL<1.5$, $HCAL<6(4)$, $TrkIso=0$} \\\\  " << endl;
		outFile << "\\multicolumn{7}{|l|}{ 1PhotonA\\_DblThr: $Et<$A, $ECAL<3.0$, $HCAL<12(8)$, $TrkIso\\leq2$} \\\\  " << endl;
		outFile << "\\multicolumn{7}{|l|}{ 1PhotonA\\_Open: $Et<$A, $TrkIso=0$} \\\\ " << endl;
		outFile << "\\multicolumn{7}{|l|}{ 1PhotonANoTrkIso: $Et<$A} \\\\ \\hline" << endl;

    //outFile << "\\end{tabular}" << endl;
    outFile << "\\end{longtable}" << endl;
    outFile << "\\end{footnotesize}" << endl;
    outFile << "\\clearpage" << endl;

  outFile << "\\end{document}";
  outFile.close();


	TString Command = TString("latex ") + texFile + TString("; dvips -t A4 -f ") + dviFile + TString(" -o ") + psFile;

	cout << "Executing the following latex command: " << endl;
	cout << Command << endl;
	// do it again to fix for column size shift within the table
	//for(int i=0;i<2;i++) system("latex hltTable_2.0e+31.tex; dvips -t A4 -f hltTable_2.0e+31.dvi -o hltTable_2.0e+31.ps");
	for(int i=0;i<2;i++) system(Command.Data());

    //gSystem->Exec("latex hltTable.tex");
    //gSystem->Exec("dvips -t letter -f hltTable.dvi -o hltTable.ps");
    //gSystem->Exec("dvips -t A4 -f hltTable.dvi -o hltTable.ps");


}


void BookOHltMenu_2e31_v03(OHltMenu*  menu_2e31_v03, double &iLumi, double &nBunches) {

		iLumi = 2E31;
		nBunches = 156;

		menu_2e31_v03->AddHlt("AptHLT1MuonLevel1Open","SingleMuOpen",1.5E3,"-"," "); //1000; 
		menu_2e31_v03->AddHlt("AptHLT1MuonLevel1","SingleMu7,DoubleMu3",2E2,"-"," "); //1000; 
		//menu_2e31_v03->AddHlt("AptHLT1MuonIso3_DefThr","SingleMu3",1,"3"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1MuonIso5_DefThr","SingleMu5",1,"5"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1MuonIso7_DefThr","SingleMu7",1,"7"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1MuonIso9_DefThr","SingleMu7",1,"9"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1MuonIso11_DefThr","SingleMu7",1,"11"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Muon3_DefThr","SingleMu3",1,"3"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Muon5_DefThr","SingleMu5",1,"5"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Muon7_DefThr","SingleMu7",1,"7"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Muon9_DefThr","SingleMu7",1,"9"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Muon11_DefThr","SingleMu7",1,"11"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Muon13_DefThr","SingleMu7",1,"13"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Muon15_DefThr","SingleMu7",1,"15"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1MuonIso3","SingleMu3",1E6,"3"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1MuonIso5","SingleMu5",1E6,"5"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1MuonIso7","SingleMu7",1E6,"7"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1MuonIso9","SingleMu7",1,"9"," "); //5E2; 
		//menu_2e31_v03->AddHlt("AptHLT1MuonIso11","SingleMu7",1,"11"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Muon3","SingleMu3",5E2,"3"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Muon5","SingleMu5",5E2,"5"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Muon7","SingleMu7",1E2,"7"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Muon9","SingleMu7",5E1,"9"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Muon11","SingleMu7",1,"11"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Muon13","SingleMu7",1,"13"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Muon15","SingleMu7",1,"15"," "); 
		menu_2e31_v03->AddHlt("AptHLT1L2Muon11","SingleMu7",1E1,"11"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1L2Muon16","SingleMu7",1E1,"16"," "); 
		menu_2e31_v03->AddHlt("AptHLT2Muon3","DobuleMu3",1,"(3,3)"," "); 
		menu_2e31_v03->AddHlt("HLT1MuonIso","SingleMu7",1,"11"," "); 
		menu_2e31_v03->AddHlt("HLT1MuonNonIso","SingleMu7",1,"16"," "); 
		menu_2e31_v03->AddHlt("HLT2MuonIso","DoubleMu3",1,"(3,3)"," "); 
		menu_2e31_v03->AddHlt("HLT2MuonNonIso","DoubleMu3",1,"(3,3)"," "); 
		menu_2e31_v03->AddHlt("HLT2MuonJPsi","DoubleMu3",1,"(3,3)"," "); 
		menu_2e31_v03->AddHlt("HLT2MuonUpsilon","DoubleMu3",1,"(3,3)"," "); 
		menu_2e31_v03->AddHlt("HLT2MuonZ","DoubleMu3",1,"(7,7)"," "); 
		menu_2e31_v03->AddHlt("HLTNMuonNonIso","DoubleMu3",1,"(3,3,3) "," "); 
		menu_2e31_v03->AddHlt("HLT2MuonSameSign","DoubleMu3",1,"(3,3)"," "); 
		//menu_2e31_v03->AddHlt("HLT1MuonPrescalePt3","SingleMu3",1E6,"3"," "); 
		//menu_2e31_v03->AddHlt("HLT1MuonPrescalePt5","SingleMu5",1E6,"5"," "); 
		//menu_2e31_v03->AddHlt("HLT1MuonPrescalePt7x7","SingleMu7",1E6,"7"," "); 
		//menu_2e31_v03->AddHlt("HLT1MuonPrescalePt7x10","SingleMu7",1E6,"7,10"," "); 
		menu_2e31_v03->AddHlt("HLTB1JetMu","Mu5_Jet15",1,"20 "," "); 
		menu_2e31_v03->AddHlt("HLTB2JetMu","Mu5_Jet15",1,"120"," "); 
		menu_2e31_v03->AddHlt("HLTB3JetMu","Mu5_Jet15",1,"70"," "); 
		menu_2e31_v03->AddHlt("HLTB4JetMu","Mu5_Jet15",1,"40"," "); 
		menu_2e31_v03->AddHlt("HLTBHTMu","HTT250",10,"300"," "); 
		menu_2e31_v03->AddHlt("HLTBJPsiMuMu","DoubleMu3",1,"(4,4)"," "); 
		menu_2e31_v03->AddHlt("HLTXMuonBJet","Mu5_Jet15",1,"(7,35)"," "); 
		menu_2e31_v03->AddHlt("HLTXMuonBJetSoftMuon","Mu5_Jet15",1,"(7,20)"," "); 
		menu_2e31_v03->AddHlt("HLTXMuonJets","Mu5_Jet15",1,"(7,40)"," "); 
		menu_2e31_v03->AddHlt("HLTXElectronMuon","\\star",1,"(8,7)"," "); 
		menu_2e31_v03->AddHlt("HLTXElectronMuonRelaxed","\\star",1,"(10,10)"," "); 
		menu_2e31_v03->AddHlt("HLTXMuonTau","Mu5_TauJet20",1,"(15,20)"," "); 

		menu_2e31_v03->AddHlt("AptHLT1Level1jet15","SingleJet15",1E5,"-"," "); 
		menu_2e31_v03->AddHlt("AptHLT1jet30","SingleJet15",1E4,"30"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1jet40","SingleJet15",1E3,"40"," "); 
		menu_2e31_v03->AddHlt("AptHLT1jet50","SingleJet30",5E2,"50"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1jet60","SingleJet30",1E2,"60"," "); 
		menu_2e31_v03->AddHlt("AptHLT1jet80","SingleJet50",1E2,"80"," "); 
		menu_2e31_v03->AddHlt("AptHLT1jet110","SingleJet70",2E1,"110"," "); 
		menu_2e31_v03->AddHlt("AptHLT1jet180","SingleJet70",1,"180"," "); 
		menu_2e31_v03->AddHlt("AptHLT1jet250","SingleJet70",1,"250"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1jet150","SingleJet100",5,"150"," "); 
		//menu_2e31_v03->AddHlt("AptHLT2jet100","SingleJet100, DoubleJet70",1,"(100,100)"," "); 
		menu_2e31_v03->AddHlt("HLT1jet","SingleJet150",1,"200"," "); 
		menu_2e31_v03->AddHlt("HLT2jet","SingleJet150, DoubleJet70",1,"150"," "); 
		menu_2e31_v03->AddHlt("HLT3jet","\\dagger",1,"85"," "); 
		menu_2e31_v03->AddHlt("HLT4jet","\\ddagger",1,"60"," "); 
		menu_2e31_v03->AddHlt("HLT1MET","ETM40",1,"65"," "); 
		menu_2e31_v03->AddHlt("HLT2jetAco","SingleJet150, DoubleJet70",1,"125"," "); 
		menu_2e31_v03->AddHlt("HLT1jet1METAco","SingleJet150",1,"(100,60)"," "); 
		menu_2e31_v03->AddHlt("HLT1jet1MET","SingleJet150",1,"(180,60)"," "); 
		menu_2e31_v03->AddHlt("HLT2jet1MET","SingleJet150",1,"(125,60)"," "); 
		menu_2e31_v03->AddHlt("HLT3jet1MET","SingleJet150",1,"(60,60)"," "); 
		menu_2e31_v03->AddHlt("HLT4jet1MET","SingleJet150",1,"(35,60)"," "); 
		menu_2e31_v03->AddHlt("HLT1MET1HT","HTT300 ",1,"(350,65)"," "); 
		menu_2e31_v03->AddHlt("HLT1SumET","ETT60",1E1,"120"," "); 
		//menu_2e31_v03->AddHlt("HLT1jetPE7","SingleJet15",1,"30"," "); 
		//menu_2e31_v03->AddHlt("HLT1jetPE5","SingleJet30",1,"60"," "); //10000; 
		//menu_2e31_v03->AddHlt("HLT1jetPE3","SingleJet70 ",1,"110"," "); //100; 
		//menu_2e31_v03->AddHlt("HLT1jetPE1","SingleJet100 ",1,"150"," "); 
		//menu_2e31_v03->AddHlt("HLT1METPre3","MinBians_HTT10",1,"15"," "); 
		//menu_2e31_v03->AddHlt("HLT1METPre2","MinBians_HTT10",1,"20"," "); 
		//menu_2e31_v03->AddHlt("HLT1METPre1","ETM20",1E3,"30"," "); //100; 
		menu_2e31_v03->AddHlt("AptHLT1MET20","ETM20",1E4,"20"," "); //100; 
		menu_2e31_v03->AddHlt("AptHLT1MET35","ETM30",1E2,"35"," "); //100; 
		menu_2e31_v03->AddHlt("AptHLT1MET50","ETM40",1E1,"50"," "); //100; 
		menu_2e31_v03->AddHlt("AptHLT1MET65","ETM40",1,"65"," "); //100; 
		menu_2e31_v03->AddHlt("AptHLT1MET75","ETM40",1,"75"," "); //100; 
		//menu_2e31_v03->AddHlt("HLT2jetAve30","SingleJet15",5E2,"30"," "); 
		//menu_2e31_v03->AddHlt("HLT2jetAve60","SingleJet30",1E2,"60"," "); 
		//menu_2e31_v03->AddHlt("HLT2jetAve110","SingleJet70",1E1,"110"," "); 
		//menu_2e31_v03->AddHlt("HLT2jetAve150","SingleJet100",1,"150"," "); 
		//menu_2e31_v03->AddHlt("HLT2jetAve200","SingleJet150",1,"200"," "); 
		//menu_2e31_v03->AddHlt("AptHLT2jetAve15","SingleJet15",5E3,"15"," "); 
		//menu_2e31_v03->AddHlt("AptHLT2jetAve30","SingleJet30",5E2,"30"," "); 
		//menu_2e31_v03->AddHlt("AptHLT2jetAve50","SingleJet50",1E2,"50"," "); 
		//menu_2e31_v03->AddHlt("AptHLT2jetAve70","SingleJet70",1E1,"70"," "); 
		//menu_2e31_v03->AddHlt("AptHLT2jetAve130","SingleJet130",5,"130"," "); 
		//menu_2e31_v03->AddHlt("AptHLT2jetAve220","SingleJet220",1,"220"," "); 
		menu_2e31_v03->AddHlt("HLT2jetvbfMET","ETM30",10,"(40,60)"," "); 
		menu_2e31_v03->AddHlt("HLTS2jet1METNV","SingleJet150",1,"(-,60)"," "); 
		menu_2e31_v03->AddHlt("HLTS2jet1METAco","SingleJet150",1,"(-,70)"," "); 
		menu_2e31_v03->AddHlt("HLTSjet1MET1Aco","SingleJet150",1,"(-,70)"," "); 
		menu_2e31_v03->AddHlt("HLTSjet2MET1Aco","SingleJet150",1,"(-,70)"," "); 
		menu_2e31_v03->AddHlt("HLTS2jetAco","SingleJet150",1,"(-,-)"," "); 
		menu_2e31_v03->AddHlt("HLTJetMETRapidityGap","IsoEG10_Jet20_ForJet10",1E2,"20"," "); //100; 
		menu_2e31_v03->AddHlt("AptHLT4jet30","QuadJet15",1E2,"30"," "); 
		menu_2e31_v03->AddHlt("AptHLTXRelMu3_3jet30","Mu3_TripleJet20",1E1,"(3,30)"," "); 
		menu_2e31_v03->AddHlt("AptHLTXRelEl5_3jet30","EG5_TripleJet20",1E1,"(5,30)"," "); 

		//menu_2e31_v03->AddHlt("AptHLT1Electron5_Open","SingleEG5",1E2,"5"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Electron7_Open","SingleEG5",1E4,"7"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Electron9_Open","SingleEG5",1E4,"9"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Electron11_Open","SingleEG8",1E4,"11"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Electron13_Open","SingleEG8",1E4,"13"," "); 

		menu_2e31_v03->AddHlt("AptHLT1Electron15_Open","SingleEG10",1E1,"15"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Electron20_Open","SingleEG15",1E1,"20"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Electron25_Open","SingleEG15",1E1,"25"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Electron30_Open","SingleEG20",1,"30"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1ElectronIso8_10_DefThr","SingleIsoEG8",1E1,"10"," "); 
		menu_2e31_v03->AddHlt("AptHLT1ElectronIso8_12_DefThr","SingleIsoEG8",1,"12"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Electron12_15_DefThr","SingleEG12",1,"15"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1ElectronIso8_10_DblThr","SingleIsoEG8",1E1,"10"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1ElectronIso8_12_DblThr","SingleIsoEG8",1,"12"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1ElectronIso8_15_DblThr","SingleIsoEG8",1,"15"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Electron12_15_DblThr","SingleEG12",1,"15"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Electron12_17_DblThr","SingleEG12",1,"17"," "); 
		menu_2e31_v03->AddHlt("AptHLT2ElectronIso8_10_Open","DoubleIsoEG8",1,"10"," "); 
		menu_2e31_v03->AddHlt("AptHLT2Electron10_12_Open","DoubleEG10",1,"12"," "); 
		menu_2e31_v03->AddHlt("HLT1Electron","SingleIsoEG12",1,"15"," "); 
		menu_2e31_v03->AddHlt("HLT1ElectronRelaxed","SingleEG15",1,"18"," "); 
		menu_2e31_v03->AddHlt("HLT2Electron","DoubleIsoEG8",1,"10"," "); 
		menu_2e31_v03->AddHlt("HLT2ElectronRelaxed","DoubleEG10",1,"12"," "); 

		//menu_2e31_v03->AddHlt("AptHLT1Photon5NoTrkIso","SingleEG5",5E3,"5"," ");
		menu_2e31_v03->AddHlt("AptHLT1Photon10NoTrkIso","SingleEG8",1E3,"10"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon15NoTrkIso","SingleEG10",1E2,"15"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Photon20NoTrkIso","SingleEG10",1E2,"20"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon25NoTrkIso","SingleEG10",1E2,"25"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Photon30NoTrkIso","SingleEG20",1E2,"30"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon35NoTrkIso","SingleEG20",1E2,"35"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Photon40NoTrkIso","SingleEG25",2,"40"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon5_Open","SingleEG5",10,"5"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon10_Open","SingleEG8",1,"10"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon12NoTrkIso","SingleEG8",1,"12"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon12_Open","SingleEG8",1,"12"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon15"_Open,"SingleEG10",1,"15"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon25"_Open,"SingleEG10",1,"25"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon30"_Open,"SingleEG20",1,"30"," "); 
		//menu_2e31_v03->AddHlt("AptHLT2Photon15"_Open,"DoubleEG5",1,"15"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1PhotonIso8_15_DefThr","SingleIsoEG8",1,"15"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1PhotonIso8_17_DefThr","SingleIsoEG8",1,"17"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1PhotonIso8_20_DefThr","SingleIsoEG8",1,"20"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1PhotonIso8_25_DefThr","SingleIsoEG8",1,"25"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1PhotonIso8_30_DefThr","SingleIsoEG8",1,"30"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon12_20_DefThr","SingleEG12",1,"20"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon12_30_DefThr","SingleEG12",1,"30"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon12_35_DefThr","SingleEG12",1,"35"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon12_40_DefThr","SingleEG12",1,"40"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1PhotonIso8_20_DblThr","SingleIsoEG8",1,"20"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1PhotonIso8_25_DblThr","SingleIsoEG8",1,"25"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1PhotonIso8_30_DblThr","SingleIsoEG8",1,"30"," "); 
		//menu_2e31_v03->AddHlt("AptHLT1Photon12_15_DblThr","SingleEG12",5,"15"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Photon12_25_DefThr","SingleEG12",1,"25"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Photon12_25_DblThr","SingleEG12",1,"25"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Photon12_30_DblThr","SingleEG12",1,"30"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Photon12_35_DblThr","SingleEG12",1,"35"," "); 
		menu_2e31_v03->AddHlt("AptHLT1Photon12_40_DblThr","SingleEG12",1,"40"," "); 
		menu_2e31_v03->AddHlt("AptHLT2PhotonIso8_15_Open","DoubleIsoEG8",1,"15"," "); 
		menu_2e31_v03->AddHlt("AptHLT2Photon10_20_Open","DoubleEG10",1,"20"," "); 
		menu_2e31_v03->AddHlt("AptHLT2Photon10_20NoTrkIso","DoubleEG10",1,"(20,20)"," "); 
		//menu_2e31_v03->AddHlt("AptHLT2Photon10_30NoTrkIso","DoubleEG10",1,"(30,30)"," "); 
		menu_2e31_v03->AddHlt("HLT1Photon","SingleIsoEG12",1,"30"," "); 
		menu_2e31_v03->AddHlt("HLT1PhotonRelaxed","SingleEG15",1,"40"," "); 
		menu_2e31_v03->AddHlt("HLT2Photon","DoubleIsoEG8",1,"(20,20)"," "); 
		menu_2e31_v03->AddHlt("HLT2PhotonRelaxed","DoubleEG10",1,"(20,20)"," "); 
		menu_2e31_v03->AddHlt("HLT1EMHighEt","SingleEG15",1,"80"," "); 
		menu_2e31_v03->AddHlt("HLT1EMVeryHighEt","SingleEG15",1,"200"," "); 
		menu_2e31_v03->AddHlt("HLT2ElectronZCounter","DoubleIsoEG8",1,"(10,10)"," "); 
		menu_2e31_v03->AddHlt("HLT2ElectronExclusive","ExclusiveDoubleIsoEG6",1,"(6,6)"," "); 
		menu_2e31_v03->AddHlt("HLT2PhotonExclusive","ExclusiveDoubleIsoEG6",1,"(10,10)"," "); 
		menu_2e31_v03->AddHlt("HLT1PhotonL1Isolated","SingleIsoEG10",1E2,"12"," "); 

		menu_2e31_v03->AddHlt("HLTB1Jet","\\Diamond ",1,"180"," "); 
		menu_2e31_v03->AddHlt("HLTB2Jet","\\Diamond ",1,"120"," "); 
		menu_2e31_v03->AddHlt("HLTB3Jet","\\Diamond ",1,"70"," "); 
		menu_2e31_v03->AddHlt("HLTB4Jet","\\Diamond ",1,"40"," "); 
		menu_2e31_v03->AddHlt("HLTBHT","\\Diamond ",1,"470"," "); 
		menu_2e31_v03->AddHlt("HLT1Tau","SingleTauJet80",1,"15"," "); 
		menu_2e31_v03->AddHlt("HLT1Tau1MET","TauJet30_ETM30",1,"15"," "); 
		menu_2e31_v03->AddHlt("HLT2TauPixel","TauJet40",1,"15"," "); 
		menu_2e31_v03->AddHlt("HLTXElectronBJet","IsoEG10_Jet20",1,"(10,35)"," "); 
		menu_2e31_v03->AddHlt("HLTXElectron1Jet","IsoEG10_Jet30",1,"(12,40)"," "); 
		menu_2e31_v03->AddHlt("HLTXElectron2Jet","IsoEG10_Jet30",1,"(12,80)"," "); 
		menu_2e31_v03->AddHlt("HLTXElectron3Jet","IsoEG10_Jet30",1,"(12,60)"," "); 
		menu_2e31_v03->AddHlt("HLTXElectron4Jet","IsoEG10_Jet30",1,"(12,35)"," "); 
		menu_2e31_v03->AddHlt("HLTXElectronTau","IsoEG10_TauJet20",1,"(12,20)"," "); 

		//menu_2e31_v03->AddHlt("HLTHcalIsolatedTrack","NA",1,"NA"," "); 
		//menu_2e31_v03->AddHlt("HLTHcalIsolatedTrackNoEcalIsol"," ",1," "," ");  
		menu_2e31_v03->AddHlt("HLTMinBiasPixel","ZeroBias",1,"- "," "); 
		//menu_2e31_v03->AddHlt("HLTMinBiasForAlignment"," ",1," "," "); 
		menu_2e31_v03->AddHlt("HLTMinBias","MinBias_HTT10",1,"-"," "); 
		menu_2e31_v03->AddHlt("HLTZeroBias","ZeroBias",1,"-"," "); 

}

void BookOHltMenu_2e31_v02(OHltMenu*  menu_2e31_v2, double &iLumi, double &nBunches) {

		iLumi = 2E31;
		nBunches = 156;

		menu_2e31_v2->AddHlt("AptHLT1MuonLevel1Open","SingleMuOpen",1.5E3,"-"," "); //1000; 
		menu_2e31_v2->AddHlt("AptHLT1MuonLevel1","SingleMu7,DoubleMu3",2E2,"-"," "); //1000; 
		//menu_2e31_v2->AddHlt("AptHLT1MuonIso3_DefThr","SingleMu3",1,"3"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1MuonIso5_DefThr","SingleMu5",1,"5"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1MuonIso7_DefThr","SingleMu7",1,"7"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1MuonIso9_DefThr","SingleMu7",1,"9"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1MuonIso11_DefThr","SingleMu7",1,"11"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Muon3_DefThr","SingleMu3",1,"3"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Muon5_DefThr","SingleMu5",1,"5"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Muon7_DefThr","SingleMu7",1,"7"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Muon9_DefThr","SingleMu7",1,"9"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Muon11_DefThr","SingleMu7",1,"11"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Muon13_DefThr","SingleMu7",1,"13"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Muon15_DefThr","SingleMu7",1,"15"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1MuonIso3","SingleMu3",1E6,"3"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1MuonIso5","SingleMu5",1E6,"5"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1MuonIso7","SingleMu7",1E6,"7"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1MuonIso9","SingleMu7",1,"9"," "); //5E2; 
		//menu_2e31_v2->AddHlt("AptHLT1MuonIso11","SingleMu7",1,"11"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Muon3","SingleMu3",5E2,"3"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Muon5","SingleMu5",5E2,"5"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Muon7","SingleMu7",1E2,"7"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Muon9","SingleMu7",5E1,"9"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Muon11","SingleMu7",1,"11"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Muon13","SingleMu7",1,"13"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Muon15","SingleMu7",1,"15"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1L2Muon11","SingleMu7",1E1,"11"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1L2Muon16","SingleMu7",1E1,"16"," "); 
		menu_2e31_v2->AddHlt("AptHLT2Muon3","DobuleMu3",1,"(3,3)"," "); 
		menu_2e31_v2->AddHlt("HLT1MuonIso","SingleMu7",1,"11"," "); 
		menu_2e31_v2->AddHlt("HLT1MuonNonIso","SingleMu7",1,"16"," "); 
		menu_2e31_v2->AddHlt("HLT2MuonIso","DoubleMu3",1,"(3,3)"," "); 
		menu_2e31_v2->AddHlt("HLT2MuonNonIso","DoubleMu3",1,"(3,3)"," "); 
		menu_2e31_v2->AddHlt("HLT2MuonJPsi","DoubleMu3",1,"(3,3)"," "); 
		menu_2e31_v2->AddHlt("HLT2MuonUpsilon","DoubleMu3",1,"(3,3)"," "); 
		menu_2e31_v2->AddHlt("HLT2MuonZ","DoubleMu3",1,"(7,7)"," "); 
		menu_2e31_v2->AddHlt("HLTNMuonNonIso","DoubleMu3",1,"(3,3,3) "," "); 
		menu_2e31_v2->AddHlt("HLT2MuonSameSign","DoubleMu3",1,"(3,3)"," "); 
		//menu_2e31_v2->AddHlt("HLT1MuonPrescalePt3","SingleMu3",1E6,"3"," "); 
		//menu_2e31_v2->AddHlt("HLT1MuonPrescalePt5","SingleMu5",1E6,"5"," "); 
		//menu_2e31_v2->AddHlt("HLT1MuonPrescalePt7x7","SingleMu7",1E6,"7"," "); 
		//menu_2e31_v2->AddHlt("HLT1MuonPrescalePt7x10","SingleMu7",1E6,"7,10"," "); 
		menu_2e31_v2->AddHlt("HLTB1JetMu","Mu5_Jet15",1,"20 "," "); 
		menu_2e31_v2->AddHlt("HLTB2JetMu","Mu5_Jet15",1,"120"," "); 
		menu_2e31_v2->AddHlt("HLTB3JetMu","Mu5_Jet15",1,"70"," "); 
		menu_2e31_v2->AddHlt("HLTB4JetMu","Mu5_Jet15",1,"40"," "); 
		menu_2e31_v2->AddHlt("HLTBHTMu","HTT250",10,"300"," "); 
		menu_2e31_v2->AddHlt("HLTBJPsiMuMu","DoubleMu3",1,"(4,4)"," "); 
		menu_2e31_v2->AddHlt("HLTXMuonBJet","Mu5_Jet15",1,"(7,35)"," "); 
		menu_2e31_v2->AddHlt("HLTXMuonBJetSoftMuon","Mu5_Jet15",1,"(7,20)"," "); 
		menu_2e31_v2->AddHlt("HLTXMuonJets","Mu5_Jet15",1,"(7,40)"," "); 
		menu_2e31_v2->AddHlt("HLTXElectronMuon","\\star",1,"(8,7)"," "); 
		menu_2e31_v2->AddHlt("HLTXElectronMuonRelaxed","\\star",1,"(10,10)"," "); 
		menu_2e31_v2->AddHlt("HLTXMuonTau","Mu5_TauJet20",1,"(15,20)"," "); 

		//menu_2e31_v2->AddHlt("AptHLT1Level1jet15","SingleJet15",1E5,"-"," "); 
		menu_2e31_v2->AddHlt("AptHLT1jet30","SingleJet15",5E3,"30"," "); 
		menu_2e31_v2->AddHlt("AptHLT1jet40","SingleJet15",1E3,"40"," "); 
		menu_2e31_v2->AddHlt("AptHLT1jet50","SingleJet30",5E2,"50"," "); 
		menu_2e31_v2->AddHlt("AptHLT1jet60","SingleJet30",1E2,"60"," "); 
		menu_2e31_v2->AddHlt("AptHLT1jet80","SingleJet50",1E2,"80"," "); 
		menu_2e31_v2->AddHlt("AptHLT1jet110","SingleJet70",1E1,"110"," "); 
		menu_2e31_v2->AddHlt("AptHLT1jet180","SingleJet70",1,"180"," "); 
		menu_2e31_v2->AddHlt("AptHLT1jet250","SingleJet70",1,"250"," "); 
		menu_2e31_v2->AddHlt("AptHLT1jet150","SingleJet100",5,"150"," "); 
		//menu_2e31_v2->AddHlt("AptHLT2jet100","SingleJet100, DoubleJet70",1,"(100,100)"," "); 
		menu_2e31_v2->AddHlt("HLT1jet","SingleJet150",1,"200"," "); 
		menu_2e31_v2->AddHlt("HLT2jet","SingleJet150, DoubleJet70",1,"150"," "); 
		menu_2e31_v2->AddHlt("HLT3jet","\\dagger",1,"85"," "); 
		menu_2e31_v2->AddHlt("HLT4jet","\\ddagger",1,"60"," "); 
		menu_2e31_v2->AddHlt("HLT1MET","ETM40",1,"65"," "); 
		menu_2e31_v2->AddHlt("HLT2jetAco","SingleJet150, DoubleJet70",1,"125"," "); 
		menu_2e31_v2->AddHlt("HLT1jet1METAco","SingleJet150",1,"(100,60)"," "); 
		menu_2e31_v2->AddHlt("HLT1jet1MET","SingleJet150",1,"(180,60)"," "); 
		menu_2e31_v2->AddHlt("HLT2jet1MET","SingleJet150",1,"(125,60)"," "); 
		menu_2e31_v2->AddHlt("HLT3jet1MET","SingleJet150",1,"(60,60)"," "); 
		menu_2e31_v2->AddHlt("HLT4jet1MET","SingleJet150",1,"(35,60)"," "); 
		menu_2e31_v2->AddHlt("HLT1MET1HT","HTT300 ",1,"(350,65)"," "); 
		menu_2e31_v2->AddHlt("HLT1SumET","ETT60",1E1,"120"," "); 
		//menu_2e31_v2->AddHlt("HLT1jetPE7","SingleJet15",1,"30"," "); 
		//menu_2e31_v2->AddHlt("HLT1jetPE5","SingleJet30",1,"60"," "); //10000; 
		//menu_2e31_v2->AddHlt("HLT1jetPE3","SingleJet70 ",1,"110"," "); //100; 
		//menu_2e31_v2->AddHlt("HLT1jetPE1","SingleJet100 ",1,"150"," "); 
		menu_2e31_v2->AddHlt("HLT1METPre3","MinBians_HTT10",1,"15"," "); 
		menu_2e31_v2->AddHlt("HLT1METPre2","MinBians_HTT10",1,"20"," "); 
		menu_2e31_v2->AddHlt("HLT1METPre1","ETM20",1E3,"30"," "); //100; 
		//menu_2e31_v2->AddHlt("AptHLT1MET20","ETM20",1E4,"20"," "); //100; 
		//menu_2e31_v2->AddHlt("AptHLT1MET35","ETM30",1E2,"35"," "); //100; 
		//menu_2e31_v2->AddHlt("AptHLT1MET50","ETM40",1E1,"50"," "); //100; 
		//menu_2e31_v2->AddHlt("AptHLT1MET65","ETM40",1,"65"," "); //100; 
		//menu_2e31_v2->AddHlt("AptHLT1MET75","ETM40",1,"75"," "); //100; 
		//menu_2e31_v2->AddHlt("HLT2jetAve30","SingleJet15",5E2,"30"," "); 
		//menu_2e31_v2->AddHlt("HLT2jetAve60","SingleJet30",1E2,"60"," "); 
		//menu_2e31_v2->AddHlt("HLT2jetAve110","SingleJet70",1E1,"110"," "); 
		//menu_2e31_v2->AddHlt("HLT2jetAve150","SingleJet100",1,"150"," "); 
		//menu_2e31_v2->AddHlt("HLT2jetAve200","SingleJet150",1,"200"," "); 
		//menu_2e31_v2->AddHlt("AptHLT2jetAve15","SingleJet15",5E3,"15"," "); 
		//menu_2e31_v2->AddHlt("AptHLT2jetAve30","SingleJet30",5E2,"30"," "); 
		//menu_2e31_v2->AddHlt("AptHLT2jetAve50","SingleJet50",1E2,"50"," "); 
		//menu_2e31_v2->AddHlt("AptHLT2jetAve70","SingleJet70",1E1,"70"," "); 
		//menu_2e31_v2->AddHlt("AptHLT2jetAve130","SingleJet130",5,"130"," "); 
		//menu_2e31_v2->AddHlt("AptHLT2jetAve220","SingleJet220",1,"220"," "); 
		menu_2e31_v2->AddHlt("HLT2jetvbfMET","ETM30",10,"(40,60)"," "); 
		menu_2e31_v2->AddHlt("HLTS2jet1METNV","SingleJet150",1,"(-,60)"," "); 
		menu_2e31_v2->AddHlt("HLTS2jet1METAco","SingleJet150",1,"(-,70)"," "); 
		menu_2e31_v2->AddHlt("HLTSjet1MET1Aco","SingleJet150",1,"(-,70)"," "); 
		menu_2e31_v2->AddHlt("HLTSjet2MET1Aco","SingleJet150",1,"(-,70)"," "); 
		menu_2e31_v2->AddHlt("HLTS2jetAco","SingleJet150",1,"(-,-)"," "); 
		menu_2e31_v2->AddHlt("HLTJetMETRapidityGap","IsoEG10_Jet20_ForJet10",1E2,"20"," "); //100; 
		menu_2e31_v2->AddHlt("AptHLT4jet30","QuadJet15",1E2,"30"," "); 
		menu_2e31_v2->AddHlt("AptHLTXRelMu3_3jet30","Mu3_TripleJet15",1E1,"(3,30)"," "); 
		menu_2e31_v2->AddHlt("AptHLTXRelEl5_3jet30","EG5_TripleJet15",1E1,"(5,30)"," "); 

		//menu_2e31_v2->AddHlt("AptHLT1Electron5_Open","SingleEG5",1E2,"5"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Electron7_Open","SingleEG5",1E4,"7"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Electron9_Open","SingleEG5",1E4,"9"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Electron11_Open","SingleEG8",1E4,"11"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Electron13_Open","SingleEG8",1E4,"13"," "); 

		menu_2e31_v2->AddHlt("AptHLT1Electron15_Open","SingleEG10",1E1,"15"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Electron20_Open","SingleEG15",1E1,"20"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Electron25_Open","SingleEG15",1E1,"25"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Electron30_Open","SingleEG20",1,"30"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1ElectronIso8_10_DefThr","SingleIsoEG8",1E1,"10"," "); 
		menu_2e31_v2->AddHlt("AptHLT1ElectronIso8_12_DefThr","SingleIsoEG8",1,"12"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1ElectronIso8_10_DblThr","SingleIsoEG8",1E1,"10"," "); 
		menu_2e31_v2->AddHlt("AptHLT1ElectronIso8_12_DblThr","SingleIsoEG8",1,"12"," "); 
		menu_2e31_v2->AddHlt("AptHLT1ElectronIso8_15_DblThr","SingleIsoEG8",1,"15"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Electron12_15_DefThr","SingleEG12",1,"15"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Electron12_15_DblThr","SingleEG12",1,"15"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Electron12_17_DblThr","SingleEG12",1,"17"," "); 
		menu_2e31_v2->AddHlt("AptHLT2ElectronIso8_10_Open","DoubleIsoEG8",1,"10"," "); 
		menu_2e31_v2->AddHlt("AptHLT2Electron10_12_Open","DoubleEG10",1,"12"," "); 
		menu_2e31_v2->AddHlt("HLT1Electron","SingleIsoEG12",1,"15"," "); 
		menu_2e31_v2->AddHlt("HLT1ElectronRelaxed","SingleEG15",1,"18"," "); 
		menu_2e31_v2->AddHlt("HLT2Electron","DoubleIsoEG8",1,"10"," "); 
		menu_2e31_v2->AddHlt("HLT2ElectronRelaxed","DoubleEG10",1,"12"," "); 

		//menu_2e31_v2->AddHlt("AptHLT1Photon5NoTrkIso","SingleEG5",5E3,"5"," ");
		menu_2e31_v2->AddHlt("AptHLT1Photon10NoTrkIso","SingleEG8",1E3,"10"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon15NoTrkIso","SingleEG10",1E2,"15"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Photon20NoTrkIso","SingleEG10",1E2,"20"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon25NoTrkIso","SingleEG10",1E2,"25"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Photon30NoTrkIso","SingleEG20",1E2,"30"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon35NoTrkIso","SingleEG20",1E2,"35"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Photon40NoTrkIso","SingleEG25",1,"40"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon5_Open","SingleEG5",10,"5"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon10_Open","SingleEG8",1,"10"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon12NoTrkIso","SingleEG8",1,"12"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon12_Open","SingleEG8",1,"12"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon15"_Open,"SingleEG10",1,"15"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon25"_Open,"SingleEG10",1,"25"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon30"_Open,"SingleEG20",1,"30"," "); 
		//menu_2e31_v2->AddHlt("AptHLT2Photon15"_Open,"DoubleEG5",1,"15"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1PhotonIso8_15_DefThr","SingleIsoEG8",1,"15"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1PhotonIso8_17_DefThr","SingleIsoEG8",1,"17"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1PhotonIso8_20_DefThr","SingleIsoEG8",1,"20"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1PhotonIso8_25_DefThr","SingleIsoEG8",1,"25"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1PhotonIso8_30_DefThr","SingleIsoEG8",1,"30"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon12_20_DefThr","SingleEG12",1,"20"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon12_30_DefThr","SingleEG12",1,"30"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon12_35_DefThr","SingleEG12",1,"35"," "); 
		//menu_2e31_v2->AddHlt("AptHLT1Photon12_40_DefThr","SingleEG12",1,"40"," "); 
		menu_2e31_v2->AddHlt("AptHLT1PhotonIso8_20_DblThr","SingleIsoEG8",1,"20"," "); 
		menu_2e31_v2->AddHlt("AptHLT1PhotonIso8_25_DblThr","SingleIsoEG8",1,"25"," "); 
		menu_2e31_v2->AddHlt("AptHLT1PhotonIso8_30_DblThr","SingleIsoEG8",1,"30"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Photon12_15_DblThr","SingleEG12",1,"15"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Photon12_25_DefThr","SingleEG12",1,"25"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Photon12_25_DblThr","SingleEG12",1,"25"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Photon12_30_DblThr","SingleEG12",1,"30"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Photon12_35_DblThr","SingleEG12",1,"35"," "); 
		menu_2e31_v2->AddHlt("AptHLT1Photon12_40_DblThr","SingleEG12",1,"40"," "); 
		menu_2e31_v2->AddHlt("AptHLT2PhotonIso8_15_Open","DoubleIsoEG8",1,"15"," "); 
		menu_2e31_v2->AddHlt("AptHLT2Photon10_20_Open","DoubleEG10",1,"20"," "); 
		menu_2e31_v2->AddHlt("AptHLT2Photon10_20NoTrkIso","DoubleEG10",1,"(20,20)"," "); 
		menu_2e31_v2->AddHlt("AptHLT2Photon10_30NoTrkIso","DoubleEG10",1,"(30,30)"," "); 
		menu_2e31_v2->AddHlt("HLT1Photon","SingleIsoEG12",1,"30"," "); 
		menu_2e31_v2->AddHlt("HLT1PhotonRelaxed","SingleEG15",1,"40"," "); 
		menu_2e31_v2->AddHlt("HLT2Photon","DoubleIsoEG8",1,"(20,20)"," "); 
		menu_2e31_v2->AddHlt("HLT2PhotonRelaxed","DoubleEG10",1,"(20,20)"," "); 
		menu_2e31_v2->AddHlt("HLT1EMHighEt","SingleEG15",1,"80"," "); 
		menu_2e31_v2->AddHlt("HLT1EMVeryHighEt","SingleEG15",1,"200"," "); 
		menu_2e31_v2->AddHlt("HLT2ElectronZCounter","DoubleIsoEG8",1,"(10,10)"," "); 
		menu_2e31_v2->AddHlt("HLT2ElectronExclusive","ExclusiveDoubleIsoEG6",1,"(6,6)"," "); 
		menu_2e31_v2->AddHlt("HLT2PhotonExclusive","ExclusiveDoubleIsoEG6",1,"(10,10)"," "); 
		menu_2e31_v2->AddHlt("HLT1PhotonL1Isolated","SingleIsoEG10",1E2,"12"," "); 

		menu_2e31_v2->AddHlt("HLTB1Jet","\\Diamond ",1,"180"," "); 
		menu_2e31_v2->AddHlt("HLTB2Jet","\\Diamond ",1,"120"," "); 
		menu_2e31_v2->AddHlt("HLTB3Jet","\\Diamond ",1,"70"," "); 
		menu_2e31_v2->AddHlt("HLTB4Jet","\\Diamond ",1,"40"," "); 
		menu_2e31_v2->AddHlt("HLTBHT","\\Diamond ",1,"470"," "); 
		menu_2e31_v2->AddHlt("HLT1Tau","SingleTauJet80",1,"15"," "); 
		menu_2e31_v2->AddHlt("HLT1Tau1MET","TauJet30_ETM30",1,"15"," "); 
		menu_2e31_v2->AddHlt("HLT2TauPixel","TauJet40",1,"15"," "); 
		menu_2e31_v2->AddHlt("HLTXElectronBJet","IsoEG10_Jet20",1,"(10,35)"," "); 
		menu_2e31_v2->AddHlt("HLTXElectron1Jet","IsoEG10_Jet30",1,"(12,40)"," "); 
		menu_2e31_v2->AddHlt("HLTXElectron2Jet","IsoEG10_Jet30",1,"(12,80)"," "); 
		menu_2e31_v2->AddHlt("HLTXElectron3Jet","IsoEG10_Jet30",1,"(12,60)"," "); 
		menu_2e31_v2->AddHlt("HLTXElectron4Jet","IsoEG10_Jet30",1,"(12,35)"," "); 
		menu_2e31_v2->AddHlt("HLTXElectronTau","IsoEG10_TauJet20",1,"(12,20)"," "); 

		menu_2e31_v2->AddHlt("HLTHcalIsolatedTrack","NA",1,"NA"," "); 
		//menu_2e31_v2->AddHlt("HLTHcalIsolatedTrackNoEcalIsol"," ",1," "," ");  
		menu_2e31_v2->AddHlt("HLTMinBiasPixel","ZeroBias",1,"- "," "); 
		//menu_2e31_v2->AddHlt("HLTMinBiasForAlignment"," ",1," "," "); 
		menu_2e31_v2->AddHlt("HLTMinBias","MinBias_HTT10",1,"-"," "); 
		menu_2e31_v2->AddHlt("HLTZeroBias","ZeroBias",1,"-"," "); 

}

void BookOHltMenu_2e30_v03(OHltMenu*  menu_2e30_v03, double &iLumi, double &nBunches) {

		iLumi = 2E30;
		nBunches = 43;

		menu_2e30_v03->AddHlt("AptHLT1MuonLevel1Open","SingleMuOpen",1.5E2,"-"," "); //1000; 
		menu_2e30_v03->AddHlt("AptHLT1MuonLevel1","SingleMu7,DoubleMu3",2E1,"-"," "); //1000; 
		//menu_2e30_v03->AddHlt("AptHLT1MuonIso3_DefThr","SingleMu3",1,"3"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1MuonIso5_DefThr","SingleMu5",1,"5"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1MuonIso7_DefThr","SingleMu7",1,"7"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1MuonIso9_DefThr","SingleMu7",1,"9"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1MuonIso11_DefThr","SingleMu7",1,"11"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Muon3_DefThr","SingleMu3",1,"3"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Muon5_DefThr","SingleMu5",1,"5"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Muon7_DefThr","SingleMu7",1,"7"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Muon9_DefThr","SingleMu7",1,"9"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Muon11_DefThr","SingleMu7",1,"11"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Muon13_DefThr","SingleMu7",1,"13"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Muon15_DefThr","SingleMu7",1,"15"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1MuonIso3","SingleMu3",1E6,"3"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1MuonIso5","SingleMu5",1E6,"5"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1MuonIso7","SingleMu7",1E6,"7"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1MuonIso9","SingleMu7",1,"9"," "); //5E2; 
		//menu_2e30_v03->AddHlt("AptHLT1MuonIso11","SingleMu7",1,"11"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Muon3","SingleMu3",10,"3"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Muon5","SingleMu5",1,"5"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Muon7","SingleMu7",1,"7"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Muon9","SingleMu7",1,"9"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Muon11","SingleMu7",1,"11"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Muon13","SingleMu7",1,"13"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Muon15","SingleMu7",1,"15"," "); 
		menu_2e30_v03->AddHlt("AptHLT1L2Muon11","SingleMu7",1,"11"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1L2Muon16","SingleMu7",1,"16"," "); 
		menu_2e30_v03->AddHlt("AptHLT2Muon3","DobuleMu3",1,"(3,3)"," "); 
		//menu_2e30_v03->AddHlt("HLT1MuonLevel1","SingleMu0",1,"-"," "); 
		menu_2e30_v03->AddHlt("HLT1MuonIso","SingleMu7",1,"11"," "); 
		menu_2e30_v03->AddHlt("HLT1MuonNonIso","SingleMu7",1,"16"," "); 
		menu_2e30_v03->AddHlt("HLT2MuonIso","DoubleMu3",1,"(3,3)"," "); 
		menu_2e30_v03->AddHlt("HLT2MuonNonIso","DoubleMu3",1,"(3,3)"," "); 
		menu_2e30_v03->AddHlt("HLT2MuonJPsi","DoubleMu3",1,"(3,3)"," "); 
		menu_2e30_v03->AddHlt("HLT2MuonUpsilon","DoubleMu3",1,"(3,3)"," "); 
		menu_2e30_v03->AddHlt("HLT2MuonZ","DoubleMu3",1,"(7,7)"," "); 
		menu_2e30_v03->AddHlt("HLTNMuonNonIso","DoubleMu3",1,"(3,3,3) "," "); 
		menu_2e30_v03->AddHlt("HLT2MuonSameSign","DoubleMu3",1,"(3,3)"," "); 
		//menu_2e30_v03->AddHlt("HLT1MuonPrescalePt3","SingleMu3",1E6,"3"," "); 
		//menu_2e30_v03->AddHlt("HLT1MuonPrescalePt5","SingleMu5",1E6,"5"," "); 
		//menu_2e30_v03->AddHlt("HLT1MuonPrescalePt7x7","SingleMu7",1E6,"7"," "); 
		//menu_2e30_v03->AddHlt("HLT1MuonPrescalePt7x10","SingleMu7",1E6,"7,10"," "); 
		menu_2e30_v03->AddHlt("HLTB1JetMu","Mu5_Jet15",1,"20 "," "); 
		menu_2e30_v03->AddHlt("HLTB2JetMu","Mu5_Jet15",1,"120"," "); 
		menu_2e30_v03->AddHlt("HLTB3JetMu","Mu5_Jet15",1,"70"," "); 
		menu_2e30_v03->AddHlt("HLTB4JetMu","Mu5_Jet15",1,"40"," "); 
		menu_2e30_v03->AddHlt("HLTBHTMu","HTT250",1,"300"," "); 
		menu_2e30_v03->AddHlt("HLTBJPsiMuMu","DoubleMu3",1,"(4,4)"," "); 
		menu_2e30_v03->AddHlt("HLTXMuonBJet","Mu5_Jet15",1,"(7,35)"," "); 
		menu_2e30_v03->AddHlt("HLTXMuonBJetSoftMuon","Mu5_Jet15",1,"(7,20)"," "); 
		menu_2e30_v03->AddHlt("HLTXMuonJets","Mu5_Jet15",1,"(7,40)"," "); 
		menu_2e30_v03->AddHlt("HLTXElectronMuon","\\star",1,"(8,7)"," "); 
		menu_2e30_v03->AddHlt("HLTXElectronMuonRelaxed","\\star",1,"(10,10)"," "); 
		menu_2e30_v03->AddHlt("HLTXMuonTau","Mu5_TauJet20",1,"(15,20)"," "); 

		menu_2e30_v03->AddHlt("AptHLT1Level1jet15","SingleJet15",1E4,"-"," "); 
		menu_2e30_v03->AddHlt("AptHLT1jet30","SingleJet15",2E2,"30"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1jet40","SingleJet15",1E2,"40"," "); 
		menu_2e30_v03->AddHlt("AptHLT1jet50","SingleJet30",2E1,"50"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1jet60","SingleJet30",1E1,"60"," "); 
		menu_2e30_v03->AddHlt("AptHLT1jet80","SingleJet50",1E1,"80"," "); 
		menu_2e30_v03->AddHlt("AptHLT1jet110","SingleJet70",1,"110"," "); 
		menu_2e30_v03->AddHlt("AptHLT1jet180","SingleJet70",1,"180"," "); 
		menu_2e30_v03->AddHlt("AptHLT1jet250","SingleJet70",1,"250"," "); 
		menu_2e30_v03->AddHlt("AptHLT2jet100","SingleJet100, DoubleJet70",1,"(100,100)"," "); 
		menu_2e30_v03->AddHlt("HLT1jet","SingleJet150",1,"200"," "); 
		menu_2e30_v03->AddHlt("HLT2jet","SingleJet150, DoubleJet70",1,"150"," "); 
		menu_2e30_v03->AddHlt("HLT3jet","\\dagger",1,"85"," "); 
		menu_2e30_v03->AddHlt("HLT4jet","\\ddagger",1,"60"," "); 
		menu_2e30_v03->AddHlt("HLT1MET","ETM40",1,"65"," "); 
		menu_2e30_v03->AddHlt("HLT2jetAco","SingleJet150, DoubleJet70",1,"125"," "); 
		menu_2e30_v03->AddHlt("HLT1jet1METAco","SingleJet150",1,"(100,60)"," "); 
		menu_2e30_v03->AddHlt("HLT1jet1MET","SingleJet150",1,"(180,60)"," "); 
		menu_2e30_v03->AddHlt("HLT2jet1MET","SingleJet150",1,"(125,60)"," "); 
		menu_2e30_v03->AddHlt("HLT3jet1MET","SingleJet150",1,"(60,60)"," "); 
		menu_2e30_v03->AddHlt("HLT4jet1MET","SingleJet150",1,"(35,60)"," "); 
		menu_2e30_v03->AddHlt("HLT1MET1HT","HTT300 ",1,"(350,65)"," "); 
		menu_2e30_v03->AddHlt("HLT1SumET","ETT60",1,"120"," "); 
		//menu_2e30_v03->AddHlt("HLT1jetPE7","SingleJet15",1,"30"," "); 
		//menu_2e30_v03->AddHlt("HLT1jetPE5","SingleJet30",1,"60"," "); //10000; 
		//menu_2e30_v03->AddHlt("HLT1jetPE3","SingleJet70 ",1,"110"," "); //100; 
		//menu_2e30_v03->AddHlt("HLT1jetPE1","SingleJet100 ",1,"150"," "); 
		//menu_2e30_v03->AddHlt("HLT1METPre3","MinBians_HTT10",1,"15"," "); 
		//menu_2e30_v03->AddHlt("HLT1METPre2","MinBians_HTT10",1,"20"," "); 
		//menu_2e30_v03->AddHlt("HLT1METPre1","ETM20",1E2,"30"," "); //100; 
		menu_2e30_v03->AddHlt("AptHLT1MET20","ETM20",1E3,"20"," "); //100; 
		menu_2e30_v03->AddHlt("AptHLT1MET35","ETM30",1E1,"35"," "); //100; 
		menu_2e30_v03->AddHlt("AptHLT1MET50","ETM40",1,"50"," "); //100; 
		menu_2e30_v03->AddHlt("AptHLT1MET65","ETM40",1,"65"," "); //100; 
		menu_2e30_v03->AddHlt("AptHLT1MET75","ETM40",1,"75"," "); //100; 
		//menu_2e30_v03->AddHlt("HLT2jetAve30","SingleJet15",1,"30"," "); 
		//menu_2e30_v03->AddHlt("HLT2jetAve60","SingleJet30",1,"60"," "); 
		//menu_2e30_v03->AddHlt("HLT2jetAve110","SingleJet70",1,"110"," "); 
		//menu_2e30_v03->AddHlt("HLT2jetAve150","SingleJet100",1,"150"," "); 
		//menu_2e30_v03->AddHlt("HLT2jetAve200","SingleJet150",1,"200"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2jetAve15","SingleJet15",4E2,"15"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2jetAve30","SingleJet30",2E1,"30"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2jetAve50","SingleJet50",4E1,"50"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2jetAve70","SingleJet70",1,"70"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2jetAve130","SingleJet130",1,"130"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2jetAve220","SingleJet220",1,"220"," "); 
		menu_2e30_v03->AddHlt("HLT2jetvbfMET","ETM30",1,"(40,60)"," "); 
		menu_2e30_v03->AddHlt("HLTS2jet1METNV","SingleJet150",1,"(-,60)"," "); 
		menu_2e30_v03->AddHlt("HLTS2jet1METAco","SingleJet150",1,"(-,70)"," "); 
		menu_2e30_v03->AddHlt("HLTSjet1MET1Aco","SingleJet150",1,"(-,70)"," "); 
		menu_2e30_v03->AddHlt("HLTSjet2MET1Aco","SingleJet150",1,"(-,70)"," "); 
		menu_2e30_v03->AddHlt("HLTS2jetAco","SingleJet150",1,"(-,-)"," "); 
		menu_2e30_v03->AddHlt("HLTJetMETRapidityGap","IsoEG10_Jet20_ForJet10",1,"20"," "); //100; 
		menu_2e30_v03->AddHlt("AptHLT4jet30","QuadJet15",1E1,"30"," "); 
		menu_2e30_v03->AddHlt("AptHLTXRelMu3_3jet30","Mu3_TripleJet20",1,"(3,30)"," "); 
		menu_2e30_v03->AddHlt("AptHLTXRelEl5_3jet30","EG5_TripleJet20",1,"(5,30)"," "); 

		menu_2e30_v03->AddHlt("AptHLT1Electron5_Open","SingleEG5",1E1,"5"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Electron8_DefThr","SingleEG5",1,"8"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Electron8_Open","SingleEG5",1,"8"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Electron10_Open","SingleEG8",1,"10"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Electron7_Open","SingleEG5",1E4,"7"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Electron9_Open","SingleEG5",1E4,"9"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Electron11_Open","SingleEG8",1E4,"11"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Electron13_Open","SingleEG8",1E4,"13"," "); 

		menu_2e30_v03->AddHlt("AptHLT1Electron15_Open","SingleEG10",1,"15"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Electron20_Open","SingleEG15",1,"20"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Electron25_Open","SingleEG15",1,"25"," "); 
		menu_2e30_v03->AddHlt("AptHLT2Electron5_Open","DoubleEG5",1,"5"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Electron30_Open","SingleEG20",1,"30"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1ElectronIso8_10_DefThr","SingleIsoEG8",1E1,"10"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1ElectronIso8_12_DefThr","SingleIsoEG8",1,"12"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1ElectronIso8_10_DblThr","SingleIsoEG8",1E1,"10"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1ElectronIso8_12_DblThr","SingleIsoEG8",1,"12"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1ElectronIso8_15_DblThr","SingleIsoEG8",1,"15"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Electron12_15_DefThr","SingleEG12",1,"15"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Electron12_15_DblThr","SingleEG12",1,"15"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Electron12_17_DblThr","SingleEG12",1,"17"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2ElectronIso8_10_Open","DoubleIsoEG8",1,"10"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2Electron10_12_Open","DoubleEG10",1,"12"," "); 
		menu_2e30_v03->AddHlt("HLT1Electron","SingleIsoEG12",1,"15"," "); 
		menu_2e30_v03->AddHlt("HLT1ElectronRelaxed","SingleEG15",1,"18"," "); 
		menu_2e30_v03->AddHlt("HLT2Electron","DoubleIsoEG8",1,"10"," "); 
		menu_2e30_v03->AddHlt("HLT2ElectronRelaxed","DoubleEG10",1,"12"," "); 

		//menu_2e30_v03->AddHlt("AptHLT1Photon5NoTrkIso","SingleEG5",5E3,"5"," ");
		menu_2e30_v03->AddHlt("AptHLT1Photon10_DefThr","SingleEG8",1E1,"10"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Photon10NoTrkIso","SingleEG8",1E2,"10"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon15NoTrkIso","SingleEG10",1E2,"15"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Photon20NoTrkIso","SingleEG10",1,"20"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Photon25NoTrkIso","SingleEG15",1,"25"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Photon30NoTrkIso","SingleEG20",1,"30"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon35NoTrkIso","SingleEG20",1E2,"35"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon40NoTrkIso","SingleEG25",2,"40"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon5_Open","SingleEG5",10,"5"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon10_Open","SingleEG8",1,"10"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon12NoTrkIso","SingleEG8",1,"12"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon12_Open","SingleEG8",1,"12"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Photon15_Open","SingleEG10",1,"15"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Photon25_Open","SingleEG10",1,"25"," "); 
		menu_2e30_v03->AddHlt("AptHLT1Photon30_Open","SingleEG20",1,"30"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2Photon15_Open","DoubleEG5",1,"15"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1PhotonIso8_15_DefThr","SingleIsoEG8",1,"15"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1PhotonIso8_17_DefThr","SingleIsoEG8",1,"17"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1PhotonIso8_20_DefThr","SingleIsoEG8",1,"20"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1PhotonIso8_25_DefThr","SingleIsoEG8",1,"25"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1PhotonIso8_30_DefThr","SingleIsoEG8",1,"30"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon12_20_DefThr","SingleEG12",1,"20"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon12_30_DefThr","SingleEG12",1,"30"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon12_35_DefThr","SingleEG12",1,"35"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon12_40_DefThr","SingleEG12",1,"40"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1PhotonIso8_20_DblThr","SingleIsoEG8",1,"20"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1PhotonIso8_25_DblThr","SingleIsoEG8",1,"25"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1PhotonIso8_30_DblThr","SingleIsoEG8",1,"30"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon12_15_DblThr","SingleEG12",5,"15"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon12_25_DefThr","SingleEG12",1,"25"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon12_25_DblThr","SingleEG12",1,"25"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon12_30_DblThr","SingleEG12",1,"30"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon12_35_DblThr","SingleEG12",1,"35"," "); 
		//menu_2e30_v03->AddHlt("AptHLT1Photon12_40_DblThr","SingleEG12",1,"40"," "); 
		menu_2e30_v03->AddHlt("AptHLT2Photon10NoTrkIso","DoubleEG5",1,"(10,10)"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2PhotonIso8_15_Open","DoubleIsoEG8",1,"15"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2Photon10_20_Open","DoubleEG10",1,"20"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2Photon10_20NoTrkIso","DoubleEG10",1,"(20,20)"," "); 
		//menu_2e30_v03->AddHlt("AptHLT2Photon10_30NoTrkIso","DoubleEG10",1,"(30,30)"," "); 
		menu_2e30_v03->AddHlt("HLT1Photon","SingleIsoEG12",1,"30"," "); 
		menu_2e30_v03->AddHlt("HLT1PhotonRelaxed","SingleEG15",1,"40"," "); 
		menu_2e30_v03->AddHlt("HLT2Photon","DoubleIsoEG8",1,"(20,20)"," "); 
		menu_2e30_v03->AddHlt("HLT2PhotonRelaxed","DoubleEG10",1,"(20,20)"," "); 
		menu_2e30_v03->AddHlt("HLT1EMHighEt","SingleEG15",1,"80"," "); 
		menu_2e30_v03->AddHlt("HLT1EMVeryHighEt","SingleEG15",1,"200"," "); 
		menu_2e30_v03->AddHlt("HLT2ElectronZCounter","DoubleIsoEG8",1,"(10,10)"," "); 
		menu_2e30_v03->AddHlt("HLT2ElectronExclusive","ExclusiveDoubleIsoEG6",1,"(6,6)"," "); 
		menu_2e30_v03->AddHlt("HLT2PhotonExclusive","ExclusiveDoubleIsoEG6",1,"(10,10)"," "); 
		menu_2e30_v03->AddHlt("HLT1PhotonL1Isolated","SingleIsoEG10",1E1,"12"," "); 

		menu_2e30_v03->AddHlt("HLTB1Jet","\\Diamond ",1,"180"," "); 
		menu_2e30_v03->AddHlt("HLTB2Jet","\\Diamond ",1,"120"," "); 
		menu_2e30_v03->AddHlt("HLTB3Jet","\\Diamond ",1,"70"," "); 
		menu_2e30_v03->AddHlt("HLTB4Jet","\\Diamond ",1,"40"," "); 
		menu_2e30_v03->AddHlt("HLTBHT","\\Diamond ",1,"470"," "); 
		menu_2e30_v03->AddHlt("HLT1Tau","SingleTauJet80",1,"15"," "); 
		menu_2e30_v03->AddHlt("HLT1Tau1MET","TauJet30_ETM30",1,"15"," "); 
		menu_2e30_v03->AddHlt("HLT2TauPixel","TauJet40",1,"15"," "); 
		menu_2e30_v03->AddHlt("HLTXElectronBJet","IsoEG10_Jet20",1,"(10,35)"," "); 
		menu_2e30_v03->AddHlt("HLTXElectron1Jet","IsoEG10_Jet30",1,"(12,40)"," "); 
		menu_2e30_v03->AddHlt("HLTXElectron2Jet","IsoEG10_Jet30",1,"(12,80)"," "); 
		menu_2e30_v03->AddHlt("HLTXElectron3Jet","IsoEG10_Jet30",1,"(12,60)"," "); 
		menu_2e30_v03->AddHlt("HLTXElectron4Jet","IsoEG10_Jet30",1,"(12,35)"," "); 
		menu_2e30_v03->AddHlt("HLTXElectronTau","IsoEG10_TauJet20",1,"(12,20)"," "); 

		//menu_2e30_v03->AddHlt("HLTHcalIsolatedTrack","NA",1,"NA"," "); 
		//menu_2e30_v03->AddHlt("HLTHcalIsolatedTrackNoEcalIsol"," ",1," "," ");  
		menu_2e30_v03->AddHlt("HLTMinBiasPixel","ZeroBias",1,"- "," "); 
		//menu_2e30_v03->AddHlt("HLTMinBiasForAlignment"," ",1," "," "); 
		menu_2e30_v03->AddHlt("HLTMinBias","MinBias_HTT10",1,"-"," "); 
		menu_2e30_v03->AddHlt("HLTZeroBias","ZeroBias",1,"-"," "); 

}

void BookOHltMenu_2e30_v02(OHltMenu*  menu_2e30_v02, double &iLumi, double &nBunches) {

		iLumi = 2E30;
		nBunches = 43;

		menu_2e30_v02->AddHlt("AptHLT1MuonLevel1Open","SingleMuOpen",1.5E2,"-"," "); //1000; 
		menu_2e30_v02->AddHlt("AptHLT1MuonLevel1","SingleMu7,DoubleMu3",2E1,"-"," "); //1000; 
		//menu_2e30_v02->AddHlt("AptHLT1MuonIso3_DefThr","SingleMu3",1,"3"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1MuonIso5_DefThr","SingleMu5",1,"5"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1MuonIso7_DefThr","SingleMu7",1,"7"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1MuonIso9_DefThr","SingleMu7",1,"9"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1MuonIso11_DefThr","SingleMu7",1,"11"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Muon3_DefThr","SingleMu3",1,"3"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Muon5_DefThr","SingleMu5",1,"5"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Muon7_DefThr","SingleMu7",1,"7"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Muon9_DefThr","SingleMu7",1,"9"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Muon11_DefThr","SingleMu7",1,"11"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Muon13_DefThr","SingleMu7",1,"13"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Muon15_DefThr","SingleMu7",1,"15"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1MuonIso3","SingleMu3",1E6,"3"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1MuonIso5","SingleMu5",1E6,"5"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1MuonIso7","SingleMu7",1E6,"7"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1MuonIso9","SingleMu7",1,"9"," "); //5E2; 
		//menu_2e30_v02->AddHlt("AptHLT1MuonIso11","SingleMu7",1,"11"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Muon3","SingleMu3",10,"3"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Muon5","SingleMu5",1,"5"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Muon7","SingleMu7",1,"7"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Muon9","SingleMu7",1,"9"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Muon11","SingleMu7",1,"11"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Muon13","SingleMu7",1,"13"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Muon15","SingleMu7",1,"15"," "); 
		menu_2e30_v02->AddHlt("AptHLT1L2Muon11","SingleMu7",1,"11"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1L2Muon16","SingleMu7",1,"16"," "); 
		menu_2e30_v02->AddHlt("AptHLT2Muon3","DobuleMu3",1,"(3,3)"," "); 
		//menu_2e30_v02->AddHlt("HLT1MuonLevel1","SingleMu0",1,"-"," "); 
		menu_2e30_v02->AddHlt("HLT1MuonIso","SingleMu7",1,"11"," "); 
		menu_2e30_v02->AddHlt("HLT1MuonNonIso","SingleMu7",1,"16"," "); 
		menu_2e30_v02->AddHlt("HLT2MuonIso","DoubleMu3",1,"(3,3)"," "); 
		menu_2e30_v02->AddHlt("HLT2MuonNonIso","DoubleMu3",1,"(3,3)"," "); 
		menu_2e30_v02->AddHlt("HLT2MuonJPsi","DoubleMu3",1,"(3,3)"," "); 
		menu_2e30_v02->AddHlt("HLT2MuonUpsilon","DoubleMu3",1,"(3,3)"," "); 
		menu_2e30_v02->AddHlt("HLT2MuonZ","DoubleMu3",1,"(7,7)"," "); 
		menu_2e30_v02->AddHlt("HLTNMuonNonIso","DoubleMu3",1,"(3,3,3) "," "); 
		menu_2e30_v02->AddHlt("HLT2MuonSameSign","DoubleMu3",1,"(3,3)"," "); 
		//menu_2e30_v02->AddHlt("HLT1MuonPrescalePt3","SingleMu3",1E6,"3"," "); 
		//menu_2e30_v02->AddHlt("HLT1MuonPrescalePt5","SingleMu5",1E6,"5"," "); 
		//menu_2e30_v02->AddHlt("HLT1MuonPrescalePt7x7","SingleMu7",1E6,"7"," "); 
		//menu_2e30_v02->AddHlt("HLT1MuonPrescalePt7x10","SingleMu7",1E6,"7,10"," "); 
		menu_2e30_v02->AddHlt("HLTB1JetMu","Mu5_Jet15",1,"20 "," "); 
		menu_2e30_v02->AddHlt("HLTB2JetMu","Mu5_Jet15",1,"120"," "); 
		menu_2e30_v02->AddHlt("HLTB3JetMu","Mu5_Jet15",1,"70"," "); 
		menu_2e30_v02->AddHlt("HLTB4JetMu","Mu5_Jet15",1,"40"," "); 
		menu_2e30_v02->AddHlt("HLTBHTMu","HTT250",1,"300"," "); 
		menu_2e30_v02->AddHlt("HLTBJPsiMuMu","DoubleMu3",1,"(4,4)"," "); 
		menu_2e30_v02->AddHlt("HLTXMuonBJet","Mu5_Jet15",1,"(7,35)"," "); 
		menu_2e30_v02->AddHlt("HLTXMuonBJetSoftMuon","Mu5_Jet15",1,"(7,20)"," "); 
		menu_2e30_v02->AddHlt("HLTXMuonJets","Mu5_Jet15",1,"(7,40)"," "); 
		menu_2e30_v02->AddHlt("HLTXElectronMuon","\\star",1,"(8,7)"," "); 
		menu_2e30_v02->AddHlt("HLTXElectronMuonRelaxed","\\star",1,"(10,10)"," "); 
		menu_2e30_v02->AddHlt("HLTXMuonTau","Mu5_TauJet20",1,"(15,20)"," "); 

		//menu_2e30_v02->AddHlt("AptHLT1Level1jet15","SingleJet15",1E4,"-"," "); 
		menu_2e30_v02->AddHlt("AptHLT1jet30","SingleJet15",2E2,"30"," "); 
		menu_2e30_v02->AddHlt("AptHLT1jet40","SingleJet15",1E2,"40"," "); 
		menu_2e30_v02->AddHlt("AptHLT1jet50","SingleJet30",2E1,"50"," "); 
		menu_2e30_v02->AddHlt("AptHLT1jet60","SingleJet30",1E1,"60"," "); 
		menu_2e30_v02->AddHlt("AptHLT1jet80","SingleJet50",1,"80"," "); 
		menu_2e30_v02->AddHlt("AptHLT1jet110","SingleJet70",1,"110"," "); 
		menu_2e30_v02->AddHlt("AptHLT1jet150","SingleJet70",1,"150"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1jet180","SingleJet70",1,"180"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1jet250","SingleJet70",1,"250"," "); 
		menu_2e30_v02->AddHlt("AptHLT2jet100","SingleJet100, DoubleJet70",1,"(100,100)"," "); 
		menu_2e30_v02->AddHlt("HLT1jet","SingleJet150",1,"200"," "); 
		menu_2e30_v02->AddHlt("HLT2jet","SingleJet150, DoubleJet70",1,"150"," "); 
		menu_2e30_v02->AddHlt("HLT3jet","\\dagger",1,"85"," "); 
		menu_2e30_v02->AddHlt("HLT4jet","\\ddagger",1,"60"," "); 
		menu_2e30_v02->AddHlt("HLT1MET","ETM40",1,"65"," "); 
		menu_2e30_v02->AddHlt("HLT2jetAco","SingleJet150, DoubleJet70",1,"125"," "); 
		menu_2e30_v02->AddHlt("HLT1jet1METAco","SingleJet150",1,"(100,60)"," "); 
		menu_2e30_v02->AddHlt("HLT1jet1MET","SingleJet150",1,"(180,60)"," "); 
		menu_2e30_v02->AddHlt("HLT2jet1MET","SingleJet150",1,"(125,60)"," "); 
		menu_2e30_v02->AddHlt("HLT3jet1MET","SingleJet150",1,"(60,60)"," "); 
		menu_2e30_v02->AddHlt("HLT4jet1MET","SingleJet150",1,"(35,60)"," "); 
		menu_2e30_v02->AddHlt("HLT1MET1HT","HTT300 ",1,"(350,65)"," "); 
		menu_2e30_v02->AddHlt("HLT1SumET","ETT60",1,"120"," "); 
		//menu_2e30_v02->AddHlt("HLT1jetPE7","SingleJet15",1,"30"," "); 
		//menu_2e30_v02->AddHlt("HLT1jetPE5","SingleJet30",1,"60"," "); //10000; 
		//menu_2e30_v02->AddHlt("HLT1jetPE3","SingleJet70 ",1,"110"," "); //100; 
		//menu_2e30_v02->AddHlt("HLT1jetPE1","SingleJet100 ",1,"150"," "); 
		menu_2e30_v02->AddHlt("HLT1METPre3","MinBians_HTT10",1,"15"," "); 
		menu_2e30_v02->AddHlt("HLT1METPre2","MinBians_HTT10",1,"20"," "); 
		menu_2e30_v02->AddHlt("HLT1METPre1","ETM20",1E2,"30"," "); //100; 
		//menu_2e30_v02->AddHlt("AptHLT1MET20","ETM20",1E3,"20"," "); //100; 
		//menu_2e30_v02->AddHlt("AptHLT1MET35","ETM30",1E1,"35"," "); //100; 
		//menu_2e30_v02->AddHlt("AptHLT1MET50","ETM40",1,"50"," "); //100; 
		//menu_2e30_v02->AddHlt("AptHLT1MET65","ETM40",1,"65"," "); //100; 
		//menu_2e30_v02->AddHlt("AptHLT1MET75","ETM40",1,"75"," "); //100; 
		//menu_2e30_v02->AddHlt("HLT2jetAve30","SingleJet15",1,"30"," "); 
		//menu_2e30_v02->AddHlt("HLT2jetAve60","SingleJet30",1,"60"," "); 
		//menu_2e30_v02->AddHlt("HLT2jetAve110","SingleJet70",1,"110"," "); 
		//menu_2e30_v02->AddHlt("HLT2jetAve150","SingleJet100",1,"150"," "); 
		//menu_2e30_v02->AddHlt("HLT2jetAve200","SingleJet150",1,"200"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2jetAve15","SingleJet15",4E2,"15"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2jetAve30","SingleJet30",2E1,"30"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2jetAve50","SingleJet50",4E1,"50"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2jetAve70","SingleJet70",1,"70"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2jetAve130","SingleJet130",1,"130"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2jetAve220","SingleJet220",1,"220"," "); 
		menu_2e30_v02->AddHlt("HLT2jetvbfMET","ETM30",1,"(40,60)"," "); 
		menu_2e30_v02->AddHlt("HLTS2jet1METNV","SingleJet150",1,"(-,60)"," "); 
		menu_2e30_v02->AddHlt("HLTS2jet1METAco","SingleJet150",1,"(-,70)"," "); 
		menu_2e30_v02->AddHlt("HLTSjet1MET1Aco","SingleJet150",1,"(-,70)"," "); 
		menu_2e30_v02->AddHlt("HLTSjet2MET1Aco","SingleJet150",1,"(-,70)"," "); 
		menu_2e30_v02->AddHlt("HLTS2jetAco","SingleJet150",1,"(-,-)"," "); 
		menu_2e30_v02->AddHlt("HLTJetMETRapidityGap","IsoEG10_Jet20_ForJet10",1,"20"," "); //100; 
		menu_2e30_v02->AddHlt("AptHLT4jet30","QuadJet15",1E1,"30"," "); 
		menu_2e30_v02->AddHlt("AptHLTXRelMu3_3jet30","Mu3_TripleJet20",1,"(3,30)"," "); 
		menu_2e30_v02->AddHlt("AptHLTXRelEl5_3jet30","EG5_TripleJet20",1,"(5,30)"," "); 

		menu_2e30_v02->AddHlt("AptHLT1Electron5_Open","SingleEG5",1E1,"5"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Electron8_DefThr","SingleEG5",1,"8"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Electron8_Open","SingleEG5",1,"8"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Electron10_Open","SingleEG8",1,"10"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Electron7_Open","SingleEG5",1E4,"7"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Electron9_Open","SingleEG5",1E4,"9"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Electron11_Open","SingleEG8",1E4,"11"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Electron13_Open","SingleEG8",1E4,"13"," "); 

		menu_2e30_v02->AddHlt("AptHLT1Electron15_Open","SingleEG10",1,"15"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Electron20_Open","SingleEG15",1,"20"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Electron25_Open","SingleEG15",1,"25"," "); 
		menu_2e30_v02->AddHlt("AptHLT2Electron5_Open","DoubleEG5",1,"5"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Electron30_Open","SingleEG20",1,"30"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1ElectronIso8_10_DefThr","SingleIsoEG8",1E1,"10"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1ElectronIso8_12_DefThr","SingleIsoEG8",1,"12"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1ElectronIso8_10_DblThr","SingleIsoEG8",1E1,"10"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1ElectronIso8_12_DblThr","SingleIsoEG8",1,"12"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1ElectronIso8_15_DblThr","SingleIsoEG8",1,"15"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Electron12_15_DefThr","SingleEG12",1,"15"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Electron12_15_DblThr","SingleEG12",1,"15"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Electron12_17_DblThr","SingleEG12",1,"17"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2ElectronIso8_10_Open","DoubleIsoEG8",1,"10"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2Electron10_12_Open","DoubleEG10",1,"12"," "); 
		menu_2e30_v02->AddHlt("HLT1Electron","SingleIsoEG12",1,"15"," "); 
		menu_2e30_v02->AddHlt("HLT1ElectronRelaxed","SingleEG15",1,"18"," "); 
		menu_2e30_v02->AddHlt("HLT2Electron","DoubleIsoEG8",1,"10"," "); 
		menu_2e30_v02->AddHlt("HLT2ElectronRelaxed","DoubleEG10",1,"12"," "); 

		//menu_2e30_v02->AddHlt("AptHLT1Photon5NoTrkIso","SingleEG5",5E3,"5"," ");
		menu_2e30_v02->AddHlt("AptHLT1Photon10_DefThr","SingleEG8",1E1,"10"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Photon10NoTrkIso","SingleEG8",1E2,"10"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon15NoTrkIso","SingleEG10",1E2,"15"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Photon20NoTrkIso","SingleEG10",1,"20"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Photon25NoTrkIso","SingleEG15",1,"25"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Photon30NoTrkIso","SingleEG20",1,"30"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon35NoTrkIso","SingleEG20",1E2,"35"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon40NoTrkIso","SingleEG25",2,"40"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon5_Open","SingleEG5",10,"5"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon10_Open","SingleEG8",1,"10"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon12NoTrkIso","SingleEG8",1,"12"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon12_Open","SingleEG8",1,"12"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Photon15_Open","SingleEG10",1,"15"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Photon25_Open","SingleEG10",1,"25"," "); 
		menu_2e30_v02->AddHlt("AptHLT1Photon30_Open","SingleEG20",1,"30"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2Photon15_Open","DoubleEG5",1,"15"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1PhotonIso8_15_DefThr","SingleIsoEG8",1,"15"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1PhotonIso8_17_DefThr","SingleIsoEG8",1,"17"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1PhotonIso8_20_DefThr","SingleIsoEG8",1,"20"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1PhotonIso8_25_DefThr","SingleIsoEG8",1,"25"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1PhotonIso8_30_DefThr","SingleIsoEG8",1,"30"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon12_20_DefThr","SingleEG12",1,"20"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon12_30_DefThr","SingleEG12",1,"30"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon12_35_DefThr","SingleEG12",1,"35"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon12_40_DefThr","SingleEG12",1,"40"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1PhotonIso8_20_DblThr","SingleIsoEG8",1,"20"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1PhotonIso8_25_DblThr","SingleIsoEG8",1,"25"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1PhotonIso8_30_DblThr","SingleIsoEG8",1,"30"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon12_15_DblThr","SingleEG12",5,"15"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon12_25_DefThr","SingleEG12",1,"25"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon12_25_DblThr","SingleEG12",1,"25"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon12_30_DblThr","SingleEG12",1,"30"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon12_35_DblThr","SingleEG12",1,"35"," "); 
		//menu_2e30_v02->AddHlt("AptHLT1Photon12_40_DblThr","SingleEG12",1,"40"," "); 
		menu_2e30_v02->AddHlt("AptHLT2Photon10NoTrkIso","DoubleEG5",1,"(10,10)"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2PhotonIso8_15_Open","DoubleIsoEG8",1,"15"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2Photon10_20_Open","DoubleEG10",1,"20"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2Photon10_20NoTrkIso","DoubleEG10",1,"(20,20)"," "); 
		//menu_2e30_v02->AddHlt("AptHLT2Photon10_30NoTrkIso","DoubleEG10",1,"(30,30)"," "); 
		menu_2e30_v02->AddHlt("HLT1Photon","SingleIsoEG12",1,"30"," "); 
		menu_2e30_v02->AddHlt("HLT1PhotonRelaxed","SingleEG15",1,"40"," "); 
		menu_2e30_v02->AddHlt("HLT2Photon","DoubleIsoEG8",1,"(20,20)"," "); 
		menu_2e30_v02->AddHlt("HLT2PhotonRelaxed","DoubleEG10",1,"(20,20)"," "); 
		menu_2e30_v02->AddHlt("HLT1EMHighEt","SingleEG15",1,"80"," "); 
		menu_2e30_v02->AddHlt("HLT1EMVeryHighEt","SingleEG15",1,"200"," "); 
		menu_2e30_v02->AddHlt("HLT2ElectronZCounter","DoubleIsoEG8",1,"(10,10)"," "); 
		menu_2e30_v02->AddHlt("HLT2ElectronExclusive","ExclusiveDoubleIsoEG6",1,"(6,6)"," "); 
		menu_2e30_v02->AddHlt("HLT2PhotonExclusive","ExclusiveDoubleIsoEG6",1,"(10,10)"," "); 
		menu_2e30_v02->AddHlt("HLT1PhotonL1Isolated","SingleIsoEG10",1E1,"12"," "); 

		menu_2e30_v02->AddHlt("HLTB1Jet","\\Diamond ",1,"180"," "); 
		menu_2e30_v02->AddHlt("HLTB2Jet","\\Diamond ",1,"120"," "); 
		menu_2e30_v02->AddHlt("HLTB3Jet","\\Diamond ",1,"70"," "); 
		menu_2e30_v02->AddHlt("HLTB4Jet","\\Diamond ",1,"40"," "); 
		menu_2e30_v02->AddHlt("HLTBHT","\\Diamond ",1,"470"," "); 
		menu_2e30_v02->AddHlt("HLT1Tau","SingleTauJet80",1,"15"," "); 
		menu_2e30_v02->AddHlt("HLT1Tau1MET","TauJet30_ETM30",1,"15"," "); 
		menu_2e30_v02->AddHlt("HLT2TauPixel","TauJet40",1,"15"," "); 
		menu_2e30_v02->AddHlt("HLTXElectronBJet","IsoEG10_Jet20",1,"(10,35)"," "); 
		menu_2e30_v02->AddHlt("HLTXElectron1Jet","IsoEG10_Jet30",1,"(12,40)"," "); 
		menu_2e30_v02->AddHlt("HLTXElectron2Jet","IsoEG10_Jet30",1,"(12,80)"," "); 
		menu_2e30_v02->AddHlt("HLTXElectron3Jet","IsoEG10_Jet30",1,"(12,60)"," "); 
		menu_2e30_v02->AddHlt("HLTXElectron4Jet","IsoEG10_Jet30",1,"(12,35)"," "); 
		menu_2e30_v02->AddHlt("HLTXElectronTau","IsoEG10_TauJet20",1,"(12,20)"," "); 

		//menu_2e30_v02->AddHlt("HLTHcalIsolatedTrack","NA",1,"NA"," "); 
		menu_2e30_v02->AddHlt("HLTHcalIsolatedTrackNoEcalIsol"," ",1," "," ");  
		menu_2e30_v02->AddHlt("HLTMinBiasPixel","ZeroBias",1,"- "," "); 
		//menu_2e30_v02->AddHlt("HLTMinBiasForAlignment"," ",1," "," "); 
		menu_2e30_v02->AddHlt("HLTMinBias","MinBias_HTT10",1,"-"," "); 
		menu_2e30_v02->AddHlt("HLTZeroBias","ZeroBias",1,"-"," "); 

}

void BookOHltMenu_Validate(OHltMenu*  menu_validate, double &iLumi, double &nBunches) {

		iLumi = 2E30;
		nBunches = 43;

		//menu_validate->AddHlt("AptHLT1MuonLevel1Open","SingleMuOpen",1.5E2,"-"," "); //1000; 
		//menu_validate->AddHlt("AptHLT1MuonLevel1","SingleMu7,DoubleMu3",2E1,"-"," "); //1000; 
		//menu_validate->AddHlt("AptHLT1MuonIso3_DefThr","SingleMu3",1,"3"," "); 
		//menu_validate->AddHlt("AptHLT1MuonIso5_DefThr","SingleMu5",1,"5"," "); 
		//menu_validate->AddHlt("AptHLT1MuonIso7_DefThr","SingleMu7",1,"7"," "); 
		//menu_validate->AddHlt("AptHLT1MuonIso9_DefThr","SingleMu7",1,"9"," "); 
		//menu_validate->AddHlt("AptHLT1MuonIso11_DefThr","SingleMu7",1,"11"," "); 
		//menu_validate->AddHlt("AptHLT1Muon3_DefThr","SingleMu3",1,"3"," "); 
		//menu_validate->AddHlt("AptHLT1Muon5_DefThr","SingleMu5",1,"5"," "); 
		//menu_validate->AddHlt("AptHLT1Muon7_DefThr","SingleMu7",1,"7"," "); 
		//menu_validate->AddHlt("AptHLT1Muon9_DefThr","SingleMu7",1,"9"," "); 
		//menu_validate->AddHlt("AptHLT1Muon11_DefThr","SingleMu7",1,"11"," "); 
		//menu_validate->AddHlt("AptHLT1Muon13_DefThr","SingleMu7",1,"13"," "); 
		//menu_validate->AddHlt("AptHLT1Muon15_DefThr","SingleMu7",1,"15"," "); 
		//menu_validate->AddHlt("AptHLT1MuonIso3","SingleMu3",1E6,"3"," "); 
		//menu_validate->AddHlt("AptHLT1MuonIso5","SingleMu5",1E6,"5"," "); 
		//menu_validate->AddHlt("AptHLT1MuonIso7","SingleMu7",1E6,"7"," "); 
		//menu_validate->AddHlt("AptHLT1MuonIso9","SingleMu7",1,"9"," "); //5E2; 
		//menu_validate->AddHlt("AptHLT1MuonIso11","SingleMu7",1,"11"," "); 
		//menu_validate->AddHlt("AptHLT1Muon3","SingleMu3",10,"3"," "); 
		//menu_validate->AddHlt("AptHLT1Muon5","SingleMu5",1,"5"," "); 
		//menu_validate->AddHlt("AptHLT1Muon7","SingleMu7",1,"7"," "); 
		//menu_validate->AddHlt("AptHLT1Muon9","SingleMu7",1,"9"," "); 
		//menu_validate->AddHlt("AptHLT1Muon11","SingleMu7",1,"11"," "); 
		//menu_validate->AddHlt("AptHLT1Muon13","SingleMu7",1,"13"," "); 
		//menu_validate->AddHlt("AptHLT1Muon15","SingleMu7",1,"15"," "); 
		//menu_validate->AddHlt("AptHLT1L2Muon11","SingleMu7",1,"11"," "); 
		////menu_validate->AddHlt("AptHLT1L2Muon16","SingleMu7",1,"16"," "); 
		//menu_validate->AddHlt("AptHLT2Muon3","DobuleMu3",1,"(3,3)"," "); 
		//menu_validate->AddHlt("HLT1MuonLevel1","SingleMu0",1,"-"," "); 
		menu_validate->AddHlt("HLT1MuonIso","SingleMu7",1,"11"," "); 
		menu_validate->AddHlt("OpenHLT1MuonIso","SingleMu7",1,"11"," "); 
		menu_validate->AddHlt("HLT1Electron","SingleIsoEG12",1,"15"," "); 
		menu_validate->AddHlt("OpenHLT1Electron","SingleIsoEG12",1,"15"," "); 
		/*
		menu_validate->AddHlt("HLT1MuonNonIso","SingleMu7",1,"16"," "); 
		menu_validate->AddHlt("HLT2MuonIso","DoubleMu3",1,"(3,3)"," "); 
		menu_validate->AddHlt("HLT2MuonNonIso","DoubleMu3",1,"(3,3)"," "); 
		menu_validate->AddHlt("HLT2MuonJPsi","DoubleMu3",1,"(3,3)"," "); 
		menu_validate->AddHlt("HLT2MuonUpsilon","DoubleMu3",1,"(3,3)"," "); 
		menu_validate->AddHlt("HLT2MuonZ","DoubleMu3",1,"(7,7)"," "); 
		menu_validate->AddHlt("HLTNMuonNonIso","DoubleMu3",1,"(3,3,3) "," "); 
		menu_validate->AddHlt("HLT2MuonSameSign","DoubleMu3",1,"(3,3)"," "); 
		//menu_validate->AddHlt("HLT1MuonPrescalePt3","SingleMu3",1E6,"3"," "); 
		//menu_validate->AddHlt("HLT1MuonPrescalePt5","SingleMu5",1E6,"5"," "); 
		//menu_validate->AddHlt("HLT1MuonPrescalePt7x7","SingleMu7",1E6,"7"," "); 
		//menu_validate->AddHlt("HLT1MuonPrescalePt7x10","SingleMu7",1E6,"7,10"," "); 
		menu_validate->AddHlt("HLTB1JetMu","Mu5_Jet15",1,"20 "," "); 
		menu_validate->AddHlt("HLTB2JetMu","Mu5_Jet15",1,"120"," "); 
		menu_validate->AddHlt("HLTB3JetMu","Mu5_Jet15",1,"70"," "); 
		menu_validate->AddHlt("HLTB4JetMu","Mu5_Jet15",1,"40"," "); 
		menu_validate->AddHlt("HLTBHTMu","HTT250",1,"300"," "); 
		menu_validate->AddHlt("HLTBJPsiMuMu","DoubleMu3",1,"(4,4)"," "); 
		menu_validate->AddHlt("HLTXMuonBJet","Mu5_Jet15",1,"(7,35)"," "); 
		menu_validate->AddHlt("HLTXMuonBJetSoftMuon","Mu5_Jet15",1,"(7,20)"," "); 
		menu_validate->AddHlt("HLTXMuonJets","Mu5_Jet15",1,"(7,40)"," "); 
		menu_validate->AddHlt("HLTXElectronMuon","\\star",1,"(8,7)"," "); 
		menu_validate->AddHlt("HLTXElectronMuonRelaxed","\\star",1,"(10,10)"," "); 
		menu_validate->AddHlt("HLTXMuonTau","Mu5_TauJet20",1,"(15,20)"," "); 

		//menu_validate->AddHlt("AptHLT1Level1jet15","SingleJet15",1E4,"-"," "); 
		menu_validate->AddHlt("AptHLT1jet30","SingleJet15",2E2,"30"," "); 
		menu_validate->AddHlt("AptHLT1jet40","SingleJet15",1E2,"40"," "); 
		menu_validate->AddHlt("AptHLT1jet50","SingleJet30",2E1,"50"," "); 
		menu_validate->AddHlt("AptHLT1jet60","SingleJet30",1E1,"60"," "); 
		menu_validate->AddHlt("AptHLT1jet80","SingleJet50",1,"80"," "); 
		menu_validate->AddHlt("AptHLT1jet110","SingleJet70",1,"110"," "); 
		menu_validate->AddHlt("AptHLT1jet150","SingleJet70",1,"150"," "); 
		//menu_validate->AddHlt("AptHLT1jet180","SingleJet70",1,"180"," "); 
		//menu_validate->AddHlt("AptHLT1jet250","SingleJet70",1,"250"," "); 
		menu_validate->AddHlt("AptHLT2jet100","SingleJet100, DoubleJet70",1,"(100,100)"," "); 
		menu_validate->AddHlt("HLT1jet","SingleJet150",1,"200"," "); 
		menu_validate->AddHlt("HLT2jet","SingleJet150, DoubleJet70",1,"150"," "); 
		menu_validate->AddHlt("HLT3jet","\\dagger",1,"85"," "); 
		menu_validate->AddHlt("HLT4jet","\\ddagger",1,"60"," "); 
		menu_validate->AddHlt("HLT1MET","ETM40",1,"65"," "); 
		menu_validate->AddHlt("HLT2jetAco","SingleJet150, DoubleJet70",1,"125"," "); 
		menu_validate->AddHlt("HLT1jet1METAco","SingleJet150",1,"(100,60)"," "); 
		menu_validate->AddHlt("HLT1jet1MET","SingleJet150",1,"(180,60)"," "); 
		menu_validate->AddHlt("HLT2jet1MET","SingleJet150",1,"(125,60)"," "); 
		menu_validate->AddHlt("HLT3jet1MET","SingleJet150",1,"(60,60)"," "); 
		menu_validate->AddHlt("HLT4jet1MET","SingleJet150",1,"(35,60)"," "); 
		menu_validate->AddHlt("HLT1MET1HT","HTT300 ",1,"(350,65)"," "); 
		menu_validate->AddHlt("HLT1SumET","ETT60",1,"120"," "); 
		//menu_validate->AddHlt("HLT1jetPE7","SingleJet15",1,"30"," "); 
		//menu_validate->AddHlt("HLT1jetPE5","SingleJet30",1,"60"," "); //10000; 
		//menu_validate->AddHlt("HLT1jetPE3","SingleJet70 ",1,"110"," "); //100; 
		//menu_validate->AddHlt("HLT1jetPE1","SingleJet100 ",1,"150"," "); 
		menu_validate->AddHlt("HLT1METPre3","MinBians_HTT10",1,"15"," "); 
		menu_validate->AddHlt("HLT1METPre2","MinBians_HTT10",1,"20"," "); 
		menu_validate->AddHlt("HLT1METPre1","ETM20",1E2,"30"," "); //100; 
		//menu_validate->AddHlt("AptHLT1MET20","ETM20",1E3,"20"," "); //100; 
		//menu_validate->AddHlt("AptHLT1MET35","ETM30",1E1,"35"," "); //100; 
		//menu_validate->AddHlt("AptHLT1MET50","ETM40",1,"50"," "); //100; 
		//menu_validate->AddHlt("AptHLT1MET65","ETM40",1,"65"," "); //100; 
		//menu_validate->AddHlt("AptHLT1MET75","ETM40",1,"75"," "); //100; 
		//menu_validate->AddHlt("HLT2jetAve30","SingleJet15",1,"30"," "); 
		//menu_validate->AddHlt("HLT2jetAve60","SingleJet30",1,"60"," "); 
		//menu_validate->AddHlt("HLT2jetAve110","SingleJet70",1,"110"," "); 
		//menu_validate->AddHlt("HLT2jetAve150","SingleJet100",1,"150"," "); 
		//menu_validate->AddHlt("HLT2jetAve200","SingleJet150",1,"200"," "); 
		//menu_validate->AddHlt("AptHLT2jetAve15","SingleJet15",4E2,"15"," "); 
		//menu_validate->AddHlt("AptHLT2jetAve30","SingleJet30",2E1,"30"," "); 
		//menu_validate->AddHlt("AptHLT2jetAve50","SingleJet50",4E1,"50"," "); 
		//menu_validate->AddHlt("AptHLT2jetAve70","SingleJet70",1,"70"," "); 
		//menu_validate->AddHlt("AptHLT2jetAve130","SingleJet130",1,"130"," "); 
		//menu_validate->AddHlt("AptHLT2jetAve220","SingleJet220",1,"220"," "); 
		menu_validate->AddHlt("HLT2jetvbfMET","ETM30",1,"(40,60)"," "); 
		menu_validate->AddHlt("HLTS2jet1METNV","SingleJet150",1,"(-,60)"," "); 
		menu_validate->AddHlt("HLTS2jet1METAco","SingleJet150",1,"(-,70)"," "); 
		menu_validate->AddHlt("HLTSjet1MET1Aco","SingleJet150",1,"(-,70)"," "); 
		menu_validate->AddHlt("HLTSjet2MET1Aco","SingleJet150",1,"(-,70)"," "); 
		menu_validate->AddHlt("HLTS2jetAco","SingleJet150",1,"(-,-)"," "); 
		menu_validate->AddHlt("HLTJetMETRapidityGap","IsoEG10_Jet20_ForJet10",1,"20"," "); //100; 
		menu_validate->AddHlt("AptHLT4jet30","QuadJet15",1E1,"30"," "); 
		menu_validate->AddHlt("AptHLTXRelMu3_3jet30","Mu3_TripleJet20",1,"(3,30)"," "); 
		menu_validate->AddHlt("AptHLTXRelEl5_3jet30","EG5_TripleJet20",1,"(5,30)"," "); 

		menu_validate->AddHlt("AptHLT1Electron5_Open","SingleEG5",1E1,"5"," "); 
		menu_validate->AddHlt("AptHLT1Electron8_DefThr","SingleEG5",1,"8"," "); 
		menu_validate->AddHlt("AptHLT1Electron8_Open","SingleEG5",1,"8"," "); 
		menu_validate->AddHlt("AptHLT1Electron10_Open","SingleEG8",1,"10"," "); 
		//menu_validate->AddHlt("AptHLT1Electron7_Open","SingleEG5",1E4,"7"," "); 
		//menu_validate->AddHlt("AptHLT1Electron9_Open","SingleEG5",1E4,"9"," "); 
		//menu_validate->AddHlt("AptHLT1Electron11_Open","SingleEG8",1E4,"11"," "); 
		//menu_validate->AddHlt("AptHLT1Electron13_Open","SingleEG8",1E4,"13"," "); 

		menu_validate->AddHlt("AptHLT1Electron15_Open","SingleEG10",1,"15"," "); 
		menu_validate->AddHlt("AptHLT1Electron20_Open","SingleEG15",1,"20"," "); 
		menu_validate->AddHlt("AptHLT1Electron25_Open","SingleEG15",1,"25"," "); 
		menu_validate->AddHlt("AptHLT2Electron5_Open","DoubleEG5",1,"5"," "); 
		//menu_validate->AddHlt("AptHLT1Electron30_Open","SingleEG20",1,"30"," "); 
		//menu_validate->AddHlt("AptHLT1ElectronIso8_10_DefThr","SingleIsoEG8",1E1,"10"," "); 
		//menu_validate->AddHlt("AptHLT1ElectronIso8_12_DefThr","SingleIsoEG8",1,"12"," "); 
		//menu_validate->AddHlt("AptHLT1ElectronIso8_10_DblThr","SingleIsoEG8",1E1,"10"," "); 
		//menu_validate->AddHlt("AptHLT1ElectronIso8_12_DblThr","SingleIsoEG8",1,"12"," "); 
		//menu_validate->AddHlt("AptHLT1ElectronIso8_15_DblThr","SingleIsoEG8",1,"15"," "); 
		//menu_validate->AddHlt("AptHLT1Electron12_15_DefThr","SingleEG12",1,"15"," "); 
		//menu_validate->AddHlt("AptHLT1Electron12_15_DblThr","SingleEG12",1,"15"," "); 
		//menu_validate->AddHlt("AptHLT1Electron12_17_DblThr","SingleEG12",1,"17"," "); 
		//menu_validate->AddHlt("AptHLT2ElectronIso8_10_Open","DoubleIsoEG8",1,"10"," "); 
		//menu_validate->AddHlt("AptHLT2Electron10_12_Open","DoubleEG10",1,"12"," "); 
		menu_validate->AddHlt("HLT1Electron","SingleIsoEG12",1,"15"," "); 
		menu_validate->AddHlt("HLT1ElectronRelaxed","SingleEG15",1,"18"," "); 
		menu_validate->AddHlt("HLT2Electron","DoubleIsoEG8",1,"10"," "); 
		menu_validate->AddHlt("HLT2ElectronRelaxed","DoubleEG10",1,"12"," "); 

		//menu_validate->AddHlt("AptHLT1Photon5NoTrkIso","SingleEG5",5E3,"5"," ");
		menu_validate->AddHlt("AptHLT1Photon10_DefThr","SingleEG8",1E1,"10"," "); 
		menu_validate->AddHlt("AptHLT1Photon10NoTrkIso","SingleEG8",1E2,"10"," "); 
		//menu_validate->AddHlt("AptHLT1Photon15NoTrkIso","SingleEG10",1E2,"15"," "); 
		menu_validate->AddHlt("AptHLT1Photon20NoTrkIso","SingleEG10",1,"20"," "); 
		menu_validate->AddHlt("AptHLT1Photon25NoTrkIso","SingleEG15",1,"25"," "); 
		menu_validate->AddHlt("AptHLT1Photon30NoTrkIso","SingleEG20",1,"30"," "); 
		//menu_validate->AddHlt("AptHLT1Photon35NoTrkIso","SingleEG20",1E2,"35"," "); 
		//menu_validate->AddHlt("AptHLT1Photon40NoTrkIso","SingleEG25",2,"40"," "); 
		//menu_validate->AddHlt("AptHLT1Photon5_Open","SingleEG5",10,"5"," "); 
		//menu_validate->AddHlt("AptHLT1Photon10_Open","SingleEG8",1,"10"," "); 
		//menu_validate->AddHlt("AptHLT1Photon12NoTrkIso","SingleEG8",1,"12"," "); 
		//menu_validate->AddHlt("AptHLT1Photon12_Open","SingleEG8",1,"12"," "); 
		menu_validate->AddHlt("AptHLT1Photon15_Open","SingleEG10",1,"15"," "); 
		menu_validate->AddHlt("AptHLT1Photon25_Open","SingleEG10",1,"25"," "); 
		menu_validate->AddHlt("AptHLT1Photon30_Open","SingleEG20",1,"30"," "); 
		//menu_validate->AddHlt("AptHLT2Photon15_Open","DoubleEG5",1,"15"," "); 
		//menu_validate->AddHlt("AptHLT1PhotonIso8_15_DefThr","SingleIsoEG8",1,"15"," "); 
		//menu_validate->AddHlt("AptHLT1PhotonIso8_17_DefThr","SingleIsoEG8",1,"17"," "); 
		//menu_validate->AddHlt("AptHLT1PhotonIso8_20_DefThr","SingleIsoEG8",1,"20"," "); 
		//menu_validate->AddHlt("AptHLT1PhotonIso8_25_DefThr","SingleIsoEG8",1,"25"," "); 
		//menu_validate->AddHlt("AptHLT1PhotonIso8_30_DefThr","SingleIsoEG8",1,"30"," "); 
		//menu_validate->AddHlt("AptHLT1Photon12_20_DefThr","SingleEG12",1,"20"," "); 
		//menu_validate->AddHlt("AptHLT1Photon12_30_DefThr","SingleEG12",1,"30"," "); 
		//menu_validate->AddHlt("AptHLT1Photon12_35_DefThr","SingleEG12",1,"35"," "); 
		//menu_validate->AddHlt("AptHLT1Photon12_40_DefThr","SingleEG12",1,"40"," "); 
		//menu_validate->AddHlt("AptHLT1PhotonIso8_20_DblThr","SingleIsoEG8",1,"20"," "); 
		//menu_validate->AddHlt("AptHLT1PhotonIso8_25_DblThr","SingleIsoEG8",1,"25"," "); 
		//menu_validate->AddHlt("AptHLT1PhotonIso8_30_DblThr","SingleIsoEG8",1,"30"," "); 
		//menu_validate->AddHlt("AptHLT1Photon12_15_DblThr","SingleEG12",5,"15"," "); 
		//menu_validate->AddHlt("AptHLT1Photon12_25_DefThr","SingleEG12",1,"25"," "); 
		//menu_validate->AddHlt("AptHLT1Photon12_25_DblThr","SingleEG12",1,"25"," "); 
		//menu_validate->AddHlt("AptHLT1Photon12_30_DblThr","SingleEG12",1,"30"," "); 
		//menu_validate->AddHlt("AptHLT1Photon12_35_DblThr","SingleEG12",1,"35"," "); 
		//menu_validate->AddHlt("AptHLT1Photon12_40_DblThr","SingleEG12",1,"40"," "); 
		menu_validate->AddHlt("AptHLT2Photon10NoTrkIso","DoubleEG5",1,"(10,10)"," "); 
		//menu_validate->AddHlt("AptHLT2PhotonIso8_15_Open","DoubleIsoEG8",1,"15"," "); 
		//menu_validate->AddHlt("AptHLT2Photon10_20_Open","DoubleEG10",1,"20"," "); 
		//menu_validate->AddHlt("AptHLT2Photon10_20NoTrkIso","DoubleEG10",1,"(20,20)"," "); 
		//menu_validate->AddHlt("AptHLT2Photon10_30NoTrkIso","DoubleEG10",1,"(30,30)"," "); 
		menu_validate->AddHlt("HLT1Photon","SingleIsoEG12",1,"30"," "); 
		menu_validate->AddHlt("HLT1PhotonRelaxed","SingleEG15",1,"40"," "); 
		menu_validate->AddHlt("HLT2Photon","DoubleIsoEG8",1,"(20,20)"," "); 
		menu_validate->AddHlt("HLT2PhotonRelaxed","DoubleEG10",1,"(20,20)"," "); 
		menu_validate->AddHlt("HLT1EMHighEt","SingleEG15",1,"80"," "); 
		menu_validate->AddHlt("HLT1EMVeryHighEt","SingleEG15",1,"200"," "); 
		menu_validate->AddHlt("HLT2ElectronZCounter","DoubleIsoEG8",1,"(10,10)"," "); 
		menu_validate->AddHlt("HLT2ElectronExclusive","ExclusiveDoubleIsoEG6",1,"(6,6)"," "); 
		menu_validate->AddHlt("HLT2PhotonExclusive","ExclusiveDoubleIsoEG6",1,"(10,10)"," "); 
		menu_validate->AddHlt("HLT1PhotonL1Isolated","SingleIsoEG10",1E1,"12"," "); 

		menu_validate->AddHlt("HLTB1Jet","\\Diamond ",1,"180"," "); 
		menu_validate->AddHlt("HLTB2Jet","\\Diamond ",1,"120"," "); 
		menu_validate->AddHlt("HLTB3Jet","\\Diamond ",1,"70"," "); 
		menu_validate->AddHlt("HLTB4Jet","\\Diamond ",1,"40"," "); 
		menu_validate->AddHlt("HLTBHT","\\Diamond ",1,"470"," "); 
		menu_validate->AddHlt("HLT1Tau","SingleTauJet80",1,"15"," "); 
		menu_validate->AddHlt("HLT1Tau1MET","TauJet30_ETM30",1,"15"," "); 
		menu_validate->AddHlt("HLT2TauPixel","TauJet40",1,"15"," "); 
		menu_validate->AddHlt("HLTXElectronBJet","IsoEG10_Jet20",1,"(10,35)"," "); 
		menu_validate->AddHlt("HLTXElectron1Jet","IsoEG10_Jet30",1,"(12,40)"," "); 
		menu_validate->AddHlt("HLTXElectron2Jet","IsoEG10_Jet30",1,"(12,80)"," "); 
		menu_validate->AddHlt("HLTXElectron3Jet","IsoEG10_Jet30",1,"(12,60)"," "); 
		menu_validate->AddHlt("HLTXElectron4Jet","IsoEG10_Jet30",1,"(12,35)"," "); 
		menu_validate->AddHlt("HLTXElectronTau","IsoEG10_TauJet20",1,"(12,20)"," "); 

		//menu_validate->AddHlt("HLTHcalIsolatedTrack","NA",1,"NA"," "); 
		menu_validate->AddHlt("HLTHcalIsolatedTrackNoEcalIsol"," ",1," "," ");  
		menu_validate->AddHlt("HLTMinBiasPixel","ZeroBias",1,"- "," "); 
		//menu_validate->AddHlt("HLTMinBiasForAlignment"," ",1," "," "); 
		menu_validate->AddHlt("HLTMinBias","MinBias_HTT10",1,"-"," "); 
		menu_validate->AddHlt("HLTZeroBias","ZeroBias",1,"-"," "); 
		*/

}
