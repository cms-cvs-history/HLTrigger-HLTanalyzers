/////////////////////////////////////////////////////////////////////////////////////////////////
//
//          Macro calculating overlaps between different triggers, L1 trigger individual- and
//          pure-rates, taking into account background and signal type samples
//
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include "L1Tree.h"

#include "TH1.h"
#include "TChain.h"
#include "TCut.h"

#include <map>

using namespace std;

Double_t eff(Int_t a, Int_t b);
Double_t seff(Int_t a, Int_t b);
Double_t eff(Double_t a, Double_t b);
Double_t seff(Double_t a, Double_t b);

void MakeL1Menu_1E27_1(double &ILumi, double &nFilledBunches, 
		       vector<string> &trignames, vector<int> &prescales, int Version);
void MakeL1Menu_38E29_43(double &ILumi, double &nFilledBunches, 
			 vector<string> &trignames, vector<int> &prescales, int Version);
void MakeL1Menu_17E30_43(double &ILumi, double &nFilledBunches, 
			 vector<string> &trignames, vector<int> &prescales, int Version);
void MakeL1Menu_10E32_1000(double &ILumi, double &nFilledBunches, 
			 vector<string> &trignames, vector<int> &prescales, int Version);


void MakeL1Menu_61E30_43(double &ILumi, double &nFilledBunches, 
			 vector<string> &trignames, vector<int> &prescales, int Version);
void MakeL1Menu_11E31_156(double &ILumi, double &nFilledBunches, 
			 vector<string> &trignames, vector<int> &prescales, int Version);
void MakeL1Menu_56E31_156(double &ILumi, double &nFilledBunches, 
			 vector<string> &trignames, vector<int> &prescales, int Version);




int main(int argc, char *argv[]){
  
  int NEntries = -1;
  if (argc>1) {
    NEntries = atoi(argv[1]);
  }
  int Menu = -1;
  if (argc>2) {
    Menu = atoi(argv[2]);
  }
  int Version = -1;
  if (argc>3) {
    Version = atoi(argv[3]);
  }

  float scale = 1000.; // convert to kHz
    
  ////////////////////////////////////////////////////////////
  // Instanteneous Luminosity [cm^-2 s^-1]
  
  // For accurate rate calculation
  const double bunchCrossingTime = 25.0E-09;  // 25 ns
  const double maxFilledBunches = 3557;
  
  double ILumi = 1.E27;
  double nFilledBunches = 1;

  ////////////////////////////////////////////////////////////
  // Trigger names and prescales.
  // The order is maintained during rate calculation. So are prescales. 
  // I.e., keep the order and assign to any triggername a prescale!
  vector<string> trignames; vector<int> prescales;

  // Select Trigger Menu
  if (Menu==0) MakeL1Menu_1E27_1(ILumi, nFilledBunches, trignames, prescales, Version);
  else if (Menu==1) MakeL1Menu_38E29_43(ILumi, nFilledBunches, trignames, prescales, Version);
  else if (Menu==2) MakeL1Menu_17E30_43(ILumi, nFilledBunches, trignames, prescales, Version);
  else if (Menu==3) MakeL1Menu_10E32_1000(ILumi, nFilledBunches, trignames, prescales, Version);
  else if (Menu==4) MakeL1Menu_61E30_43(ILumi, nFilledBunches, trignames, prescales, Version);
  else if (Menu==5) MakeL1Menu_11E31_156(ILumi, nFilledBunches, trignames, prescales, Version);
  else if (Menu==6) MakeL1Menu_56E31_156(ILumi, nFilledBunches, trignames, prescales, Version);
  //else MakeL1Menu_10E32_1000(ILumi, nFilledBunches, trignames, prescales);

  int Ntrig = (int)trignames.size();
  double collisionRate = (nFilledBunches / maxFilledBunches) / bunchCrossingTime ;  // Hz
  
  ////////////////////////////////////////////////////////////
  // Cross-sections [pb]
      
  vector<Double_t> xsec;
  xsec.push_back(7.923E10); // PYTHIA cross-section for MinBias (pb)
  xsec.push_back(7.923E10); // PYTHIA cross-section for MinBias (pb)
  xsec.push_back(7.923E10); // PYTHIA cross-section for MinBias (pb)
  //xsec.push_back(5.5E9); //Wrong!!!!: pp->muX (pb), for checking with old rates only!!!!!!!!!!!!!!!!!!

  vector<Double_t> skmeff; // Skim efficiencies
  skmeff.push_back(1.);
  skmeff.push_back(0.0008); // ppMuX filter efficiency
  skmeff.push_back(0.0097); // ppEleX filter efficiency

  for (unsigned int ip = 0; ip < skmeff.size(); ip++)
    {xsec[ip] *= skmeff[ip];}
  
  // Convert cross-sections to cm^2
  for (unsigned int i = 0; i < skmeff.size(); i++){xsec[i] *= 1.E-36;}
  
  ////////////////////////////////////////////////////////////
  // Files
  
  /**/
  vector<TChain*> TabChain;
  vector<bool> doMuonCut; vector<bool> doElecCut;
  vector<TString> ProcFil;

  TString SamplesDIR = "/uscms/home/chinhan/lpctau/CMSSW_1_6_0/src/UserCode/chinhan/Config/RatesAndPrescales/hlt";
  //TString SamplesDIR = "/scratch/bargassa";

  /**/
  ProcFil.push_back(SamplesDIR+"/mergedFiles_CSA07_MB-hltanal1_calotowers_1.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_CSA07_MB-hltanal1_calotowers_2.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_CSA07_MB-hltanal1_calotowers_3.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_CSA07_MB-hltanal1_calotowers_4.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_CSA07_MB-hltanal1_calotowers_5.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_CSA07_MB-hltanal1_calotowers_6.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_CSA07_MB-hltanal1_calotowers_7.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_CSA07_MB-hltanal1_calotowers_8.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_CSA07_MB-hltanal1_calotowers_9.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_CSA07_MB-hltanal1_calotowers_10.root");
  TabChain.push_back(new TChain("HltTree"));
  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
    TabChain[0]->Add(ProcFil[ipfile]);
  }
  doMuonCut.push_back(true); doElecCut.push_back(true);

  /**/
  ProcFil.clear();
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppMuXHLTanaL1_calotowers_1.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppMuXHLTanaL1_calotowers_2.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppMuXHLTanaL1_calotowers_3.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppMuXHLTanaL1_calotowers_4.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppMuXHLTanaL1_calotowers_5.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppMuXHLTanaL1_calotowers_6.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppMuXHLTanaL1_calotowers_7.root");
  TabChain.push_back(new TChain("HltTree"));
  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
    TabChain[1]->Add(ProcFil[ipfile]);
  }
  doMuonCut.push_back(false); doElecCut.push_back(true);

  ProcFil.clear();
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppEleXHLTanaL1_calotowers_1.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppEleXHLTanaL1_calotowers_2.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppEleXHLTanaL1_calotowers_3.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppEleXHLTanaL1_calotowers_4.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppEleXHLTanaL1_calotowers_5.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppEleXHLTanaL1_calotowers_6.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppEleXHLTanaL1_calotowers_7.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppEleXHLTanaL1_calotowers_8.root");
  ProcFil.push_back(SamplesDIR+"/mergedFiles_ppEleXHLTanaL1_calotowers_9.root");
  TabChain.push_back(new TChain("HltTree"));
  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){
    TabChain[2]->Add(ProcFil[ipfile]);
  }
  doMuonCut.push_back(true); doElecCut.push_back(false);
  /**/

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
  
  vector<L1Tree*> l1t;
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

    l1t.push_back(new L1Tree((TTree*)TabChain[ip],Ntrig));
    cout<<"Processing bin "<<ip<<" ( "<<l1t[ip]->fChain->GetEntries()<<" events ) "<<" ..."<<endl;
    int deno = NEntries; 
    if (NEntries <= 0) {
      deno = (int)l1t[ip]->fChain->GetEntries(); 
    }
    l1t[ip]->Loop(iCount,sPureCount,pureCount,overlapCount,trignames,prescales,
		  deno,doMuonCut[ip],doElecCut[ip]);
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
    /**/
    cout.setf(ios::floatfield,ios::fixed);
    cout<<setprecision(2);
    for (int it=0; it < Ntrig; it++){
      cout  << setw(3) << it << ")";
      cout  << setw(28) << trignames[it]  << " (";
      if (prescales[it]>=1000000) {
	cout  <<setw(5)<<scientific<<setprecision(0)<< (float)prescales[it]<<fixed;
      } else {
	cout  <<setw(7)<<prescales[it];
      }
      cout  <<  setprecision(2) << ")";
      cout<< " :   Indiv.: " << setw(5) << Rat_bin[it]/scale << " +- " << setw(5) << sqrt(sRat_bin[it])/scale 
	  << "   sPure: " << setw(5) << seqpRat_bin[it]/scale
	  << "   Pure: " << setw(5) << pRat_bin[it]/scale 
	  << "   Cumul: " << setw(6) << cRat_bin[it]/scale << "\n"<<flush;
    }
    cout << "\n"<<flush;
    cout << setw(60) << "TOTAL RATE : " << setw(5) << RTOT_bin/scale << " +- " << sRTOT_bin/scale << " kHz" << "\n";
    cout << "\n"<<flush;
    /**/
    
  }

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

  /**/
  cout << endl;
  cout << "Trigger global overlaps : " << endl;
  for (int it = 0; it != Ntrig; ++it){
    for (int jt = 0; jt != Ntrig; ++jt){
      if (jt>=it) {
	// Overlap O(ij) = T(i) x T(j) / T(j)
        cout << "i="  << setw(3)<< it << " j="  << setw(3)<< jt << "     "  << setw(6)<< eff((Onum.at(it))[jt],Odenp[jt]) << endl;   
      }
    }
  }
  /**/

  cout.setf(ios::floatfield,ios::fixed);
  cout<<setprecision(2);

  cout << "\n";
  cout << "Level-1 Trigger Rates [kHz] (prescale): " << "\n";
  cout << "--------------------------------------------------------------------------------------------------------------------\n"<<flush;
  // This is with the accurate formula: 
  for (int it=0; it < Ntrig; it++){
    cout  << setw(3) << it << ")";
    cout  << setw(28) << trignames[it]  << " (";
    if (prescales[it]>=1000000) {
      cout  <<setw(5)<<scientific<<setprecision(0)<< (float)prescales[it]<<fixed;
    } else {
      cout  <<setw(7)<<prescales[it];
    }
    cout  <<  setprecision(2) << ")" 
	  << " :   Indiv: " << setw(5) << Rat[it]/scale << " +- " << setw(5) << sqrt(sRat[it])/scale 
	  << "   sPure: " << setw(5) << seqpRat[it]/scale
	  << "   Pure: " << setw(5) << pRat[it]/scale 
	  << "   Cumul: " << setw(6) << cRat[it]/scale << "\n"<<flush;
  }
  cout << "--------------------------------------------------------------------------------------------------------------------\n"<<flush;
  cout << setw(60) << "TOTAL RATE : " << setw(5) << RTOT/scale << " +- " << sRTOT/scale << " kHz" << "\n";

}

/**** Different Beam conditions: ****/
 
//const double ILumi = 1.E27;
//const double nFilledBunches = 1;

//const double ILumi = 3.8E29;
//const double nFilledBunches = 43;

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

//const double ILumi = 2.3E31;
//const double nFilledBunches = 936;

//const double ILumi = 1.0E32;
//const double nFilledBunches = 936;

//const double ILumi = 5.0E32;
//const double nFilledBunches = 936;

//const double ILumi = 1.7E32;
//const double nFilledBunches = 2808;

//const double ILumi = 1.0E33;
//const double nFilledBunches = 2808;

//const double ILumi = 1.0E34;
//const double nFilledBunches = 2808;


/****   ***   ****/
 
// Option 0
void MakeL1Menu_1E27_1(double &ILumi, double &nFilledBunches, vector<string> &trignames, vector<int> &prescales,
		       int Version = 1) {
  ILumi = 1.0E27;
  nFilledBunches = 1;

  /**/ // New Relaxed triggers for early run periods 
  trignames.push_back("SingleMu0"); prescales.push_back(1);
  trignames.push_back("SingleMu3"); prescales.push_back(1);	
  trignames.push_back("SingleMu5"); prescales.push_back(1);    
  trignames.push_back("SingleMu7"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu3"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu5"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu7"); prescales.push_back(1);	
  trignames.push_back("L1_DoubleMu3"); prescales.push_back(1);	
  trignames.push_back("L1_TripleMu3"); prescales.push_back(1);	
  trignames.push_back("L1_Mu3_Jet15"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_Jet15"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_Jet20"); prescales.push_back(1); 
  trignames.push_back("L1_Mu3_IsoEG5"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_IsoEG10"); prescales.push_back(1);
  trignames.push_back("L1_Mu3_EG12"); prescales.push_back(1);
  trignames.push_back("SingleEG5"); prescales.push_back(1);	
  trignames.push_back("SingleEG6"); prescales.push_back(1);	
  trignames.push_back("SingleEG7"); prescales.push_back(1);	
  trignames.push_back("SingleEG8"); prescales.push_back(1);	
  trignames.push_back("SingleEG9"); prescales.push_back(1);	
  trignames.push_back("SingleEG10"); prescales.push_back(1);	
  trignames.push_back("SingleEG12"); prescales.push_back(1);	
  trignames.push_back("SingleEG15"); prescales.push_back(1);	
  trignames.push_back("SingleIsoEG12"); prescales.push_back(1);	
  trignames.push_back("SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG8"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG12"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG15"); prescales.push_back(1);	
  //trignames.push_back("SingleJet10"); prescales.push_back(1);
  //trignames.push_back("SingleJet15"); prescales.push_back(1);
  trignames.push_back("SingleJet20"); prescales.push_back(1);
  trignames.push_back("SingleJet25"); prescales.push_back(1);
  trignames.push_back("SingleJet30"); prescales.push_back(1);
  trignames.push_back("SingleJet35"); prescales.push_back(1);
  trignames.push_back("SingleJet40"); prescales.push_back(1);
  trignames.push_back("SingleJet50"); prescales.push_back(1);
  trignames.push_back("SingleJet70"); prescales.push_back(1);
  trignames.push_back("SingleJet100"); prescales.push_back(1);
  trignames.push_back("SingleJet150"); prescales.push_back(1);
  //trignames.push_back("L1_SingleJet15"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet30"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet70"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet100"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet150"); prescales.push_back(1);
  //trignames.push_back("SingleTau10"); prescales.push_back(1);	
  trignames.push_back("SingleTau20"); prescales.push_back(1);	
  trignames.push_back("SingleTau25"); prescales.push_back(1);	
  trignames.push_back("SingleTau30"); prescales.push_back(1);	
  trignames.push_back("SingleTau35"); prescales.push_back(1);	
  trignames.push_back("SingleTau40"); prescales.push_back(1);	
  trignames.push_back("SingleTau50"); prescales.push_back(1);	
  trignames.push_back("SingleTau60"); prescales.push_back(1);	
  trignames.push_back("SingleTau70"); prescales.push_back(1);	
  trignames.push_back("SingleTau80"); prescales.push_back(1);	
  trignames.push_back("SingleTau100"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet80"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet100"); prescales.push_back(1);	
  trignames.push_back("L1_HTT200"); prescales.push_back(1);
  trignames.push_back("L1_HTT250"); prescales.push_back(1);
  trignames.push_back("L1_HTT300"); prescales.push_back(1);
  trignames.push_back("ETM20"); prescales.push_back(1);
  trignames.push_back("ETM25"); prescales.push_back(1);	
  trignames.push_back("ETM30"); prescales.push_back(1);	
  trignames.push_back("ETM35"); prescales.push_back(1);	
  trignames.push_back("ETM40"); prescales.push_back(1);	
  trignames.push_back("ETM50"); prescales.push_back(1);	
  trignames.push_back("ETM60"); prescales.push_back(1);	
  trignames.push_back("L1_ETM20"); prescales.push_back(1);
  trignames.push_back("L1_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_DoubleIsoEG8"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG10"); prescales.push_back(1);
  trignames.push_back("L1_DoubleJet70"); prescales.push_back(1);
  trignames.push_back("L1_DoubleJet100"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet40"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet15"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet30"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet70"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet20_ETM20"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_HTT100_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TripleJet50"); prescales.push_back(1);
  trignames.push_back("L1_QuadJet30"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveDoubleIsoEG6"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveDoubleJet60"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveJet25_Gap_Jet25"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet20_ForJet10"); prescales.push_back(1);
  trignames.push_back("L1_MinBias_HTT10"); prescales.push_back(1);
  trignames.push_back("L1_ZeroBias"); prescales.push_back(1);
  // New Minbias triggers
  trignames.push_back("MinBias_SingleHF1"); prescales.push_back(1);
  trignames.push_back("MinBias_DoubleHF1"); prescales.push_back(1);
  /**/

}


// Option 1
void MakeL1Menu_38E29_43(double &ILumi, double &nFilledBunches, vector<string> &trignames, vector<int> &prescales,
		       int Version = 1) {
  ILumi = 3.8E29;
  nFilledBunches = 43;

  /**/ // New Relaxed triggers for early run periods 
  if (Version==10) { trignames.push_back("OrAllMu"); prescales.push_back(1); }

  trignames.push_back("SingleMu0"); prescales.push_back(1);
  trignames.push_back("SingleMu3"); prescales.push_back(1);	
  trignames.push_back("SingleMu5"); prescales.push_back(1);    
  //trignames.push_back("SingleMu7"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleMu3"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleMu5"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu7"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu14"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu25"); prescales.push_back(1);	

  trignames.push_back("L1_DoubleMu3"); prescales.push_back(1);	
  trignames.push_back("L1_TripleMu3"); prescales.push_back(1);	

  trignames.push_back("L1_Mu3_Jet15"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_Jet15"); prescales.push_back(1); 
  trignames.push_back("L1_Mu3_Jet70"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_Jet20"); prescales.push_back(1);
  
  trignames.push_back("L1_Mu3_IsoEG5"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_IsoEG10"); prescales.push_back(1);
  trignames.push_back("L1_Mu3_EG12"); prescales.push_back(1);
  
  trignames.push_back("SingleEG5"); prescales.push_back(1);	
  //trignames.push_back("SingleEG6"); prescales.push_back(1);	
  //trignames.push_back("SingleEG7"); prescales.push_back(1);	
  trignames.push_back("SingleEG8"); prescales.push_back(1);	
  //trignames.push_back("SingleEG9"); prescales.push_back(1);	
  trignames.push_back("SingleEG10"); prescales.push_back(1);	
  trignames.push_back("SingleEG12"); prescales.push_back(1);	
  //trignames.push_back("SingleEG15"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG5"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG8"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG10"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG12"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG25"); prescales.push_back(1);	

  //trignames.push_back("SingleIsoEG12"); prescales.push_back(1);	
  //trignames.push_back("SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG8"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG12"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG25"); prescales.push_back(1);	

  if (Version==7) { trignames.push_back("SingleJet10"); prescales.push_back(1); }
  if (Version==7) { trignames.push_back("SingleJet15"); prescales.push_back(1); }
  if (Version==8) { trignames.push_back("SingleJet15"); prescales.push_back(1); }
  if (Version==9) { trignames.push_back("SingleJet15"); prescales.push_back(1); }
  if (Version==10) { trignames.push_back("SingleJet15"); prescales.push_back(1); }
  trignames.push_back("SingleJet20"); prescales.push_back(1);
  //trignames.push_back("SingleJet25"); prescales.push_back(1);
  trignames.push_back("SingleJet30"); prescales.push_back(1);
  //trignames.push_back("SingleJet35"); prescales.push_back(1);
  //trignames.push_back("SingleJet40"); prescales.push_back(1);
  trignames.push_back("SingleJet50"); prescales.push_back(1);
  trignames.push_back("SingleJet70"); prescales.push_back(1);
  //trignames.push_back("SingleJet100"); prescales.push_back(1);
  //trignames.push_back("SingleJet150"); prescales.push_back(1);
  if (Version==7) { trignames.push_back("L1_SingleJet15"); prescales.push_back(1);}
  if (Version==8) { trignames.push_back("L1_SingleJet15"); prescales.push_back(1);}
  //trignames.push_back("L1_SingleJet30"); prescales.push_back(1);
  //trignames.push_back("L1_SingleJet70"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet100"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet150"); prescales.push_back(1);
  
  //trignames.push_back("SingleTauJet10"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet20"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet25"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet30"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet35"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet40"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet50"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet60"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet70"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet80"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet100"); prescales.push_back(1);
  //trignames.push_back("L1_SingleTauJet10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet30"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet40"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet60"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleTauJet80"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet100"); prescales.push_back(1);
  
  trignames.push_back("L1_HTT200"); prescales.push_back(1);
  trignames.push_back("L1_HTT250"); prescales.push_back(1);
  trignames.push_back("L1_HTT300"); prescales.push_back(1);
  trignames.push_back("L1_HTT400"); prescales.push_back(1);
  //trignames.push_back("L1_HTT500"); prescales.push_back(1);

  trignames.push_back("ETM20"); prescales.push_back(1);
  trignames.push_back("ETM25"); prescales.push_back(1);	
  //trignames.push_back("ETM30"); prescales.push_back(1);	
  //trignames.push_back("ETM35"); prescales.push_back(1);	
  //trignames.push_back("ETM40"); prescales.push_back(1);	
  //trignames.push_back("ETM50"); prescales.push_back(1);	
  //trignames.push_back("ETM60"); prescales.push_back(1);	
  //trignames.push_back("L1_ETM20"); prescales.push_back(1);
  trignames.push_back("L1_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_ETM50"); prescales.push_back(1);
  trignames.push_back("L1_ETM60"); prescales.push_back(1);
  
  trignames.push_back("L1_DoubleIsoEG8"); prescales.push_back(1);
  trignames.push_back("L1_DoubleIsoEG10"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG5"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG10"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG15"); prescales.push_back(1);

  trignames.push_back("L1_DoubleJet70"); prescales.push_back(1);
  trignames.push_back("L1_DoubleJet100"); prescales.push_back(1);

  trignames.push_back("L1_DoubleTauJet20"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet30"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet35"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet40"); prescales.push_back(1);
  
  trignames.push_back("L1_IsoEG10_Jet15"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet30"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet70"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet20_ETM20"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_HTT100_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TripleJet50"); prescales.push_back(1);
  trignames.push_back("L1_QuadJet30"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveDoubleIsoEG6"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveDoubleJet60"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveJet25_Gap_Jet25"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet20_ForJet10"); prescales.push_back(1);
  trignames.push_back("L1_MinBias_HTT10"); prescales.push_back(1);
  trignames.push_back("L1_ZeroBias"); prescales.push_back(1);

  // New Minbias triggers
  if (Version==1) { trignames.push_back("MinBias_SingleHF1"); prescales.push_back(1);  }
  if (Version==1) { trignames.push_back("MinBias_DoubleHF1"); prescales.push_back(1);  }
  if (Version==2) { trignames.push_back("MinBias_SingleHF1"); prescales.push_back(3);  }
  if (Version==2) { trignames.push_back("MinBias_DoubleHF1"); prescales.push_back(2);  }
  if (Version==3) { trignames.push_back("MinBias_SingleHF2"); prescales.push_back(1);  }
  if (Version==3) { trignames.push_back("MinBias_DoubleHF2"); prescales.push_back(1);  }
  if (Version==4) { trignames.push_back("MinBias_SingleHF1"); prescales.push_back(10); }
  if (Version==4) { trignames.push_back("MinBias_DoubleHF1"); prescales.push_back(10); } 
  if (Version==5) { trignames.push_back("MinBias_SingleHF1"); prescales.push_back(3);  }
  if (Version==5) { trignames.push_back("MinBias_DoubleHF1"); prescales.push_back(5);  }
  if (Version==6) { trignames.push_back("MinBias_SingleHF2"); prescales.push_back(1);  }
  if (Version==6) { trignames.push_back("MinBias_DoubleHF2"); prescales.push_back(1);  }
  if (Version==7) { trignames.push_back("MinBias_SingleHF2"); prescales.push_back(1);  }
  if (Version==7) { trignames.push_back("MinBias_DoubleHF2"); prescales.push_back(1);  }
  if (Version==8) { trignames.push_back("MinBias_SingleHF2"); prescales.push_back(1);  }
  if (Version==8) { trignames.push_back("MinBias_DoubleHF2"); prescales.push_back(1);  }
  if (Version==9) { trignames.push_back("MinBias_SingleHF2"); prescales.push_back(1);  }
  if (Version==9) { trignames.push_back("MinBias_DoubleHF2"); prescales.push_back(1);  }
  if (Version==10) { trignames.push_back("MinBias_SingleHF2"); prescales.push_back(1);  }
  if (Version==10) { trignames.push_back("MinBias_DoubleHF2"); prescales.push_back(1);  }

  /**/

}


// Option 2
void MakeL1Menu_17E30_43(double &ILumi, double &nFilledBunches, vector<string> &trignames, vector<int> &prescales,
		       int Version = 1) {
  ILumi = 1.7E30;
  nFilledBunches = 43;

  /**/ // New Relaxed triggers for early run periods 
  trignames.push_back("SingleMu0"); prescales.push_back(1);
  trignames.push_back("SingleMu3"); prescales.push_back(1);	
  trignames.push_back("SingleMu5"); prescales.push_back(1);    
  //trignames.push_back("SingleMu7"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleMu3"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleMu5"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu7"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu14"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu25"); prescales.push_back(1);	

  trignames.push_back("L1_DoubleMu3"); prescales.push_back(1);	
  trignames.push_back("L1_TripleMu3"); prescales.push_back(1);	

  trignames.push_back("L1_Mu3_Jet15"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_Jet15"); prescales.push_back(1); 
  trignames.push_back("L1_Mu3_Jet70"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_Jet20"); prescales.push_back(1); 

  trignames.push_back("L1_Mu3_IsoEG5"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_IsoEG10"); prescales.push_back(1);
  trignames.push_back("L1_Mu3_EG12"); prescales.push_back(1);

  trignames.push_back("SingleEG5"); prescales.push_back(1);	
  //trignames.push_back("SingleEG6"); prescales.push_back(1);	
  //trignames.push_back("SingleEG7"); prescales.push_back(1);	
  trignames.push_back("SingleEG8"); prescales.push_back(1);	
  //trignames.push_back("SingleEG9"); prescales.push_back(1);	
  trignames.push_back("SingleEG10"); prescales.push_back(1);	
  trignames.push_back("SingleEG12"); prescales.push_back(1);	
  //trignames.push_back("SingleEG15"); prescales.push_back(1);
  //trignames.push_back("L1_SingleEG5"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG8"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG10"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG12"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG25"); prescales.push_back(1);	

  //trignames.push_back("SingleIsoEG12"); prescales.push_back(1);	
  //trignames.push_back("SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG8"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG12"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG25"); prescales.push_back(1);	

  //if (Version==1) { trignames.push_back("SingleJet10"); prescales.push_back(1);   }
  //if (Version==1) { trignames.push_back("SingleJet15"); prescales.push_back(1);	  }
  if (Version==1) { trignames.push_back("SingleJet20"); prescales.push_back(1);	  }
		  								  
  //if (Version==2) { trignames.push_back("SingleJet10"); prescales.push_back(10);  }
  if (Version==2) { trignames.push_back("SingleJet20"); prescales.push_back(1);	  }
  //if (Version==4) { trignames.push_back("SingleJet10"); prescales.push_back(10);  }
  if (Version==4) { trignames.push_back("SingleJet20"); prescales.push_back(1);	  }

  //if (Version==3) { trignames.push_back("SingleJet10"); prescales.push_back(100); }
  if (Version==3) { trignames.push_back("SingleJet20"); prescales.push_back(10);  }
  //if (Version==5) { trignames.push_back("SingleJet10"); prescales.push_back(100); }
  if (Version==5) { trignames.push_back("SingleJet20"); prescales.push_back(10);  }   

  //if (Version==6) { trignames.push_back("SingleJet10"); prescales.push_back(10);  }
  if (Version==6) { trignames.push_back("SingleJet20"); prescales.push_back(1);	  }

  //if (Version==7) { trignames.push_back("SingleJet10"); prescales.push_back(10);  }
  if (Version==7) { trignames.push_back("SingleJet20"); prescales.push_back(1);	  }

  if (Version==8) { trignames.push_back("SingleJet20"); prescales.push_back(10);	  }

  if (Version==9) { trignames.push_back("SingleJet15"); prescales.push_back(15);}
  //if (Version==9) { trignames.push_back("SingleJet20"); prescales.push_back(10);}
  
    //trignames.push_back("SingleJet25"); prescales.push_back(1);
  trignames.push_back("SingleJet30"); prescales.push_back(1);
  //trignames.push_back("SingleJet35"); prescales.push_back(1);
  //trignames.push_back("SingleJet40"); prescales.push_back(1);
  trignames.push_back("SingleJet50"); prescales.push_back(1);
  trignames.push_back("SingleJet70"); prescales.push_back(1);
  //trignames.push_back("SingleJet100"); prescales.push_back(1);
  //trignames.push_back("SingleJet150"); prescales.push_back(1);
  //trignames.push_back("L1_SingleJet30"); prescales.push_back(1);
  //trignames.push_back("L1_SingleJet70"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet100"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet150"); prescales.push_back(1);
  
  //if (Version==1) { trignames.push_back("SingleTauJet10"); prescales.push_back(1);	 }
  if (Version==1) { trignames.push_back("SingleTauJet20"); prescales.push_back(1);	 }
  if (Version==1) { trignames.push_back("SingleTauJet25"); prescales.push_back(1);	 }

  if (Version==2) { trignames.push_back("SingleTauJet20"); prescales.push_back(1);	 }
  if (Version==2) { trignames.push_back("SingleTauJet25"); prescales.push_back(1);	 }
  if (Version==4) { trignames.push_back("SingleTauJet20"); prescales.push_back(1);	 }
  if (Version==4) { trignames.push_back("SingleTauJet25"); prescales.push_back(1);	 }

  if (Version==6) { trignames.push_back("SingleTauJet20"); prescales.push_back(1);	 }
  if (Version==6) { trignames.push_back("SingleTauJet25"); prescales.push_back(1);	 }

  if (Version==7) { trignames.push_back("SingleTauJet20"); prescales.push_back(1);	 }
  if (Version==7) { trignames.push_back("SingleTauJet25"); prescales.push_back(1);	 }

  if (Version==8) { trignames.push_back("SingleTauJet20"); prescales.push_back(1);	 }
  if (Version==8) { trignames.push_back("SingleTauJet25"); prescales.push_back(1);	 }

  if (Version==9) { trignames.push_back("SingleTauJet20"); prescales.push_back(1);	 }
  if (Version==9) { trignames.push_back("SingleTauJet25"); prescales.push_back(1);	 }

  trignames.push_back("SingleTauJet30"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet35"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet40"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet50"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet60"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet70"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet80"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet100"); prescales.push_back(1);
  //trignames.push_back("L1_SingleTauJet10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet30"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet40"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet60"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleTauJet80"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet100"); prescales.push_back(1);
  
  trignames.push_back("L1_HTT200"); prescales.push_back(1);
  trignames.push_back("L1_HTT250"); prescales.push_back(1);
  trignames.push_back("L1_HTT300"); prescales.push_back(1);
  trignames.push_back("L1_HTT400"); prescales.push_back(1);
  //trignames.push_back("L1_HTT500"); prescales.push_back(1);

  if (Version<8) { trignames.push_back("ETM20"); prescales.push_back(1); }
  if (Version<8) { trignames.push_back("ETM25"); prescales.push_back(1); }

  trignames.push_back("ETM30"); prescales.push_back(1);	
  //trignames.push_back("ETM35"); prescales.push_back(1);	
  //trignames.push_back("ETM40"); prescales.push_back(1);	
  //trignames.push_back("ETM50"); prescales.push_back(1);	
  //trignames.push_back("ETM60"); prescales.push_back(1);	
  if (Version<8) { trignames.push_back("L1_ETM20"); prescales.push_back(1); }
  if (Version<8) { trignames.push_back("L1_ETM25"); prescales.push_back(1); }
  trignames.push_back("L1_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_ETM50"); prescales.push_back(1);
  trignames.push_back("L1_ETM60"); prescales.push_back(1);
  
  trignames.push_back("L1_DoubleIsoEG8"); prescales.push_back(1);
  trignames.push_back("L1_DoubleIsoEG10"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG5"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG10"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG15"); prescales.push_back(1);

  trignames.push_back("L1_DoubleJet70"); prescales.push_back(1);
  trignames.push_back("L1_DoubleJet100"); prescales.push_back(1);

  trignames.push_back("L1_DoubleTauJet20"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet30"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet35"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet40"); prescales.push_back(1);
  
  trignames.push_back("L1_IsoEG10_Jet15"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet30"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet70"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet20_ETM20"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_HTT100_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TripleJet50"); prescales.push_back(1);
  trignames.push_back("L1_QuadJet30"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveDoubleIsoEG6"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveDoubleJet60"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveJet25_Gap_Jet25"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet20_ForJet10"); prescales.push_back(1);
  trignames.push_back("L1_MinBias_HTT10"); prescales.push_back(1);
  trignames.push_back("L1_ZeroBias"); prescales.push_back(1);

  // New Minbias triggers
  if (Version==1) { trignames.push_back("MinBias_SingleHF1"); prescales.push_back(1);    }
  if (Version==1) { trignames.push_back("MinBias_DoubleHF1"); prescales.push_back(1);	 }
  if (Version==2) { trignames.push_back("MinBias_SingleHF1"); prescales.push_back(10);	 }
  if (Version==2) { trignames.push_back("MinBias_DoubleHF1"); prescales.push_back(10);	 }
  if (Version==3) { trignames.push_back("MinBias_SingleHF2"); prescales.push_back(1);	 }
  if (Version==3) { trignames.push_back("MinBias_DoubleHF2"); prescales.push_back(1);	 }
  if (Version==4) { trignames.push_back("MinBias_SingleHF1"); prescales.push_back(100);	 }
  if (Version==4) { trignames.push_back("MinBias_DoubleHF1"); prescales.push_back(100);	 }
  if (Version==5) { trignames.push_back("MinBias_SingleHF2"); prescales.push_back(10);	 }
  if (Version==5) { trignames.push_back("MinBias_DoubleHF2"); prescales.push_back(10);   } 
  if (Version==6) { trignames.push_back("MinBias_SingleHF2"); prescales.push_back(3);	 }
  if (Version==6) { trignames.push_back("MinBias_DoubleHF2"); prescales.push_back(3);	 }
  if (Version==7) { trignames.push_back("MinBias_SingleHF1"); prescales.push_back(10);	 }
  if (Version==8) { trignames.push_back("MinBias_SingleHF3"); prescales.push_back(1);	 }
  if (Version==8) { trignames.push_back("MinBias_DoubleHF3"); prescales.push_back(1);	 }

  if (Version==9) { trignames.push_back("MinBias_SingleHF3"); prescales.push_back(1);	 }
  if (Version==9) { trignames.push_back("MinBias_DoubleHF3"); prescales.push_back(1);	 }

  /**/

}

// Option 3
void MakeL1Menu_10E32_1000(double &ILumi, double &nFilledBunches, vector<string> &trignames, vector<int> &prescales,
		       int Version = 1) {
  ILumi = 1.0E32;
  //nFilledBunches = 1000;
  nFilledBunches = 936;

  /**/  // 160 table for lumi = 1E32
  //trignames.push_back("SingleMu0"); prescales.push_back(1);
  //trignames.push_back("SingleMu3"); prescales.push_back(4000);	
  //trignames.push_back("SingleMu5"); prescales.push_back(2000);    
  //trignames.push_back("SingleMu7"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu3"); prescales.push_back(4000);	
  trignames.push_back("L1_SingleMu5"); prescales.push_back(2000);	
  trignames.push_back("L1_SingleMu7"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu10"); prescales.push_back(1);
  trignames.push_back("L1_SingleMu14"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu25"); prescales.push_back(1);	

  trignames.push_back("L1_DoubleMu3"); prescales.push_back(1);	
  trignames.push_back("L1_TripleMu3"); prescales.push_back(1);
  
  trignames.push_back("L1_Mu3_Jet15"); prescales.push_back(20); 
  trignames.push_back("L1_Mu5_Jet15"); prescales.push_back(1); 
  trignames.push_back("L1_Mu3_Jet70"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_Jet20"); prescales.push_back(1); 

  trignames.push_back("L1_Mu3_IsoEG5"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_IsoEG10"); prescales.push_back(1);
  trignames.push_back("L1_Mu3_EG12"); prescales.push_back(1);

  //trignames.push_back("SingleEG5"); prescales.push_back(1);	
  //trignames.push_back("SingleEG6"); prescales.push_back(1);	
  //trignames.push_back("SingleEG7"); prescales.push_back(1);	
  //trignames.push_back("SingleEG8"); prescales.push_back(1);	
  //trignames.push_back("SingleEG9"); prescales.push_back(1);	
  //trignames.push_back("SingleEG10"); prescales.push_back(1);	
  //trignames.push_back("SingleEG12"); prescales.push_back(1);	
  //trignames.push_back("SingleEG15"); prescales.push_back(1);
  trignames.push_back("L1_SingleEG5"); prescales.push_back(10000);	
  trignames.push_back("L1_SingleEG8"); prescales.push_back(1000);	
  trignames.push_back("L1_SingleEG10"); prescales.push_back(100);	
  trignames.push_back("L1_SingleEG12"); prescales.push_back(100);	
  trignames.push_back("L1_SingleEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG25"); prescales.push_back(1);	

  //trignames.push_back("SingleIsoEG12"); prescales.push_back(1);	
  //trignames.push_back("SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG5"); prescales.push_back(10000);	
  trignames.push_back("L1_SingleIsoEG8"); prescales.push_back(1000);	
  trignames.push_back("L1_SingleIsoEG10"); prescales.push_back(100);	
  trignames.push_back("L1_SingleIsoEG12"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG25"); prescales.push_back(1);	

  //trignames.push_back("SingleJet10"); prescales.push_back(1);
  //trignames.push_back("SingleJet15"); prescales.push_back(1);
  //trignames.push_back("SingleJet20"); prescales.push_back(1);
  //trignames.push_back("SingleJet25"); prescales.push_back(1);
  //trignames.push_back("SingleJet30"); prescales.push_back(1);
  //trignames.push_back("SingleJet35"); prescales.push_back(1);
  //trignames.push_back("SingleJet40"); prescales.push_back(1);
  //trignames.push_back("SingleJet50"); prescales.push_back(1);
  //trignames.push_back("SingleJet70"); prescales.push_back(100);
  //trignames.push_back("SingleJet100"); prescales.push_back(1);
  //trignames.push_back("SingleJet150"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet15"); prescales.push_back(100000);
  //trignames.push_back("L1_SingleJet20"); prescales.push_back(999999999);
  trignames.push_back("L1_SingleJet30"); prescales.push_back(10000);
  //trignames.push_back("L1_SingleJet50"); prescales.push_back(999999999);
  trignames.push_back("L1_SingleJet70"); prescales.push_back(100);
  trignames.push_back("L1_SingleJet100"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet150"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet200"); prescales.push_back(1);
  
  //trignames.push_back("SingleTauJet10"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet20"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet25"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet30"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet35"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet40"); prescales.push_back(1000);	
  //trignames.push_back("SingleTauJet50"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet60"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet70"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet80"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet100"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet10"); prescales.push_back(100000);	
  trignames.push_back("L1_SingleTauJet20"); prescales.push_back(100000);	
  trignames.push_back("L1_SingleTauJet30"); prescales.push_back(10000);	
  trignames.push_back("L1_SingleTauJet40"); prescales.push_back(1000);	
  trignames.push_back("L1_SingleTauJet60"); prescales.push_back(1000);	
  trignames.push_back("L1_SingleTauJet80"); prescales.push_back(1);	

  
  trignames.push_back("L1_HTT100"); prescales.push_back(10000);
  trignames.push_back("L1_HTT200"); prescales.push_back(1000);
  
  if (Version<5) {trignames.push_back("L1_HTT250"); prescales.push_back(1);}
  if (Version<5) {trignames.push_back("L1_HTT300"); prescales.push_back(1);}
  if (Version<5) {trignames.push_back("L1_HTT400"); prescales.push_back(1);}

  if (Version==5) {trignames.push_back("L1_HTT250"); prescales.push_back(100);}
  if (Version==5) {trignames.push_back("L1_HTT300"); prescales.push_back(10);}
  if (Version==5) {trignames.push_back("L1_HTT400"); prescales.push_back(1);}

  if (Version==6) {trignames.push_back("L1_HTT250"); prescales.push_back(100);}
  if (Version==6) {trignames.push_back("L1_HTT300"); prescales.push_back(1);}
  if (Version==6) {trignames.push_back("L1_HTT400"); prescales.push_back(1);}

  if (Version==7) {trignames.push_back("L1_HTT250"); prescales.push_back(100);}
  if (Version==7) {trignames.push_back("L1_HTT300"); prescales.push_back(1);}
  if (Version==7) {trignames.push_back("L1_HTT400"); prescales.push_back(1);}

  if (Version==8) {trignames.push_back("L1_HTT250"); prescales.push_back(100);}
  if (Version==8) {trignames.push_back("L1_HTT300"); prescales.push_back(1);}
  if (Version==8) {trignames.push_back("L1_HTT400"); prescales.push_back(1);}

  if (Version==9) {trignames.push_back("L1_HTT250"); prescales.push_back(100);}
  if (Version==9) {trignames.push_back("L1_HTT300"); prescales.push_back(1);}
  if (Version==9) {trignames.push_back("L1_HTT400"); prescales.push_back(1);}
  if (Version==9) {trignames.push_back("L1_HTT500"); prescales.push_back(1);}


  if (Version<8) {trignames.push_back("ETM20"); prescales.push_back(1000000);}

  //trignames.push_back("ETM25"); prescales.push_back(1);	
  //trignames.push_back("ETM35"); prescales.push_back(1);	
  if (Version==1) {trignames.push_back("ETM30"); prescales.push_back(1);}
  if (Version==2) {trignames.push_back("ETM30"); prescales.push_back(100000);}
  if (Version==2) {trignames.push_back("ETM40"); prescales.push_back(1);}
  if (Version==3) {trignames.push_back("ETM30"); prescales.push_back(100000);}
  if (Version==3) {trignames.push_back("ETM40"); prescales.push_back(100000);}

  if (Version==4) {trignames.push_back("ETM30"); prescales.push_back(100000);}
  if (Version==4) {trignames.push_back("ETM40"); prescales.push_back(100000);}
  if (Version==4) {trignames.push_back("ETM40_Jet30"); prescales.push_back(1);}
  if (Version==4) {trignames.push_back("ETM50_Jet40"); prescales.push_back(1);}

  if (Version==5) {trignames.push_back("ETM30"); prescales.push_back(100000);}
  if (Version==5) {trignames.push_back("ETM40"); prescales.push_back(1);}
  if (Version==5) {trignames.push_back("ETM40_Jet30"); prescales.push_back(1);}
  if (Version==5) {trignames.push_back("ETM50_Jet40"); prescales.push_back(1);}

  if (Version==6) {trignames.push_back("ETM30"); prescales.push_back(100000);}
  if (Version==6) {trignames.push_back("ETM40"); prescales.push_back(100000);}
  if (Version==6) {trignames.push_back("ETM45"); prescales.push_back(1);}
  if (Version==6) {trignames.push_back("ETM45_Jet30"); prescales.push_back(1);}

  if (Version==7) {trignames.push_back("ETM30"); prescales.push_back(100000);}
  if (Version==7) {trignames.push_back("ETM40"); prescales.push_back(1);}
  if (Version==7) {trignames.push_back("ETM45"); prescales.push_back(1);}
  if (Version==7) {trignames.push_back("ETM45_Jet30"); prescales.push_back(1);}

  if (Version==8) {trignames.push_back("ETM10"); prescales.push_back(10000);}
  if (Version==8) {trignames.push_back("ETM15"); prescales.push_back(5000);}
  if (Version==8) {trignames.push_back("ETM20"); prescales.push_back(1000);}
  if (Version==8) {trignames.push_back("ETM25"); prescales.push_back(1000);}
  if (Version==8) {trignames.push_back("ETM30"); prescales.push_back(1000);}
  if (Version==8) {trignames.push_back("ETM40"); prescales.push_back(1000);}
  if (Version==8) {trignames.push_back("ETM45"); prescales.push_back(1);}
  if (Version==8) {trignames.push_back("ETM45_Jet30"); prescales.push_back(1);}

  if (Version==8) {trignames.push_back("L1_ETM20"); prescales.push_back(1000);}
  if (Version==8) {trignames.push_back("L1_ETM30"); prescales.push_back(1000);}
  if (Version==8) {trignames.push_back("L1_ETM40"); prescales.push_back(1000);}

  if (Version==8) {trignames.push_back("L1_ETM20"); prescales.push_back(10000);}
  if (Version==8) {trignames.push_back("L1_ETM30"); prescales.push_back(1000);}
  if (Version==8) {trignames.push_back("L1_ETM40"); prescales.push_back(1000);}

  trignames.push_back("L1_ETM50"); prescales.push_back(1);	
  trignames.push_back("L1_ETM60"); prescales.push_back(1);	
  //trignames.push_back("ETM60"); prescales.push_back(1);	
  //trignames.push_back("ETM50_Jet40"); prescales.push_back(1);

  if (Version<8) {trignames.push_back("L1_ETM20"); prescales.push_back(1000000);}
  if (Version==1) {trignames.push_back("L1_ETM30"); prescales.push_back(1);}
  if (Version==2) {trignames.push_back("L1_ETM30"); prescales.push_back(100000);}
  if (Version==2) {trignames.push_back("L1_ETM40"); prescales.push_back(1);}
  if (Version==3) {trignames.push_back("L1_ETM30"); prescales.push_back(100000);}
  if (Version==3) {trignames.push_back("L1_ETM40"); prescales.push_back(100000);}
  if (Version==4) {trignames.push_back("L1_ETM30"); prescales.push_back(100000);}
  if (Version==4) {trignames.push_back("L1_ETM40"); prescales.push_back(100000);}
  if (Version==5) {trignames.push_back("L1_ETM30"); prescales.push_back(100000);}
  if (Version==5) {trignames.push_back("L1_ETM40"); prescales.push_back(1);}
  if (Version==6) {trignames.push_back("L1_ETM30"); prescales.push_back(100000);}
  if (Version==6) {trignames.push_back("L1_ETM40"); prescales.push_back(100000);}
  if (Version==7) {trignames.push_back("L1_ETM30"); prescales.push_back(100000);}
  if (Version==7) {trignames.push_back("L1_ETM40"); prescales.push_back(1);}
  
  //trignames.push_back("L1_ETM50"); prescales.push_back(1);
  //trignames.push_back("L1_ETM60"); prescales.push_back(1);
  trignames.push_back("L1_DoubleIsoEG8"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG10"); prescales.push_back(1);
  trignames.push_back("L1_DoubleJet70"); prescales.push_back(1);
  trignames.push_back("L1_DoubleJet100"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet40"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet15"); prescales.push_back(20);
  trignames.push_back("L1_IsoEG10_Jet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet30"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet70"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet30"); prescales.push_back(1);
  //trignames.push_back("L1_TauJet20_ETM20"); prescales.push_back(100000);
  trignames.push_back("L1_TauJet30_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_HTT100_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TripleJet50"); prescales.push_back(1);
  if (Version==1) {trignames.push_back("L1_QuadJet30"); prescales.push_back(1);}
  if (Version==2) {trignames.push_back("QuadJet40"); prescales.push_back(1);}
  if (Version==3) {trignames.push_back("QuadJet40"); prescales.push_back(1);}
  if (Version==3) {trignames.push_back("QuadJet50"); prescales.push_back(1);}
  if (Version==4) {trignames.push_back("QuadJet40"); prescales.push_back(1);}
  if (Version==4) {trignames.push_back("QuadJet50"); prescales.push_back(1);}

  if (Version==5) {trignames.push_back("QuadJet40"); prescales.push_back(1);}
  if (Version==5) {trignames.push_back("QuadJet50"); prescales.push_back(1);}

  if (Version==6) {trignames.push_back("QuadJet40"); prescales.push_back(1);}
  if (Version==6) {trignames.push_back("QuadJet50"); prescales.push_back(1);}

  if (Version==7) {trignames.push_back("QuadJet40"); prescales.push_back(1);}
  if (Version==7) {trignames.push_back("QuadJet50"); prescales.push_back(1);}

  if (Version==8) {trignames.push_back("QuadJet40"); prescales.push_back(1);}
  if (Version==8) {trignames.push_back("QuadJet50"); prescales.push_back(1);}

  if (Version==9) {trignames.push_back("QuadJet40"); prescales.push_back(1);}
  if (Version==9) {trignames.push_back("QuadJet50"); prescales.push_back(1);}

  trignames.push_back("L1_ExclusiveDoubleIsoEG6"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveDoubleJet60"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveJet25_Gap_Jet25"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet20_ForJet10"); prescales.push_back(1);
  trignames.push_back("L1_MinBias_HTT10"); prescales.push_back(1);
  trignames.push_back("L1_ZeroBias"); prescales.push_back(1);
  /**/
}


////



// Option 4
void MakeL1Menu_61E30_43(double &ILumi, double &nFilledBunches, 
			 vector<string> &trignames, vector<int> &prescales, int Version = 1) {
  ILumi = 6.1E30;
  nFilledBunches = 43;

  /**/ // New Relaxed triggers for early run periods 
  trignames.push_back("SingleMu0"); prescales.push_back(1);
  trignames.push_back("SingleMu3"); prescales.push_back(1);	
  trignames.push_back("SingleMu5"); prescales.push_back(1);    
  //trignames.push_back("SingleMu7"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleMu3"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleMu5"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu7"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu14"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu25"); prescales.push_back(1);	

  trignames.push_back("L1_DoubleMu3"); prescales.push_back(1);	
  trignames.push_back("L1_TripleMu3"); prescales.push_back(1);
  
  trignames.push_back("L1_Mu3_Jet15"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_Jet15"); prescales.push_back(1); 
  trignames.push_back("L1_Mu3_Jet70"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_Jet20"); prescales.push_back(1);
  
  trignames.push_back("L1_Mu3_IsoEG5"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_IsoEG10"); prescales.push_back(1);
  trignames.push_back("L1_Mu3_EG12"); prescales.push_back(1);
  
  if (Version==1) { trignames.push_back("SingleEG5"); prescales.push_back(5); };	
  if (Version==2) { trignames.push_back("SingleEG5"); prescales.push_back(10); };	
  if (Version==2) { trignames.push_back("SingleEG6"); prescales.push_back(1); };
  if (Version==3) { trignames.push_back("SingleEG5"); prescales.push_back(10); };	
  if (Version==3) { trignames.push_back("SingleEG6"); prescales.push_back(1); };
  
  //trignames.push_back("SingleEG7"); prescales.push_back(1);	
  trignames.push_back("SingleEG8"); prescales.push_back(1);	
  //trignames.push_back("SingleEG9"); prescales.push_back(1);	
  trignames.push_back("SingleEG10"); prescales.push_back(1);	
  trignames.push_back("SingleEG12"); prescales.push_back(1);	
  //trignames.push_back("SingleEG15"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG5"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG8"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG10"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG12"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG25"); prescales.push_back(1);	

  //trignames.push_back("SingleIsoEG12"); prescales.push_back(1);	
  //trignames.push_back("SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG8"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG12"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG25"); prescales.push_back(1);	

  //if (Version==1) { trignames.push_back("SingleJet10"); prescales.push_back(10);  }
  if (Version==1) { trignames.push_back("SingleJet20"); prescales.push_back(100);	  }
  if (Version==1) { trignames.push_back("SingleJet30"); prescales.push_back(10); }

  if (Version==2) { trignames.push_back("SingleJet20"); prescales.push_back(100);	  }
  if (Version==2) { trignames.push_back("SingleJet30"); prescales.push_back(10); }

  if (Version==3) { trignames.push_back("SingleJet15"); prescales.push_back(100);	  }
  if (Version==3) { trignames.push_back("SingleJet20"); prescales.push_back(100);	  }
  if (Version==3) { trignames.push_back("SingleJet30"); prescales.push_back(10); }

  trignames.push_back("SingleJet40"); prescales.push_back(1);
  trignames.push_back("SingleJet50"); prescales.push_back(1);
  trignames.push_back("SingleJet70"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet100"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet150"); prescales.push_back(1);

  //if (Version==1) { trignames.push_back("SingleTauJet20"); prescales.push_back(10);	 }
  //if (Version==1) { trignames.push_back("SingleTauJet25"); prescales.push_back(1);	 }

  //trignames.push_back("SingleTauJet10"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet20"); prescales.push_back(100);	
  //trignames.push_back("SingleTauJet25"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet30"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet35"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet40"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet50"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet60"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet70"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet80"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet100"); prescales.push_back(1);
  //trignames.push_back("L1_SingleTauJet10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet30"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet40"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet60"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleTauJet80"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet100"); prescales.push_back(1);

  trignames.push_back("L1_HTT200"); prescales.push_back(1);
  trignames.push_back("L1_HTT250"); prescales.push_back(1);
  trignames.push_back("L1_HTT300"); prescales.push_back(1);
  trignames.push_back("L1_HTT400"); prescales.push_back(1);
  //trignames.push_back("L1_HTT500"); prescales.push_back(1);

  if (Version==1) { trignames.push_back("ETM20"); prescales.push_back(1); };
  if (Version==2) { trignames.push_back("ETM20"); prescales.push_back(100000); };

  if (Version==1) { trignames.push_back("ETM25"); prescales.push_back(1);}
  if (Version==2) { trignames.push_back("ETM25"); prescales.push_back(100000);}

  //trignames.push_back("ETM20"); prescales.push_back(1);
  //trignames.push_back("ETM25"); prescales.push_back(1);	
  //trignames.push_back("ETM30"); prescales.push_back(1);	
  //trignames.push_back("ETM35"); prescales.push_back(1);	
  //trignames.push_back("ETM40"); prescales.push_back(1);	
  //trignames.push_back("ETM50"); prescales.push_back(1);	
  //trignames.push_back("ETM60"); prescales.push_back(1);	
  if (Version==1) { trignames.push_back("L1_ETM20"); prescales.push_back(1); };
  if (Version==2) { trignames.push_back("L1_ETM20"); prescales.push_back(100000); };
  trignames.push_back("L1_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_ETM50"); prescales.push_back(1);
  trignames.push_back("L1_ETM60"); prescales.push_back(1);  
  
  trignames.push_back("L1_DoubleIsoEG8"); prescales.push_back(1);
  trignames.push_back("L1_DoubleIsoEG10"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG5"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG10"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG15"); prescales.push_back(1);

  trignames.push_back("L1_DoubleJet70"); prescales.push_back(1);
  trignames.push_back("L1_DoubleJet100"); prescales.push_back(1);

  trignames.push_back("L1_DoubleTauJet20"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet30"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet35"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet40"); prescales.push_back(1);
  
  trignames.push_back("L1_IsoEG10_Jet15"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet30"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet70"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet20_ETM20"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_HTT100_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TripleJet50"); prescales.push_back(1);
  trignames.push_back("L1_QuadJet30"); prescales.push_back(1);
  //if (Version==2) {trignames.push_back("QuadJet40"); prescales.push_back(1);}
  trignames.push_back("L1_ExclusiveDoubleIsoEG6"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveDoubleJet60"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveJet25_Gap_Jet25"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet20_ForJet10"); prescales.push_back(1);
  trignames.push_back("L1_MinBias_HTT10"); prescales.push_back(1);
  trignames.push_back("L1_ZeroBias"); prescales.push_back(1);

  // New Minbias triggers
  //if (Version==1) { trignames.push_back("MinBias_SingleHF2"); prescales.push_back(50);	 }
  if (Version==1) { trignames.push_back("MinBias_DoubleHF2"); prescales.push_back(10);	 }
  if (Version==2) { trignames.push_back("MinBias_SingleHF3"); prescales.push_back(30);	 }
  if (Version==2) { trignames.push_back("MinBias_DoubleHF3"); prescales.push_back(10);	 }
  if (Version==3) { trignames.push_back("MinBias_SingleHF3"); prescales.push_back(20);	 }
  if (Version==3) { trignames.push_back("MinBias_DoubleHF3"); prescales.push_back(10);	 }

  /**/

}

// Option 5
void MakeL1Menu_11E31_156(double &ILumi, double &nFilledBunches, 
			 vector<string> &trignames, vector<int> &prescales, int Version = 1) {
  ILumi = 1.1E31;
  nFilledBunches = 156;

  /**/ // New Relaxed triggers for early run periods 
  trignames.push_back("SingleMu0"); prescales.push_back(1);
  trignames.push_back("SingleMu3"); prescales.push_back(1);	
  trignames.push_back("SingleMu5"); prescales.push_back(1);    
  //trignames.push_back("SingleMu7"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleMu3"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleMu5"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu7"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu14"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu25"); prescales.push_back(1);	

  trignames.push_back("L1_DoubleMu3"); prescales.push_back(1);	
  trignames.push_back("L1_TripleMu3"); prescales.push_back(1);
  
  trignames.push_back("L1_Mu3_Jet15"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_Jet15"); prescales.push_back(1); 
  trignames.push_back("L1_Mu3_Jet70"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_Jet20"); prescales.push_back(1);
  
  trignames.push_back("L1_Mu3_IsoEG5"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_IsoEG10"); prescales.push_back(1);
  trignames.push_back("L1_Mu3_EG12"); prescales.push_back(1);
  

  if (Version==1) { trignames.push_back("SingleEG5"); prescales.push_back(1); };	

  if (Version==2) { trignames.push_back("SingleEG5"); prescales.push_back(100); };	
  if (Version==2) { trignames.push_back("SingleEG6"); prescales.push_back(100); };	
  if (Version==2) { trignames.push_back("SingleEG7"); prescales.push_back(1); };	

  if (Version==3) { trignames.push_back("SingleEG5"); prescales.push_back(100); };	
  if (Version==3) { trignames.push_back("SingleEG6"); prescales.push_back(100); };	
  if (Version==3) { trignames.push_back("SingleEG7"); prescales.push_back(1); };	

  trignames.push_back("SingleEG8"); prescales.push_back(1);	
  //trignames.push_back("SingleEG9"); prescales.push_back(1);	
  trignames.push_back("SingleEG10"); prescales.push_back(1);	
  trignames.push_back("SingleEG12"); prescales.push_back(1);	
  //trignames.push_back("SingleEG15"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG5"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG8"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG10"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleEG12"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG25"); prescales.push_back(1);	

    //trignames.push_back("SingleIsoEG12"); prescales.push_back(1);	
  //trignames.push_back("SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG8"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG12"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG25"); prescales.push_back(1);	

  //if (Version==1) { trignames.push_back("SingleJet10"); prescales.push_back(10);  }
  if (Version==1) { trignames.push_back("SingleJet20"); prescales.push_back(1000);	  }
  if (Version==1) { trignames.push_back("SingleJet30"); prescales.push_back(100); }
  if (Version==1) { trignames.push_back("SingleJet40"); prescales.push_back(1); }

  if (Version==2) { trignames.push_back("SingleJet20"); prescales.push_back(1000);	  }
  if (Version==2) { trignames.push_back("SingleJet30"); prescales.push_back(100); }
  if (Version==2) { trignames.push_back("SingleJet40"); prescales.push_back(10); }

  if (Version==3) { trignames.push_back("SingleJet15"); prescales.push_back(1000);	  }
  if (Version==3) { trignames.push_back("SingleJet20"); prescales.push_back(1000);	  }
  if (Version==3) { trignames.push_back("SingleJet30"); prescales.push_back(100); }
  if (Version==3) { trignames.push_back("SingleJet40"); prescales.push_back(10); }

  //trignames.push_back("SingleJet35"); prescales.push_back(1);
  //trignames.push_back("SingleJet40"); prescales.push_back(1);
  trignames.push_back("SingleJet50"); prescales.push_back(1);
  trignames.push_back("SingleJet70"); prescales.push_back(1);
  //trignames.push_back("SingleJet100"); prescales.push_back(1);
  //trignames.push_back("SingleJet150"); prescales.push_back(1);
  //trignames.push_back("L1_SingleJet30"); prescales.push_back(1);
  //trignames.push_back("L1_SingleJet70"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet100"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet150"); prescales.push_back(1);


  //trignames.push_back("SingleTauJet10"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet20"); prescales.push_back(1000);	
  //trignames.push_back("SingleTauJet25"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet30"); prescales.push_back(100);	
  //trignames.push_back("SingleTauJet35"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet40"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet50"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet60"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet70"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet80"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet100"); prescales.push_back(1);
  //trignames.push_back("L1_SingleTauJet10"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet30"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet40"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet60"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleTauJet80"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet100"); prescales.push_back(1);

  trignames.push_back("L1_HTT200"); prescales.push_back(1);
  trignames.push_back("L1_HTT250"); prescales.push_back(1);
  trignames.push_back("L1_HTT300"); prescales.push_back(1);
  trignames.push_back("L1_HTT400"); prescales.push_back(1);
  //trignames.push_back("L1_HTT500"); prescales.push_back(1);

  if (Version==1) {trignames.push_back("ETM20"); prescales.push_back(100);}
  if (Version==2) {trignames.push_back("ETM20"); prescales.push_back(10000);}
  if (Version==3) {trignames.push_back("ETM20"); prescales.push_back(1000);}

  trignames.push_back("L1_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_ETM50"); prescales.push_back(1);
  trignames.push_back("L1_ETM60"); prescales.push_back(1);
  
  trignames.push_back("L1_DoubleIsoEG8"); prescales.push_back(1);
  trignames.push_back("L1_DoubleIsoEG10"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG5"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG10"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG15"); prescales.push_back(1);

  trignames.push_back("L1_DoubleJet70"); prescales.push_back(1);
  trignames.push_back("L1_DoubleJet100"); prescales.push_back(1);

  trignames.push_back("L1_DoubleTauJet20"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet30"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet35"); prescales.push_back(1);
  trignames.push_back("L1_DoubleTauJet40"); prescales.push_back(1);
  
  trignames.push_back("L1_IsoEG10_Jet15"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet30"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet70"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet20_ETM20"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_HTT100_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TripleJet50"); prescales.push_back(1);
  trignames.push_back("L1_QuadJet30"); prescales.push_back(1);
  //if (Version==2) {trignames.push_back("QuadJet40"); prescales.push_back(1);}
  trignames.push_back("L1_ExclusiveDoubleIsoEG6"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveDoubleJet60"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveJet25_Gap_Jet25"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet20_ForJet10"); prescales.push_back(1);
  trignames.push_back("L1_MinBias_HTT10"); prescales.push_back(1);
  trignames.push_back("L1_ZeroBias"); prescales.push_back(1);

  // New Minbias triggers
  //if (Version==1) { trignames.push_back("MinBias_SingleHF2"); prescales.push_back(50);	 }
  if (Version==1) { trignames.push_back("MinBias_DoubleHF2"); prescales.push_back(30);	 }

  if (Version==2) { trignames.push_back("MinBias_SingleHF3"); prescales.push_back(50);	 }
  if (Version==2) { trignames.push_back("MinBias_DoubleHF3"); prescales.push_back(50);	 }

  if (Version==3) { trignames.push_back("MinBias_SingleHF3"); prescales.push_back(50);	 }
  if (Version==3) { trignames.push_back("MinBias_DoubleHF3"); prescales.push_back(50);	 }

  /**/

}

// Option 6
void MakeL1Menu_56E31_156(double &ILumi, double &nFilledBunches, 
			 vector<string> &trignames, vector<int> &prescales, int Version = 1) {
  ILumi = 5.6E31;
  nFilledBunches = 156;

  /**/  // 160 table for lumi = 1E32

  //trignames.push_back("SingleMu0"); prescales.push_back(1);
  trignames.push_back("SingleMu3"); prescales.push_back(4000);	
  trignames.push_back("SingleMu5"); prescales.push_back(2000);    
  //trignames.push_back("SingleMu7"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleMu3"); prescales.push_back(4000);	
  //trignames.push_back("L1_SingleMu5"); prescales.push_back(2000);	
  trignames.push_back("L1_SingleMu7"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu10"); prescales.push_back(1);
  trignames.push_back("L1_SingleMu14"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleMu25"); prescales.push_back(1);	
  
  trignames.push_back("L1_DoubleMu3"); prescales.push_back(1);	
  trignames.push_back("L1_TripleMu3"); prescales.push_back(1);
  
  trignames.push_back("L1_Mu3_Jet15"); prescales.push_back(20); 
  trignames.push_back("L1_Mu5_Jet15"); prescales.push_back(1); 
  trignames.push_back("L1_Mu3_Jet70"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_Jet20"); prescales.push_back(1);
  
  trignames.push_back("L1_Mu3_IsoEG5"); prescales.push_back(1); 
  trignames.push_back("L1_Mu5_IsoEG10"); prescales.push_back(1);
  trignames.push_back("L1_Mu3_EG12"); prescales.push_back(1);
  
  //trignames.push_back("SingleEG5"); prescales.push_back(1);	
  //trignames.push_back("SingleEG6"); prescales.push_back(1);	
  //trignames.push_back("SingleEG7"); prescales.push_back(1);	
  //trignames.push_back("SingleEG8"); prescales.push_back(1);	
  //trignames.push_back("SingleEG9"); prescales.push_back(1);	
  //trignames.push_back("SingleEG10"); prescales.push_back(1);	
  //trignames.push_back("SingleEG12"); prescales.push_back(1);	
  //trignames.push_back("SingleEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG5"); prescales.push_back(10000);	
  trignames.push_back("L1_SingleEG8"); prescales.push_back(1000);	
  trignames.push_back("L1_SingleEG10"); prescales.push_back(100);	
  trignames.push_back("L1_SingleEG12"); prescales.push_back(100);	
  trignames.push_back("L1_SingleEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleEG25"); prescales.push_back(1);	

  //trignames.push_back("SingleIsoEG12"); prescales.push_back(1);	
  //trignames.push_back("SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG5"); prescales.push_back(10000);	
  trignames.push_back("L1_SingleIsoEG8"); prescales.push_back(1000);	
  trignames.push_back("L1_SingleIsoEG10"); prescales.push_back(10);	
  trignames.push_back("L1_SingleIsoEG12"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG15"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG20"); prescales.push_back(1);	
  trignames.push_back("L1_SingleIsoEG25"); prescales.push_back(1);

  if (Version==3) {
    //trignames.push_back("SingleJet10"); prescales.push_back(1);
    trignames.push_back("SingleJet15"); prescales.push_back(100000);
    trignames.push_back("SingleJet20"); prescales.push_back(100000);
    //trignames.push_back("SingleJet25"); prescales.push_back(1);
    trignames.push_back("SingleJet30"); prescales.push_back(10000);
    //trignames.push_back("SingleJet35"); prescales.push_back(1);
    //trignames.push_back("SingleJet40"); prescales.push_back(1);
    trignames.push_back("SingleJet50"); prescales.push_back(10000);
    trignames.push_back("SingleJet70"); prescales.push_back(100);
  }
  //trignames.push_back("SingleJet100"); prescales.push_back(1);
  //trignames.push_back("SingleJet150"); prescales.push_back(1);
  //trignames.push_back("L1_SingleJet15"); prescales.push_back(100000);
  //trignames.push_back("L1_SingleJet20"); prescales.push_back(999999999);
  //trignames.push_back("L1_SingleJet30"); prescales.push_back(10000);
  //trignames.push_back("L1_SingleJet50"); prescales.push_back(999999999);
  //trignames.push_back("L1_SingleJet70"); prescales.push_back(100);
  //trignames.push_back("L1_SingleJet100"); prescales.push_back(1);
  //trignames.push_back("L1_SingleJet150"); prescales.push_back(1);
  //trignames.push_back("L1_SingleJet200"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet100"); prescales.push_back(1);
  trignames.push_back("L1_SingleJet150"); prescales.push_back(1);

  trignames.push_back("SingleTauJet10"); prescales.push_back(100000);	
  trignames.push_back("SingleTauJet20"); prescales.push_back(10000);	
  //trignames.push_back("SingleTauJet25"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet30"); prescales.push_back(1000);	
  //trignames.push_back("SingleTauJet35"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet40"); prescales.push_back(1000);	
  //trignames.push_back("SingleTauJet50"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet60"); prescales.push_back(1000);	
  //trignames.push_back("SingleTauJet70"); prescales.push_back(1);	
  trignames.push_back("SingleTauJet80"); prescales.push_back(1);	
  //trignames.push_back("SingleTauJet100"); prescales.push_back(1);
  //trignames.push_back("L1_SingleTauJet10"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleTauJet20"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleTauJet30"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleTauJet40"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleTauJet60"); prescales.push_back(1);	
  //trignames.push_back("L1_SingleTauJet80"); prescales.push_back(1);	
  trignames.push_back("L1_SingleTauJet100"); prescales.push_back(1);
  
  trignames.push_back("L1_HTT100"); prescales.push_back(10000);
  trignames.push_back("L1_HTT200"); prescales.push_back(1000);
  
  if (Version<2) {trignames.push_back("L1_HTT250"); prescales.push_back(1);}
  if (Version<2) {trignames.push_back("L1_HTT300"); prescales.push_back(1);}
  if (Version<2) {trignames.push_back("L1_HTT400"); prescales.push_back(1);}

  if (Version==2) {trignames.push_back("L1_HTT250"); prescales.push_back(100);}
  if (Version==2) {trignames.push_back("L1_HTT300"); prescales.push_back(1);}
  if (Version==2) {trignames.push_back("L1_HTT400"); prescales.push_back(1);}

  if (Version==3) {trignames.push_back("L1_HTT250"); prescales.push_back(100);}
  if (Version==3) {trignames.push_back("L1_HTT300"); prescales.push_back(1);}
  if (Version==3) {trignames.push_back("L1_HTT400"); prescales.push_back(1);}

  //trignames.push_back("ETM20"); prescales.push_back(10000);
  //trignames.push_back("ETM25"); prescales.push_back(1);	
  //trignames.push_back("ETM35"); prescales.push_back(1);	
  if (Version==1) {trignames.push_back("ETM30"); prescales.push_back(1);}

  if (Version==2) {trignames.push_back("ETM30"); prescales.push_back(100000);}
  if (Version==2) {trignames.push_back("ETM40"); prescales.push_back(1);}
  if (Version==2) {trignames.push_back("ETM45"); prescales.push_back(1);}
  if (Version==2) {trignames.push_back("ETM45_Jet30"); prescales.push_back(1);}

  if (Version==3) {trignames.push_back("ETM20"); prescales.push_back(100000);}
  if (Version==3) {trignames.push_back("ETM30"); prescales.push_back(10000);}
  if (Version==3) {trignames.push_back("ETM40"); prescales.push_back(1);}
  if (Version==3) {trignames.push_back("ETM45"); prescales.push_back(1);}
  if (Version==3) {trignames.push_back("ETM45_Jet30"); prescales.push_back(1);}

  //trignames.push_back("L1_ETM30"); prescales.push_back(1);
  //trignames.push_back("L1_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_ETM50"); prescales.push_back(1);
  trignames.push_back("L1_ETM60"); prescales.push_back(1);
  //trignames.push_back("ETM50_Jet40"); prescales.push_back(1);

  //trignames.push_back("L1_ETM20"); prescales.push_back(100000);
  if (Version==1) {trignames.push_back("L1_ETM30"); prescales.push_back(1);}
  if (Version==2) {trignames.push_back("L1_ETM30"); prescales.push_back(100000);}
  if (Version==2) {trignames.push_back("L1_ETM40"); prescales.push_back(1);}
  trignames.push_back("L1_ETM50"); prescales.push_back(1);
  
  trignames.push_back("L1_DoubleIsoEG8"); prescales.push_back(1);
  trignames.push_back("L1_DoubleIsoEG10"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG5"); prescales.push_back(1000);
  trignames.push_back("L1_DoubleEG10"); prescales.push_back(1);
  trignames.push_back("L1_DoubleEG15"); prescales.push_back(1);

  trignames.push_back("L1_DoubleJet70"); prescales.push_back(1);
  trignames.push_back("L1_DoubleJet100"); prescales.push_back(1);

  trignames.push_back("L1_DoubleTauJet20"); prescales.push_back(1000);
  trignames.push_back("L1_DoubleTauJet30"); prescales.push_back(100);
  trignames.push_back("L1_DoubleTauJet35"); prescales.push_back(100);
  trignames.push_back("L1_DoubleTauJet40"); prescales.push_back(1);
  
  trignames.push_back("L1_IsoEG10_Jet15"); prescales.push_back(20);
  trignames.push_back("L1_IsoEG10_Jet20"); prescales.push_back(10);
  trignames.push_back("L1_IsoEG10_Jet30"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet70"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet20"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_TauJet30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet20_ETM20"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TauJet30_ETM40"); prescales.push_back(1);
  trignames.push_back("L1_HTT100_ETM30"); prescales.push_back(1);
  trignames.push_back("L1_TripleJet50"); prescales.push_back(1);
  trignames.push_back("L1_QuadJet30"); prescales.push_back(1);
  //if (Version==2) {trignames.push_back("QuadJet40"); prescales.push_back(1);}
  trignames.push_back("L1_ExclusiveDoubleIsoEG6"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveDoubleJet60"); prescales.push_back(1);
  trignames.push_back("L1_ExclusiveJet25_Gap_Jet25"); prescales.push_back(1);
  trignames.push_back("L1_IsoEG10_Jet20_ForJet10"); prescales.push_back(1);
  trignames.push_back("L1_MinBias_HTT10"); prescales.push_back(1);
  trignames.push_back("L1_ZeroBias"); prescales.push_back(1);

  if (Version==2) { trignames.push_back("MinBias_SingleHF3"); prescales.push_back(300);	 }
  if (Version==2) { trignames.push_back("MinBias_DoubleHF3"); prescales.push_back(300);	 }

  if (Version==3) { trignames.push_back("MinBias_SingleHF3"); prescales.push_back(300);	 }
  if (Version==3) { trignames.push_back("MinBias_DoubleHF3"); prescales.push_back(300);	 }
  /**/

}




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

