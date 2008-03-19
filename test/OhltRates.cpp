/////////////////////////////////////////////////////////////////////////////////////////////////
//
//          Macro calculating overlaps between different triggers, Hltrigger individual- and
//          pure-rates, taking into account background and signal type samples
//
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include "OHltTree.h"

#include "TH1.h"
#include "TChain.h"
#include "TCut.h"

#include <map>

using namespace std;

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


int main(int argc, char *argv[]){
  
  int NEntries = -1;
  if (argc>1) {
    NEntries = atoi(argv[1]);
  }

  ////////////////////////////////////////////////////////////

  // The order is maintained during rate calculation. 
  // I.e., keep the order and assign to any triggername!
  vector<TString> trignames; vector<int> prescales;

  // Store trignames and prescales in a map.
  // I.e., the order is not important.
  map<TString,int> map_TrigPrescls; 

  trignames.push_back("HLT1MuonLevel1");                   		 map_TrigPrescls["HLT1MuonLevel1"] = 1; //1000; 
  trignames.push_back("HLT1MuonIso");                          map_TrigPrescls["HLT1MuonIso"] = 1; 
  trignames.push_back("HLT1MuonNonIso");                       map_TrigPrescls["HLT1MuonNonIso"] = 1; 
  trignames.push_back("HLT2MuonIso");                      		 map_TrigPrescls["HLT2MuonIso"] = 1; 
  trignames.push_back("HLT2MuonNonIso");                       map_TrigPrescls["HLT2MuonNonIso"] = 1; 
  trignames.push_back("HLT2MuonJPsi");                         map_TrigPrescls["HLT2MuonJPsi"] = 1; 
  trignames.push_back("HLT2MuonUpsilon");                      map_TrigPrescls["HLT2MuonUpsilon"] = 1; 
  trignames.push_back("HLT2MuonZ");                            map_TrigPrescls["HLT2MuonZ"] = 1; 
  trignames.push_back("HLTNMuonNonIso");                       map_TrigPrescls["HLTNMuonNonIso"] = 1; 
  trignames.push_back("HLT2MuonSameSign");                     map_TrigPrescls["HLT2MuonSameSign"] = 1; 
  trignames.push_back("HLT1MuonPrescalePt3");              map_TrigPrescls["HLT1MuonPrescalePt3"] = 1; 
  trignames.push_back("HLT1MuonPrescalePt5");              map_TrigPrescls["HLT1MuonPrescalePt5"] = 1; 
  trignames.push_back("HLT1MuonPrescalePt7x7");            map_TrigPrescls["HLT1MuonPrescalePt7x7"] = 1; 
  trignames.push_back("HLT1MuonPrescalePt7x10");           map_TrigPrescls["HLT1MuonPrescalePt7x10"] = 1; 
  trignames.push_back("HLTB1JetMu");                           map_TrigPrescls["HLTB1JetMu"] = 1; 
  trignames.push_back("HLTB2JetMu");                           map_TrigPrescls["HLTB2JetMu"] = 1; 
  trignames.push_back("HLTB3JetMu");                           map_TrigPrescls["HLTB3JetMu"] = 1; 
  trignames.push_back("HLTB4JetMu");                           map_TrigPrescls["HLTB4JetMu"] = 1; 
  trignames.push_back("HLTBHTMu");                             map_TrigPrescls["HLTBHTMu"] = 1; 
  trignames.push_back("HLTBJPsiMuMu");                         map_TrigPrescls["HLTBJPsiMuMu"] = 1; 
  trignames.push_back("HLTXMuonBJet");                         map_TrigPrescls["HLTXMuonBJet"] = 1; 
  trignames.push_back("HLTXMuonBJetSoftMuon");                 map_TrigPrescls["HLTXMuonBJetSoftMuon"] = 1; 
  trignames.push_back("HLTXMuonJets");                         map_TrigPrescls["HLTXMuonJets"] = 1; 
  trignames.push_back("HLTXElectronMuon");                     map_TrigPrescls["HLTXElectronMuon"] = 1; 
  trignames.push_back("HLTXElectronMuonRelaxed");              map_TrigPrescls["HLTXElectronMuonRelaxed"] = 1; 
  trignames.push_back("HLTXMuonTau");                          map_TrigPrescls["HLTXMuonTau"] = 1; 

  trignames.push_back("HLT1jet");                              map_TrigPrescls["HLT1jet"] = 1; 
  trignames.push_back("HLT2jet");                              map_TrigPrescls["HLT2jet"] = 1; 
  trignames.push_back("HLT3jet");                              map_TrigPrescls["HLT3jet"] = 1; 
  trignames.push_back("HLT4jet");                              map_TrigPrescls["HLT4jet"] = 1; 
  trignames.push_back("HLT1MET");                              map_TrigPrescls["HLT1MET"] = 1; 
  trignames.push_back("HLT2jetAco");                           map_TrigPrescls["HLT2jetAco"] = 1; 
  trignames.push_back("HLT1jet1METAco");                       map_TrigPrescls["HLT1jet1METAco"] = 1; 
  trignames.push_back("HLT1jet1MET");                          map_TrigPrescls["HLT1jet1MET"] = 1; 
  trignames.push_back("HLT2jet1MET");                          map_TrigPrescls["HLT2jet1MET"] = 1; 
  trignames.push_back("HLT3jet1MET");                          map_TrigPrescls["HLT3jet1MET"] = 1; 
  trignames.push_back("HLT4jet1MET");                          map_TrigPrescls["HLT4jet1MET"] = 1; 
  trignames.push_back("HLT1MET1HT");                           map_TrigPrescls["HLT1MET1HT"] = 1; 
  trignames.push_back("HLT1SumET");                        map_TrigPrescls["HLT1SumET"] = 1; 
  trignames.push_back("HLT1jetPE7");                     map_TrigPrescls["HLT1jetPE7"] = 1; 
  trignames.push_back("HLT1jetPE5");                     map_TrigPrescls["HLT1jetPE5"] = 1; //10000; 
  trignames.push_back("HLT1jetPE3");                     map_TrigPrescls["HLT1jetPE3"] = 1; //100; 
  trignames.push_back("HLT1jetPE1");                     map_TrigPrescls["HLT1jetPE1"] = 1; 
  trignames.push_back("HLT1METPre3");                      map_TrigPrescls["HLT1METPre3"] = 1; 
  trignames.push_back("HLT1METPre2");                      map_TrigPrescls["HLT1METPre2"] = 1; 
  trignames.push_back("HLT1METPre1");                      map_TrigPrescls["HLT1METPre1"] = 1; //100; 
  trignames.push_back("HLT2jetAve30");                     map_TrigPrescls["HLT2jetAve30"] = 1; 
  trignames.push_back("HLT2jetAve60");                     map_TrigPrescls["HLT2jetAve60"] = 1; 
  trignames.push_back("HLT2jetAve110");                    map_TrigPrescls["HLT2jetAve110"] = 1; 
  trignames.push_back("HLT2jetAve150");                    map_TrigPrescls["HLT2jetAve150"] = 1; 
  trignames.push_back("HLT2jetAve200");                    map_TrigPrescls["HLT2jetAve200"] = 1; 
  trignames.push_back("HLT2jetvbfMET");                        map_TrigPrescls["HLT2jetvbfMET"] = 1; 
  trignames.push_back("HLTS2jet1METNV");                       map_TrigPrescls["HLTS2jet1METNV"] = 1; 
  trignames.push_back("HLTS2jet1METAco");                      map_TrigPrescls["HLTS2jet1METAco"] = 1; 
  trignames.push_back("HLTSjet1MET1Aco");                  map_TrigPrescls["HLTSjet1MET1Aco"] = 1; 
  trignames.push_back("HLTSjet2MET1Aco");                  map_TrigPrescls["HLTSjet2MET1Aco"] = 1; 
  trignames.push_back("HLTS2jetAco");                      map_TrigPrescls["HLTS2jetAco"] = 1; 
  trignames.push_back("HLTJetMETRapidityGap");             map_TrigPrescls["HLTJetMETRapidityGap"] = 1; //100; 

  trignames.push_back("HLT1Electron");                         map_TrigPrescls["HLT1Electron"] = 1; 
  trignames.push_back("HLT1ElectronRelaxed");                  map_TrigPrescls["HLT1ElectronRelaxed"] = 1; 
  trignames.push_back("HLT2Electron");                         map_TrigPrescls["HLT2Electron"] = 1; 
  trignames.push_back("HLT2ElectronRelaxed");                  map_TrigPrescls["HLT2ElectronRelaxed"] = 1; 

  trignames.push_back("HLT1Photon");                           map_TrigPrescls["HLT1Photon"] = 1; 
  trignames.push_back("HLT1PhotonRelaxed");                    map_TrigPrescls["HLT1PhotonRelaxed"] = 1; 
  trignames.push_back("HLT2Photon");                           map_TrigPrescls["HLT2Photon"] = 1; 
  trignames.push_back("HLT2PhotonRelaxed");                    map_TrigPrescls["HLT2PhotonRelaxed"] = 1; 
  trignames.push_back("HLT1EMHighEt");                         map_TrigPrescls["HLT1EMHighEt"] = 1; 
  trignames.push_back("HLT1EMVeryHighEt");                     map_TrigPrescls["HLT1EMVeryHighEt"] = 1; 
  trignames.push_back("HLT2ElectronZCounter");             map_TrigPrescls["HLT2ElectronZCounter"] = 1; 
  trignames.push_back("HLT2ElectronExclusive");            map_TrigPrescls["HLT2ElectronExclusive"] = 1; 
  trignames.push_back("HLT2PhotonExclusive");              map_TrigPrescls["HLT2PhotonExclusive"] = 1; 
  trignames.push_back("HLT1PhotonL1Isolated");             map_TrigPrescls["HLT1PhotonL1Isolated"] = 1; 

  trignames.push_back("HLTB1Jet");                             map_TrigPrescls["HLTB1Jet"] = 1; 
  trignames.push_back("HLTB2Jet");                             map_TrigPrescls["HLTB2Jet"] = 1; 
  trignames.push_back("HLTB3Jet");                             map_TrigPrescls["HLTB3Jet"] = 1; 
  trignames.push_back("HLTB4Jet");                             map_TrigPrescls["HLTB4Jet"] = 1; 
  trignames.push_back("HLTBHT");                               map_TrigPrescls["HLTBHT"] = 1; 
  trignames.push_back("HLT1Tau");                              map_TrigPrescls["HLT1Tau"] = 1; 
  trignames.push_back("HLT1Tau1MET");                          map_TrigPrescls["HLT1Tau1MET"] = 1; 
  trignames.push_back("HLT2TauPixel");                         map_TrigPrescls["HLT2TauPixel"] = 1; 
  trignames.push_back("HLTXElectronBJet");                     map_TrigPrescls["HLTXElectronBJet"] = 1; 
  trignames.push_back("HLTXElectron1Jet");                     map_TrigPrescls["HLTXElectron1Jet"] = 1; 
  trignames.push_back("HLTXElectron2Jet");                     map_TrigPrescls["HLTXElectron2Jet"] = 1; 
  trignames.push_back("HLTXElectron3Jet");                     map_TrigPrescls["HLTXElectron3Jet"] = 1; 
  trignames.push_back("HLTXElectron4Jet");                     map_TrigPrescls["HLTXElectron4Jet"] = 1; 
  trignames.push_back("HLTXElectronTau");                      map_TrigPrescls["HLTXElectronTau"] = 1; 

  trignames.push_back("HLTHcalIsolatedTrack");             map_TrigPrescls["HLTHcalIsolatedTrack"] = 1; 
  //trignames.push_back("HLTHcalIsolatedTrackNoEcalIsol"); //map_TrigPrescls["HLTHcalIsolatedTrackNoEcalIsol"] = 1;  
  trignames.push_back("HLTMinBiasPixel");                      map_TrigPrescls["HLTMinBiasPixel"] = 1; 
  //trignames.push_back("HLTMinBiasForAlignment");         //map_TrigPrescls["HLTMinBiasForAlignment"] = 1; 
  trignames.push_back("HLTMinBias");                           map_TrigPrescls["HLTMinBias"] = 1; 
  trignames.push_back("HLTZeroBias");                          map_TrigPrescls["HLTZeroBias"] = 1; 
  trignames.push_back("HLTMinBiasPixel");                      map_TrigPrescls["HLTMinBiasPixel"] = 1; 
  trignames.push_back("HLTMinBias");                           map_TrigPrescls["HLTMinBias"] = 1; 
  trignames.push_back("HLTZeroBias");                          map_TrigPrescls["HLTZeroBias"] = 1; 

  int Ntrig = (int)trignames.size();

  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////




  ////////////////////////////////////////////////////////////
  // Cross-sections [pb]
      
	vector<Double_t> xsec;
	vector<Double_t> skmeff; // Skim efficiencies

	//xsec.push_back(1.01E5); // PYTHIA cross-section for QCD in 170 < ^pt < 230
	//skmeff.push_back(1.);

	/*
	*/
	xsec.push_back(7.923E10); // PYTHIA cross-section for MinBias (pb)
	skmeff.push_back(0.105);  // skmeff for sbQCD 4e29
	//skmeff.push_back(0.00208);  // for sbQCD sample 1e32

	xsec.push_back(6.338E7); // PYTHIA cross-section times filter pp->muX (7.923E10*0.0008)
	skmeff.push_back(1.); // for ppMuX 4e29
	//skmeff.push_back(0.268); // for ppMuX sample 1e32

	xsec.push_back(7.685E8); // PYTHIA cross-section times filter pp->eleX (7.923E10*0.0097)
	skmeff.push_back(0.558); //for ppEleX 4e29
	//skmeff.push_back(0.0413); //for ppEleX sample 1e32

	// Add below the x-section of any other "signal" file to be read

  
  for (unsigned int ip = 0; ip < skmeff.size(); ip++) {
    xsec[ip] *= skmeff[ip];
  }
  
  // Convert cross-sections to cm^2
  for (unsigned int i = 0; i < skmeff.size(); i++){xsec[i] *= 1.E-36;}
  
  ////////////////////////////////////////////////////////////
  // Instanteneous Luminosity [cm^-2 s^-1]
  
  // For accurate rate calculation
  const double bunchCrossingTime = 25.0E-09;  // 25 ns
  const double maxFilledBunches = 3557;
  //const double nFilledBunches = 1000;
  
  
  /**** Different Beam conditions: ****/
 
  //const double ILumi = 1.E27;
  //const double nFilledBunches = 1;
  
  const double ILumi = 2E31;
  const double nFilledBunches = 156;
  
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
 
  double collisionRate = (nFilledBunches / maxFilledBunches) / bunchCrossingTime ;  // Hz
  
  ////////////////////////////////////////////////////////////
  // Files
  
  vector<TChain*> TabChain;
  vector<bool> doMuonCut; 
	vector<bool> doElecCut;

  vector<TString> ProcFil;
	ProcFil.push_back("/afs/hep.wisc.edu/cms/rekovic/3.8e29/HLTAnalyzer_QCDSingleBin_L1Skim_4e29_12kHz.root");
	ProcFil.push_back("/afs/hep.wisc.edu/cms/rekovic/3.8e29/HLTAnalyzer_ppMuX_4e29_12kHz.root");
	ProcFil.push_back("/afs/hep.wisc.edu/cms/rekovic/3.8e29/HLTAnalyzer_ppEleX_L1Skim_4e29_12kHz.root");




  for (unsigned int ipfile = 0; ipfile < ProcFil.size(); ipfile++){

  	TabChain.push_back(new TChain("HltTree"));
		TabChain[ipfile]->Add(ProcFil[ipfile]);
		//TabChain[0]->Add(ProcFil[ipfile]);

		// Set apropirate flags to avoid double counting
		if(ProcFil[ipfile].Contains("QCD") ) {

			cout << "File that requires muon/electron cut is " << ProcFil[ipfile] << endl;
 			doMuonCut.push_back(true); doElecCut.push_back(true);

		}
		else {
		
 			doMuonCut.push_back(false); doElecCut.push_back(false);

		} //end if


  }  // end for ipfile 


	cout << "Number of files (datasets) to process " << TabChain.size() << endl;


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
    cout<<"Available sample "<<ip<<", file " << ProcFil[ip] <<endl;
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
    cout<<"Processing bin "<<ip<<" ( "<< deno <<" events ) "<<", file " << ProcFil[ip] <<" (has "<<hltt[ip]->fChain->GetEntries()<<" events ) "<<endl;
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
  cout<<setprecision(3);

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



}


