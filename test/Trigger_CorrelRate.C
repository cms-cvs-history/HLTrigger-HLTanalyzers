/////////////////////////////////////////////////////////////////////////////////////////////////
//
//          Macro calculating overlaps between different triggers, HLT trigger individual- and
//          pure-rates, taking into account background and signal type samples
//
//          Send comments to bargassa@cern.ch
//
// Inputs : * Nfil : Number of files of different processes to be read
//          * ProcFil : Vector of files of different processes to be read (full path).
//                      These files are the outputs of HLTrigger/HLTanalyzers
//          * xsec : Cross-sections of different processes to be read (same order)
//          * Nqcd : Number of QCD pt-hat bins to be read
//          * ILumi : Instanteneous Luminosity
//          * Ntrig : Number of triggers
//          * trigname : Vector of triggers
//
// Outputs : * Matrix of global (averaged) overlaps between all triggers present
//           * Individual rate of each trigger (as if alone)
//           * Pure rate (additional rate given already present triggers)
//             taking into account overlaps in each process, including each QCD pt-hat bin
//           * Total (pure) rate
//
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <iomanip>
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


void Trigger_CorrelRate(){

const Int_t bin=100;
const Double_t min=0.;
const Double_t max=100.;
TH1F den("den","hlt obj",bin,min,max);
TH1F num("num","hlt obj",bin,min,max);
TH1F cum("cum","hlt obj",bin,min,max);

/////////////////////////////////////////////////////////////
// Cross-sections [pb]

const Int_t Nfil = 18; // Total number of files, of different x-section, to be read
const Int_t Nqcd = 16; // Actual number of QCD files
Double_t xsec[Nfil];
xsec[0] = 5.52E10; // PYTHIA cross-section for QCD in 0 < ^pt < 15
xsec[1] = 1.46E9; // PYTHIA cross-section for QCD in 15 < ^pt < 20
xsec[2] = 6.32E8; // PYTHIA cross-section for QCD in 20 < ^pt < 30
xsec[3] = 1.63E8; // PYTHIA cross-section for QCD in 30 < ^pt < 50
xsec[4] = 2.16E7; // PYTHIA cross-section for QCD in 50 < ^pt < 80
xsec[5] = 3.08E6; // PYTHIA cross-section for QCD in 80 < ^pt < 120
xsec[6] = 4.94E5; // PYTHIA cross-section for QCD in 120 < ^pt < 170
xsec[7] = 1.01E5; // PYTHIA cross-section for QCD in 170 < ^pt < 230
xsec[8] = 2.45E4; // PYTHIA cross-section for QCD in 230 < ^pt < 300
xsec[9] = 6.24E3; // PYTHIA cross-section for QCD in 300 < ^pt < 380
xsec[10] = 1.78E3; // PYTHIA cross-section for QCD in 380 < ^pt < 470
xsec[11] = 6.83E2; // PYTHIA cross-section for QCD in 470 < ^pt < 600
xsec[12] = 2.04E2; // PYTHIA cross-section for QCD in 600 < ^pt < 800
xsec[13] = 3.51E1; // PYTHIA cross-section for QCD in 800 < ^pt < 1000
xsec[14] = 1.09E1; // PYTHIA cross-section for QCD in 1000 < ^pt < 1400
xsec[15] = 1.06; // PYTHIA cross-section for QCD in 1400 < ^pt < 1800
// xsec[16] = 1.45E-1; // PYTHIA cross-section for QCD in 1800 < ^pt < 2200
// xsec[17] = 2.38E-2; // PYTHIA cross-section for QCD in 2200 < ^pt < 2600
// xsec[18] = 4.29E-3; // PYTHIA cross-section for QCD in 2600 < ^pt < 3000
// xsec[19] = 8.44E-4; // PYTHIA cross-section for QCD in 3000 < ^pt < 3500
// xsec[20] = 1.08E-4; // PYTHIA cross-section for QCD in 3500 < ^pt < ...
// Add below the x-section of any other "signal" file to be read
xsec[16] = 7.9E3; // PYTHIA cross-section for W -> e nu
xsec[17] = 9.8E3; // PYTHIA cross-section for W -> mu nu


// Convert cross-sections in cm^2
for (int i=0; i != Nfil; ++i){xsec[i] *= 1.E-36;}

/////////////////////////////////////////////////////////////
// Instanteneous Luminosity [cm^-2 s^-1]

// const Double_t ILumi = 8.E32;
const Double_t ILumi = 2.E33;
// const Double_t ILumi = 1.E34;

/////////////////////////////////////////////////////////////
// Files

TString ProcFil[Nfil];
ProcFil[0] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-0-15.root";
ProcFil[1] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-15-20.root";
ProcFil[2] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-20-30.root";
ProcFil[3] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-30-50.root";
ProcFil[4] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-50-80.root";
ProcFil[5] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-80-120.root";
ProcFil[6] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-120-170.root";
ProcFil[7] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-170-230.root";
ProcFil[8] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-230-300.root";
ProcFil[9] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-300-380.root";
ProcFil[10] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-380-470.root";
ProcFil[11] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-470-600.root";
ProcFil[12] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-600-800.root";
ProcFil[13] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-800-1000.root";
ProcFil[14] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-1000-1400.root";
ProcFil[15] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/1_3_0/HLTanalysis_QCD-1400-1800.root";
// ProcFil[16] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/HLTanalysis_QCD-1800-2200.root";
// ProcFil[17] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/HLTanalysis_QCD-2200-2600.root";
// ProcFil[18] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/HLTanalysis_QCD-2600-3000.root";
// ProcFil[19] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/HLTanalysis_QCD-3000-3500.root";
// ProcFil[20] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/HLTanalysis_QCD-3500.root";
// Add below any other "signal" file to be read
ProcFil[16] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/HLTanalysis_W_ENU.root";
ProcFil[17] = "/afs/cern.ch/user/b/bargassa/scratch0/HLTAnalysis_inputs/HLTanalysis_W_MUNU.root";

TChain TabChain[Nfil];
for (int ip = 0; ip != Nfil; ++ip){
  TabChain[ip] = new TChain("HltTree");
  TabChain[ip].Add(ProcFil[ip]);
}

/////////////////////////////////////////////////////////////
// Triggers

// Number of triggers
const Int_t Ntrig = 25;
TString trigname[Ntrig];

// Jet/MET triggers
trigname[0] = "HLT1jet";
trigname[1] = "HLT2jet";
trigname[2] = "HLT3jet";
trigname[3] = "HLT4jet";
trigname[4] = "HLT1MET";
trigname[5] = "HLT2jetAco";
trigname[6] = "HLT1jet1METAco";
trigname[7] = "HLT1jet1MET";
trigname[8] = "HLT2jet1MET";
trigname[9] = "HLT3jet1MET";
trigname[10] = "HLT4jet1MET";
trigname[11] = "HLT1MET1HT";
trigname[12] = "HLT1jetPE1";
trigname[13] = "HLT1jetPE3";
trigname[14] = "HLT1jetPE5";
// Egamma triggers
trigname[15] = "HLT1Electron";
trigname[16] = "HLT2Electron";
trigname[17] = "HLT2ElectronRelaxed";
trigname[18] = "HLT1Photon";
trigname[19] = "HLT2Photon";
trigname[20] = "HLT2PhotonRelaxed";
// Muon triggers
trigname[21] = "HLT1MuonIso";
trigname[22] = "HLT1MuonNonIso";
trigname[23] = "HLT2MuonIso";
trigname[24] = "HLT2MuonNonIso";
// Tau triggers
// trigname[25] = "HLT1Tau";
// trigname[26] = "HLT2Tau";
// trigname[27] = "HLTPixelTau";


// Trigger cuts
TCut hc[Ntrig],hnc[Ntrig];
for (int it=0; it != Ntrig; ++it){
  hc[it] = "TRIGG_"+trigname[it]+">0";
  hnc[it] = "TRIGG_"+trigname[it]+"==0";
}

/////////////////////////////////////////////////////////////
// Code...

Double_t Rat[Ntrig],sRat[Ntrig],pRat[Ntrig],spRat[Ntrig];
Int_t Oden[Ntrig];
Int_t Onum[Ntrig][Ntrig];
for (int it = 0; it != Ntrig; ++it){
  Rat[it] = 0.;
  sRat[it] = 0.;
  pRat[it] = 0.;
  spRat[it] = 0.;
  Oden[it] = 0;
  for (int jt = 0; jt != Ntrig; ++jt){Onum[it][jt] = 0;}
}

// Loop over files of different processes
for (int ip = 0; ip != Nfil; ++ip){

  TabChain[ip] . Draw("NobjHltPart>>den");                     // Get the denominators
  Int_t deno = den.GetEntries();                               //

  // Setup the flag for special particles generated
  Int_t BsFg = 0;
//   if (ip<Nqcd){
//     TH1F *mid = new TH1F("mid","MCid",100,0.,1000000000.);
//     TabChain[ip].Draw("MCPid>>mid");
//     Int_t *vid = TabChain[ip].GetV1();
//     for (int i=0; i<mid->GetEntries(); i++){
//       if ((vid[i]==-24)||(vid[i]==24)){BsFg = 1;}      // W
// //       if ((vid[i]==-23)||(vid[i]==23)){BsFg = 1;}      // Z
// //       if ((vid[i]==-11)||(vid[i]==11)){BsFg = 1;}      // Muon
//       else{BsFg = 0;}
//     }
//   }
//   else {BsFg = 0;}

  if (BsFg==0){

// Loop over triggers
    for (int it = 0; it != Ntrig; ++it){
      
      TabChain[ip] . Draw("NobjHltPart>>num",hc[it]);            // Get the numerators
      Int_t nume = num.GetEntries();                             //
      
      //test
      if (trigname[it]=="HLT1MuonIso"){
	cout << trigname[it] << " " << ip << " " << ILumi*xsec[ip]*eff(nume,deno) << " " 
	     << ILumi*xsec[ip]*seff(nume,deno) << "\n";
      }

      if (eff(nume,deno)>=2.E-4){                                // Not exposed to stat. fluctuation : fix the limit ... !
	Rat[it] += (xsec[ip]*eff(nume,deno));                      // Single rates
	sRat[it] += (pow(xsec[ip],2.)*pow(seff(nume,deno),2.));    //
	
	Oden[it] += nume;                                          // Get global overlap denominator
	TCut all = hc[it];
	for (int jt = 0; jt != Ntrig; ++jt){
	  if (jt==it){
	    Onum[it][jt] = Oden[it];                               // Get global overlap numerator 
	  }
	  else {
	    TCut hov = hc[it] + hc[jt];                            // Get global overlap numerator 
	    TabChain[ip] . Draw("NobjHltPart>>num",hov);           //
	    Onum[it][jt] += num.GetEntries();                      //
	  }
	  if (jt<it) {all += hnc[jt];}
	}
	TabChain[ip] . Draw("NobjHltPart>>cum",all);               // Get Trigger[it] pure rate...
	int cumu = cum.GetEntries();                               // ...
	pRat[it] += (xsec[ip]*eff(cumu,deno));                     // ... summed over different processes
	spRat[it] += (pow(xsec[ip],2.)*pow(seff(cumu,deno),2.));   //
      }
      
    } // End loop over triggers

  }

}

Double_t RTOT = 0.;
Double_t sRTOT = 0.;
// Loop over triggers
for (int it = 0; it != Ntrig; ++it){
  RTOT += pRat[it];                                            // Total Rate
  sRTOT += spRat[it];
}
sRTOT = sqrt(sRTOT);

/////////////////////////////////////////////////////////////
// Results

cout.setf(ios::floatfield,ios::fixed);
cout<<setprecision(2);

cout << "\n";
cout << "Trigger global overlaps : " << "\n";
cout << "      ";
for (int it = 0; it != Ntrig; ++it){
  if (it<=9) {cout << "T" << it << "    ";}
  else {cout << "T" << it << "   ";}
  if (it==Ntrig-1) {cout << "\n";}
}
for (int it = 0; it != Ntrig; ++it){
  if (it<=9) {cout << "T" << it << "  : ";}
  else {cout << "T" << it << " : ";}
  for (int jt = 0; jt != Ntrig; ++jt){
    if (jt>=it) {cout << eff(Onum[it][jt],Oden[jt]) << "  ";}  // Overlap O(ij) = T(i) x T(j) / T(j)
    else {cout << "      ";}
  }
  cout << "\n";
}

cout << "\n";
cout << "Trigger Rates [Hz] : " << "\n";
const Int_t esp = 5;
for (int it=0; it != Ntrig; ++it){
  cout << setw(20) << trigname[it] << " : Individual : " << setw(5) << ILumi*Rat[it] << " +- " << 
    ILumi*sqrt(sRat[it]) << " Pure : " << setw(5) << ILumi*pRat[it] << "\n";
}
cout << "\n";
cout << setw(57) << "TOTAL RATE : " << setw(5) << ILumi*RTOT << " +- " << ILumi*sRTOT << " Hz" << "\n";
cout << "\n";

}
