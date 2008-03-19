//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 12 10:37:25 2007 by ROOT version 5.14/00f
// from TTree HltTree/
// found on file: /tmp/test.root
//////////////////////////////////////////////////////////

#ifndef OHltTree_h
#define OHltTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class OHltTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           NrecoJetCal;
   Int_t           NrecoJetGen;
   Int_t           NrecoTowCal;
   Float_t         recoJetCalPt[1000];   //[NrecoJetCal]
   Float_t         recoJetCalPhi[1000];   //[NrecoJetCal]
   Float_t         recoJetCalEta[1000];   //[NrecoJetCal]
   Float_t         recoJetCalEt[1000];   //[NrecoJetCal]
   Float_t         recoJetCalE[1000];   //[NrecoJetCal]
   Float_t         recoJetGenPt[1000];   //[NrecoJetGen]
   Float_t         recoJetGenPhi[1000];   //[NrecoJetGen]
   Float_t         recoJetGenEta[1000];   //[NrecoJetGen]
   Float_t         recoJetGenEt[1000];   //[NrecoJetGen]
   Float_t         recoJetGenE[1000];   //[NrecoJetGen]
   Float_t         recoTowEt[1000];   //[NrecoTowCal]
   Float_t         recoTowEta[1000];   //[NrecoTowCal]
   Float_t         recoTowPhi[1000];   //[NrecoTowCal]
   Float_t         recoTowE[1000];   //[NrecoTowCal]
   Float_t         recoTowEm[1000];   //[NrecoTowCal]
   Float_t         recoTowHad[1000];   //[NrecoTowCal]
   Float_t         recoTowOE[1000];   //[NrecoTowCal]
   Float_t         recoMetCal;
   Float_t         recoMetCalPhi;
   Float_t         recoMetCalSum;
   Float_t         recoMetGen;
   Float_t         recoMetGenPhi;
   Float_t         recoMetGenSum;
   Float_t         recoHTCal;
   Float_t         recoHTCalPhi;
   Float_t         recoHTCalSum;
   Int_t           NohTau2;
   Float_t         ohTau2Eiso[1000];   //[NohTau2]
   Float_t         ohTau2L2Tpt[1000];   //[NohTau2]
   Int_t           ohTau2L2Tiso[1000];   //[NohTau2]
   Int_t           NohTau1;
   Float_t         ohTau1Eiso[1000];   //[NohTau1]
   Float_t         ohTau1L2Tpt[1000];   //[NohTau1]
   Int_t           ohTau1L2Tiso[1000];   //[NohTau1]
   Float_t         ohTau1L3Tpt[1000];   //[NohTau1]
   Int_t           ohTau1L3Tiso[1000];   //[NohTau1]
   Int_t           NrecoElec;
   Float_t         recoElecPt[1000];   //[NrecoElec]
   Float_t         recoElecPhi[1000];   //[NrecoElec]
   Float_t         recoElecEta[1000];   //[NrecoElec]
   Float_t         recoElecEt[1000];   //[NrecoElec]
   Float_t         recoElecE[1000];   //[NrecoElec]
   Int_t           NrecoPhot;
   Float_t         recoPhotPt[1000];   //[NrecoPhot]
   Float_t         recoPhotPhi[1000];   //[NrecoPhot]
   Float_t         recoPhotEta[1000];   //[NrecoPhot]
   Float_t         recoPhotEt[1000];   //[NrecoPhot]
   Float_t         recoPhotE[1000];   //[NrecoPhot]
   Int_t           NohPhot;
   Float_t         ohPhotEt[1000];   //[NohPhot]
   Float_t         ohPhotEta[1000];   //[NohPhot]
   Float_t         ohPhotPhi[1000];   //[NohPhot]
   Float_t         ohPhotEiso[1000];   //[NohPhot]
   Float_t         ohPhotHiso[1000];   //[NohPhot]
   Float_t         ohPhotTiso[1000];   //[NohPhot]
   Int_t           ohPhotL1iso[1000];   //[NohPhot]
   Int_t           NohEle;
   Float_t         ohEleEt[1000];   //[NohEle]
   Float_t         ohEleEta[1000];   //[NohEle]
   Float_t         ohElePhi[1000];   //[NohEle]
   Float_t         ohEleE[1000];   //[NohEle]
   Float_t         ohEleP[1000];   //[NohEle]
   Float_t         ohEleHiso[1000];   //[NohEle]
   Float_t         ohEleTiso[1000];   //[NohEle]
   Int_t           ohEleL1iso[1000];   //[NohEle]
   Int_t           NrecoMuon;
   Float_t         recoMuonPt[1000];   //[NrecoMuon]
   Float_t         recoMuonPhi[1000];   //[NrecoMuon]
   Float_t         recoMuonEta[1000];   //[NrecoMuon]
   Float_t         recoMuonEt[1000];   //[NrecoMuon]
   Float_t         recoMuonE[1000];   //[NrecoMuon]
   Int_t           NohMuL2;
   Float_t         ohMuL2Pt[1000];   //[NohMuL2]
   Float_t         ohMuL2PtErr[1000];   //[NohMuL2]
   Float_t         ohMuL2Phi[1000];   //[NohMuL2]
   Float_t         ohMuL2Eta[1000];   //[NohMuL2]
   Int_t           ohMuL2Chg[1000];   //[NohMuL2]
   Int_t           ohMuL2Iso[1000];   //[NohMuL2]
   Float_t         ohMuL2Dr[1000];   //[NohMuL2]
   Float_t         ohMuL2Dz[1000];   //[NohMuL2]
   Int_t           NohMuL3;
   Float_t         ohMuL3Pt[1000];   //[NohMuL3]
   Float_t         ohMuL3PtErr[1000];   //[NohMuL3]
   Float_t         ohMuL3Phi[1000];   //[NohMuL3]
   Float_t         ohMuL3Eta[1000];   //[NohMuL3]
   Int_t           ohMuL3Chg[1000];   //[NohMuL3]
   Int_t           ohMuL3Iso[1000];   //[NohMuL3]
   Float_t         ohMuL3Dr[1000];   //[NohMuL3]
   Float_t         ohMuL3Dz[1000];   //[NohMuL3]
   Int_t           NMCpart;
   Int_t           MCpid[1000];   //[NMCpart]
   Float_t         MCvtxX[1000];   //[NMCpart]
   Float_t         MCvtxY[1000];   //[NMCpart]
   Float_t         MCvtxZ[1000];   //[NMCpart]
   Float_t         MCpt[1000];   //[NMCpart]
   Float_t         MCeta[1000];   //[NMCpart]
   Float_t         MCphi[1000];   //[NMCpart]
   Float_t         MCPtHat;
   Int_t           MCmu3;
   Int_t           MCel1;
   Int_t           MCbb;
   Int_t           MCab;
   Int_t           NL1IsolEm;
   Float_t         L1IsolEmEt[1000];   //[NL1IsolEm]
   Float_t         L1IsolEmE[1000];   //[NL1IsolEm]
   Float_t         L1IsolEmEta[1000];   //[NL1IsolEm]
   Float_t         L1IsolEmPhi[1000];   //[NL1IsolEm]
   Int_t           NL1NIsolEm;
   Float_t         L1NIsolEmEt[1000];   //[NL1NIsolEm]
   Float_t         L1NIsolEmE[1000];   //[NL1NIsolEm]
   Float_t         L1NIsolEmEta[1000];   //[NL1NIsolEm]
   Float_t         L1NIsolEmPhi[1000];   //[NL1NIsolEm]
   Int_t           NL1Mu;
   Float_t         L1MuPt[1000];   //[NL1Mu]
   Float_t         L1MuE[1000];   //[NL1Mu]
   Float_t         L1MuEta[1000];   //[NL1Mu]
   Float_t         L1MuPhi[1000];   //[NL1Mu]
   Int_t           L1MuIsol[1000];   //[NL1Mu]
   Int_t           L1MuMip[1000];   //[NL1Mu]
   Int_t           L1MuFor[1000];   //[NL1Mu]
   Int_t           L1MuRPC[1000];   //[NL1Mu]
   Int_t           L1MuQal[1000];   //[NL1Mu]
   Int_t           NL1CenJet;
   Float_t         L1CenJetEt[1000];   //[NL1CenJet]
   Float_t         L1CenJetE[1000];   //[NL1CenJet]
   Float_t         L1CenJetEta[1000];   //[NL1CenJet]
   Float_t         L1CenJetPhi[1000];   //[NL1CenJet]
   Int_t           NL1ForJet;
   Float_t         L1ForJetEt[1000];   //[NL1ForJet]
   Float_t         L1ForJetE[1000];   //[NL1ForJet]
   Float_t         L1ForJetEta[1000];   //[NL1ForJet]
   Float_t         L1ForJetPhi[1000];   //[NL1ForJet]
   Int_t           NL1Tau;
   Float_t         L1TauEt[1000];   //[NL1Tau]
   Float_t         L1TauE[1000];   //[NL1Tau]
   Float_t         L1TauEta[1000];   //[NL1Tau]
   Float_t         L1TauPhi[1000];   //[NL1Tau]
   Float_t         L1Met;
   Float_t         L1MetPhi;
   Float_t         L1MetTot;
   Float_t         L1MetHad;
   Int_t           run;
   Int_t           event;
   Int_t           HLT1jet;
   Int_t           HLT2jet;
   Int_t           HLT3jet;
   Int_t           HLT4jet;
   Int_t           HLT1MET;
   Int_t           HLT2jetAco;
   Int_t           HLT1jet1METAco;
   Int_t           HLT1jet1MET;
   Int_t           HLT2jet1MET;
   Int_t           HLT3jet1MET;
   Int_t           HLT4jet1MET;
   Int_t           HLT1MET1HT;
   Int_t           CandHLT1SumET;
   Int_t           HLT1jetPE1;
   Int_t           HLT1jetPE3;
   Int_t           HLT1jetPE5;
   Int_t           CandHLT1jetPE7;
   Int_t           CandHLT1METPre1;
   Int_t           CandHLT1METPre2;
   Int_t           CandHLT1METPre3;
   Int_t           CandHLT2jetAve30;
   Int_t           CandHLT2jetAve60;
   Int_t           CandHLT2jetAve110;
   Int_t           CandHLT2jetAve150;
   Int_t           CandHLT2jetAve200;
   Int_t           HLT2jetvbfMET;
   Int_t           HLTS2jet1METNV;
   Int_t           HLTS2jet1METAco;
   Int_t           CandHLTSjet1MET1Aco;
   Int_t           CandHLTSjet2MET1Aco;
   Int_t           CandHLTS2jetAco;
   Int_t           CandHLTJetMETRapidityGap;
   Int_t           HLT1Electron;
   Int_t           HLT1ElectronRelaxed;
   Int_t           HLT2Electron;
   Int_t           HLT2ElectronRelaxed;
   Int_t           HLT1Photon;
   Int_t           HLT1PhotonRelaxed;
   Int_t           HLT2Photon;
   Int_t           HLT2PhotonRelaxed;
   Int_t           HLT1EMHighEt;
   Int_t           HLT1EMVeryHighEt;
   Int_t           CandHLT2ElectronZCounter;
   Int_t           CandHLT2ElectronExclusive;
   Int_t           CandHLT2PhotonExclusive;
   Int_t           CandHLT1PhotonL1Isolated;
   Int_t           HLT1MuonIso;
   Int_t           HLT1MuonNonIso;
   Int_t           CandHLT2MuonIso;
   Int_t           HLT2MuonNonIso;
   Int_t           HLT2MuonJPsi;
   Int_t           HLT2MuonUpsilon;
   Int_t           HLT2MuonZ;
   Int_t           HLTNMuonNonIso;
   Int_t           HLT2MuonSameSign;
   Int_t           CandHLT1MuonPrescalePt3;
   Int_t           CandHLT1MuonPrescalePt5;
   Int_t           CandHLT1MuonPrescalePt7x7;
   Int_t           CandHLT1MuonPrescalePt7x10;
   Int_t           CandHLT1MuonLevel1;
   Int_t           HLTB1Jet;
   Int_t           HLTB2Jet;
   Int_t           HLTB3Jet;
   Int_t           HLTB4Jet;
   Int_t           HLTBHT;
   Int_t           HLTB1JetMu;
   Int_t           HLTB2JetMu;
   Int_t           HLTB3JetMu;
   Int_t           HLTB4JetMu;
   Int_t           HLTBHTMu;
   Int_t           HLTBJPsiMuMu;
   Int_t           HLT1Tau;
   Int_t           HLT1Tau1MET;
   Int_t           HLT2TauPixel;
   Int_t           HLTXElectronBJet;
   Int_t           HLTXMuonBJet;
   Int_t           HLTXMuonBJetSoftMuon;
   Int_t           HLTXElectron1Jet;
   Int_t           HLTXElectron2Jet;
   Int_t           HLTXElectron3Jet;
   Int_t           HLTXElectron4Jet;
   Int_t           HLTXMuonJets;
   Int_t           HLTXElectronMuon;
   Int_t           HLTXElectronMuonRelaxed;
   Int_t           HLTXElectronTau;
   Int_t           HLTXMuonTau;
   Int_t           CandHLTHcalIsolatedTrack;
   Int_t           HLTMinBiasPixel;
   Int_t           HLTMinBias;
   Int_t           HLTZeroBias;
   Int_t           L1_SingleMu3;
   Int_t           L1_SingleMu5;
   Int_t           L1_SingleMu7;
   Int_t           L1_SingleMu10;
   Int_t           L1_SingleMu14;
   Int_t           L1_SingleMu20;
   Int_t           L1_SingleMu25;
   Int_t           L1_SingleIsoEG5;
   Int_t           L1_SingleIsoEG8;
   Int_t           L1_SingleIsoEG10;
   Int_t           L1_SingleIsoEG12;
   Int_t           L1_SingleIsoEG15;
   Int_t           L1_SingleIsoEG20;
   Int_t           L1_SingleIsoEG25;
   Int_t           L1_SingleEG5;
   Int_t           L1_SingleEG8;
   Int_t           L1_SingleEG10;
   Int_t           L1_SingleEG12;
   Int_t           L1_SingleEG15;
   Int_t           L1_SingleEG20;
   Int_t           L1_SingleEG25;
   Int_t           L1_SingleJet15;
   Int_t           L1_SingleJet20;
   Int_t           L1_SingleJet30;
   Int_t           L1_SingleJet50;
   Int_t           L1_SingleJet70;
   Int_t           L1_SingleJet100;
   Int_t           L1_SingleJet150;
   Int_t           L1_SingleJet200;
   Int_t           L1_SingleTauJet10;
   Int_t           L1_SingleTauJet20;
   Int_t           L1_SingleTauJet30;
   Int_t           L1_SingleTauJet35;
   Int_t           L1_SingleTauJet40;
   Int_t           L1_SingleTauJet60;
   Int_t           L1_SingleTauJet80;
   Int_t           L1_SingleTauJet100;
   Int_t           L1_HTT100;
   Int_t           L1_HTT200;
   Int_t           L1_HTT250;
   Int_t           L1_HTT300;
   Int_t           L1_HTT400;
   Int_t           L1_HTT500;
   Int_t           L1_ETM20;
   Int_t           L1_ETM30;
   Int_t           L1_ETM40;
   Int_t           L1_ETM50;
   Int_t           L1_ETM60;
   Int_t           L1_DoubleMu3;
   Int_t           L1_DoubleIsoEG8;
   Int_t           L1_DoubleIsoEG10;
   Int_t           L1_DoubleEG5;
   Int_t           L1_DoubleEG10;
   Int_t           L1_DoubleEG15;
   Int_t           L1_DoubleJet70;
   Int_t           L1_DoubleJet100;
   Int_t           L1_DoubleTauJet20;
   Int_t           L1_DoubleTauJet30;
   Int_t           L1_DoubleTauJet35;
   Int_t           L1_DoubleTauJet40;
   Int_t           L1_Mu3_IsoEG5;
   Int_t           L1_Mu5_IsoEG10;
   Int_t           L1_Mu3_EG12;
   Int_t           L1_Mu3_Jet15;
   Int_t           L1_Mu5_Jet15;
   Int_t           L1_Mu3_Jet70;
   Int_t           L1_Mu5_Jet20;
   Int_t           L1_Mu5_TauJet20;
   Int_t           L1_Mu5_TauJet30;
   Int_t           L1_IsoEG10_EG10;
   Int_t           L1_IsoEG10_Jet15;
   Int_t           L1_IsoEG10_Jet20;
   Int_t           L1_IsoEG10_Jet30;
   Int_t           L1_IsoEG10_Jet70;
   Int_t           L1_IsoEG10_TauJet20;
   Int_t           L1_IsoEG10_TauJet30;
   Int_t           L1_EG10_Jet15;
   Int_t           L1_EG12_Jet20;
   Int_t           L1_EG12_Jet70;
   Int_t           L1_EG12_TauJet40;
   Int_t           L1_Jet70_TauJet40;
   Int_t           L1_Mu3_HTT200;
   Int_t           L1_IsoEG10_HTT200;
   Int_t           L1_EG12_HTT200;
   Int_t           L1_Jet70_HTT200;
   Int_t           L1_TauJet40_HTT200;
   Int_t           L1_Mu3_ETM30;
   Int_t           L1_IsoEG10_ETM30;
   Int_t           L1_EG12_ETM30;
   Int_t           L1_Jet70_ETM40;
   Int_t           L1_TauJet20_ETM20;
   Int_t           L1_TauJet30_ETM30;
   Int_t           L1_TauJet30_ETM40;
   Int_t           L1_HTT100_ETM30;
   Int_t           L1_TripleMu3;
   Int_t           L1_TripleIsoEG5;
   Int_t           L1_TripleEG10;
   Int_t           L1_TripleJet50;
   Int_t           L1_TripleTauJet40;
   Int_t           L1_DoubleMu3_IsoEG5;
   Int_t           L1_DoubleMu3_EG10;
   Int_t           L1_DoubleIsoEG5_Mu3;
   Int_t           L1_DoubleEG10_Mu3;
   Int_t           L1_DoubleMu3_HTT200;
   Int_t           L1_DoubleIsoEG5_HTT200;
   Int_t           L1_DoubleEG10_HTT200;
   Int_t           L1_DoubleJet50_HTT200;
   Int_t           L1_DoubleTauJet40_HTT200;
   Int_t           L1_DoubleMu3_ETM20;
   Int_t           L1_DoubleIsoEG5_ETM20;
   Int_t           L1_DoubleEG10_ETM20;
   Int_t           L1_DoubleJet50_ETM20;
   Int_t           L1_DoubleTauJet40_ETM20;
   Int_t           L1_QuadJet30;
   Int_t           L1_ExclusiveDoubleIsoEG6;
   Int_t           L1_ExclusiveDoubleJet60;
   Int_t           L1_ExclusiveJet25_Gap_Jet25;
   Int_t           L1_IsoEG10_Jet20_ForJet10;
   Int_t           L1_MinBias_HTT10;
   Int_t           L1_ZeroBias;

   // List of branches
   TBranch        *b_NrecoJetCal;   //!
   TBranch        *b_NrecoJetGen;   //!
   TBranch        *b_NrecoTowCal;   //!
   TBranch        *b_recoJetCalPt;   //!
   TBranch        *b_recoJetCalPhi;   //!
   TBranch        *b_recoJetCalEta;   //!
   TBranch        *b_recoJetCalEt;   //!
   TBranch        *b_recoJetCalE;   //!
   TBranch        *b_recoJetGenPt;   //!
   TBranch        *b_recoJetGenPhi;   //!
   TBranch        *b_recoJetGenEta;   //!
   TBranch        *b_recoJetGenEt;   //!
   TBranch        *b_recoJetGenE;   //!
   TBranch        *b_recoTowEt;   //!
   TBranch        *b_recoTowEta;   //!
   TBranch        *b_recoTowPhi;   //!
   TBranch        *b_recoTowE;   //!
   TBranch        *b_recoTowEm;   //!
   TBranch        *b_recoTowHad;   //!
   TBranch        *b_recoTowOE;   //!
   TBranch        *b_recoMetCal;   //!
   TBranch        *b_recoMetCalPhi;   //!
   TBranch        *b_recoMetCalSum;   //!
   TBranch        *b_recoMetGen;   //!
   TBranch        *b_recoMetGenPhi;   //!
   TBranch        *b_recoMetGenSum;   //!
   TBranch        *b_recoHTCal;   //!
   TBranch        *b_recoHTCalPhi;   //!
   TBranch        *b_recoHTCalSum;   //!
   TBranch        *b_NohTau2;   //!
   TBranch        *b_ohTau2Eiso;   //!
   TBranch        *b_ohTau2L2Tpt;   //!
   TBranch        *b_ohTau2L2Tiso;   //!
   TBranch        *b_NohTau1;   //!
   TBranch        *b_ohTau1Eiso;   //!
   TBranch        *b_ohTau1L2Tpt;   //!
   TBranch        *b_ohTau1L2Tiso;   //!
   TBranch        *b_ohTau1L3Tpt;   //!
   TBranch        *b_ohTau1L3Tiso;   //!
   TBranch        *b_NrecoElec;   //!
   TBranch        *b_recoElecPt;   //!
   TBranch        *b_recoElecPhi;   //!
   TBranch        *b_recoElecEta;   //!
   TBranch        *b_recoElecEt;   //!
   TBranch        *b_recoElecE;   //!
   TBranch        *b_NrecoPhot;   //!
   TBranch        *b_recoPhotPt;   //!
   TBranch        *b_recoPhotPhi;   //!
   TBranch        *b_recoPhotEta;   //!
   TBranch        *b_recoPhotEt;   //!
   TBranch        *b_recoPhotE;   //!
   TBranch        *b_NohPhot;   //!
   TBranch        *b_ohPhotEt;   //!
   TBranch        *b_ohPhotEta;   //!
   TBranch        *b_ohPhotPhi;   //!
   TBranch        *b_ohPhotEiso;   //!
   TBranch        *b_ohPhotHiso;   //!
   TBranch        *b_ohPhotTiso;   //!
   TBranch        *b_ohPhotL1iso;   //!
   TBranch        *b_NohEle;   //!
   TBranch        *b_ohEleEt;   //!
   TBranch        *b_ohEleEta;   //!
   TBranch        *b_ohElePhi;   //!
   TBranch        *b_ohEleE;   //!
   TBranch        *b_ohEleP;   //!
   TBranch        *b_ohEleHiso;   //!
   TBranch        *b_ohEleTiso;   //!
   TBranch        *b_ohEleL1iso;   //!
   TBranch        *b_NrecoMuon;   //!
   TBranch        *b_recoMuonPt;   //!
   TBranch        *b_recoMuonPhi;   //!
   TBranch        *b_recoMuonEta;   //!
   TBranch        *b_recoMuonEt;   //!
   TBranch        *b_recoMuonE;   //!
   TBranch        *b_NohMuL2;   //!
   TBranch        *b_ohMuL2Pt;   //!
   TBranch        *b_ohMuL2PtErr;   //!
   TBranch        *b_ohMuL2Phi;   //!
   TBranch        *b_ohMuL2Eta;   //!
   TBranch        *b_ohMuL2Chg;   //!
   TBranch        *b_ohMuL2Iso;   //!
   TBranch        *b_ohMuL2Dr;   //!
   TBranch        *b_ohMuL2Dz;   //!
   TBranch        *b_NohMuL3;   //!
   TBranch        *b_ohMuL3Pt;   //!
   TBranch        *b_ohMuL3PtErr;   //!
   TBranch        *b_ohMuL3Phi;   //!
   TBranch        *b_ohMuL3Eta;   //!
   TBranch        *b_ohMuL3Chg;   //!
   TBranch        *b_ohMuL3Iso;   //!
   TBranch        *b_ohMuL3Dr;   //!
   TBranch        *b_ohMuL3Dz;   //!
   TBranch        *b_NMCpart;   //!
   TBranch        *b_MCpid;   //!
   TBranch        *b_MCvtxX;   //!
   TBranch        *b_MCvtxY;   //!
   TBranch        *b_MCvtxZ;   //!
   TBranch        *b_MCpt;   //!
   TBranch        *b_MCeta;   //!
   TBranch        *b_MCphi;   //!
   TBranch        *b_MCPtHat;   //!
   TBranch        *b_MCmu3;   //!
   TBranch        *b_MCel1;   //!
   TBranch        *b_MCbb;   //!
   TBranch        *b_MCab;   //!
   TBranch        *b_NL1IsolEm;   //!
   TBranch        *b_L1IsolEmEt;   //!
   TBranch        *b_L1IsolEmE;   //!
   TBranch        *b_L1IsolEmEta;   //!
   TBranch        *b_L1IsolEmPhi;   //!
   TBranch        *b_NL1NIsolEm;   //!
   TBranch        *b_L1NIsolEmEt;   //!
   TBranch        *b_L1NIsolEmE;   //!
   TBranch        *b_L1NIsolEmEta;   //!
   TBranch        *b_L1NIsolEmPhi;   //!
   TBranch        *b_NL1Mu;   //!
   TBranch        *b_L1MuPt;   //!
   TBranch        *b_L1MuE;   //!
   TBranch        *b_L1MuEta;   //!
   TBranch        *b_L1MuPhi;   //!
   TBranch        *b_L1MuIsol;   //!
   TBranch        *b_L1MuMip;   //!
   TBranch        *b_L1MuFor;   //!
   TBranch        *b_L1MuRPC;   //!
   TBranch        *b_L1MuQal;   //!
   TBranch        *b_NL1CenJet;   //!
   TBranch        *b_L1CenJetEt;   //!
   TBranch        *b_L1CenJetE;   //!
   TBranch        *b_L1CenJetEta;   //!
   TBranch        *b_L1CenJetPhi;   //!
   TBranch        *b_NL1ForJet;   //!
   TBranch        *b_L1ForJetEt;   //!
   TBranch        *b_L1ForJetE;   //!
   TBranch        *b_L1ForJetEta;   //!
   TBranch        *b_L1ForJetPhi;   //!
   TBranch        *b_NL1Tau;   //!
   TBranch        *b_L1TauEt;   //!
   TBranch        *b_L1TauE;   //!
   TBranch        *b_L1TauEta;   //!
   TBranch        *b_L1TauPhi;   //!
   TBranch        *b_L1Met;   //!
   TBranch        *b_L1MetPhi;   //!
   TBranch        *b_L1MetTot;   //!
   TBranch        *b_L1MetHad;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_HLT1jet;   //!
   TBranch        *b_HLT2jet;   //!
   TBranch        *b_HLT3jet;   //!
   TBranch        *b_HLT4jet;   //!
   TBranch        *b_HLT1MET;   //!
   TBranch        *b_HLT2jetAco;   //!
   TBranch        *b_HLT1jet1METAco;   //!
   TBranch        *b_HLT1jet1MET;   //!
   TBranch        *b_HLT2jet1MET;   //!
   TBranch        *b_HLT3jet1MET;   //!
   TBranch        *b_HLT4jet1MET;   //!
   TBranch        *b_HLT1MET1HT;   //!
   TBranch        *b_CandHLT1SumET;   //!
   TBranch        *b_HLT1jetPE1;   //!
   TBranch        *b_HLT1jetPE3;   //!
   TBranch        *b_HLT1jetPE5;   //!
   TBranch        *b_CandHLT1jetPE7;   //!
   TBranch        *b_CandHLT1METPre1;   //!
   TBranch        *b_CandHLT1METPre2;   //!
   TBranch        *b_CandHLT1METPre3;   //!
   TBranch        *b_CandHLT2jetAve30;   //!
   TBranch        *b_CandHLT2jetAve60;   //!
   TBranch        *b_CandHLT2jetAve110;   //!
   TBranch        *b_CandHLT2jetAve150;   //!
   TBranch        *b_CandHLT2jetAve200;   //!
   TBranch        *b_HLT2jetvbfMET;   //!
   TBranch        *b_HLTS2jet1METNV;   //!
   TBranch        *b_HLTS2jet1METAco;   //!
   TBranch        *b_CandHLTSjet1MET1Aco;   //!
   TBranch        *b_CandHLTSjet2MET1Aco;   //!
   TBranch        *b_CandHLTS2jetAco;   //!
   TBranch        *b_CandHLTJetMETRapidityGap;   //!
   TBranch        *b_HLT1Electron;   //!
   TBranch        *b_HLT1ElectronRelaxed;   //!
   TBranch        *b_HLT2Electron;   //!
   TBranch        *b_HLT2ElectronRelaxed;   //!
   TBranch        *b_HLT1Photon;   //!
   TBranch        *b_HLT1PhotonRelaxed;   //!
   TBranch        *b_HLT2Photon;   //!
   TBranch        *b_HLT2PhotonRelaxed;   //!
   TBranch        *b_HLT1EMHighEt;   //!
   TBranch        *b_HLT1EMVeryHighEt;   //!
   TBranch        *b_CandHLT2ElectronZCounter;   //!
   TBranch        *b_CandHLT2ElectronExclusive;   //!
   TBranch        *b_CandHLT2PhotonExclusive;   //!
   TBranch        *b_CandHLT1PhotonL1Isolated;   //!
   TBranch        *b_HLT1MuonIso;   //!
   TBranch        *b_HLT1MuonNonIso;   //!
   TBranch        *b_CandHLT2MuonIso;   //!
   TBranch        *b_HLT2MuonNonIso;   //!
   TBranch        *b_HLT2MuonJPsi;   //!
   TBranch        *b_HLT2MuonUpsilon;   //!
   TBranch        *b_HLT2MuonZ;   //!
   TBranch        *b_HLTNMuonNonIso;   //!
   TBranch        *b_HLT2MuonSameSign;   //!
   TBranch        *b_CandHLT1MuonPrescalePt3;   //!
   TBranch        *b_CandHLT1MuonPrescalePt5;   //!
   TBranch        *b_CandHLT1MuonPrescalePt7x7;   //!
   TBranch        *b_CandHLT1MuonPrescalePt7x10;   //!
   TBranch        *b_CandHLT1MuonLevel1;   //!
   TBranch        *b_HLTB1Jet;   //!
   TBranch        *b_HLTB2Jet;   //!
   TBranch        *b_HLTB3Jet;   //!
   TBranch        *b_HLTB4Jet;   //!
   TBranch        *b_HLTBHT;   //!
   TBranch        *b_HLTB1JetMu;   //!
   TBranch        *b_HLTB2JetMu;   //!
   TBranch        *b_HLTB3JetMu;   //!
   TBranch        *b_HLTB4JetMu;   //!
   TBranch        *b_HLTBHTMu;   //!
   TBranch        *b_HLTBJPsiMuMu;   //!
   TBranch        *b_HLT1Tau;   //!
   TBranch        *b_HLT1Tau1MET;   //!
   TBranch        *b_HLT2TauPixel;   //!
   TBranch        *b_HLTXElectronBJet;   //!
   TBranch        *b_HLTXMuonBJet;   //!
   TBranch        *b_HLTXMuonBJetSoftMuon;   //!
   TBranch        *b_HLTXElectron1Jet;   //!
   TBranch        *b_HLTXElectron2Jet;   //!
   TBranch        *b_HLTXElectron3Jet;   //!
   TBranch        *b_HLTXElectron4Jet;   //!
   TBranch        *b_HLTXMuonJets;   //!
   TBranch        *b_HLTXElectronMuon;   //!
   TBranch        *b_HLTXElectronMuonRelaxed;   //!
   TBranch        *b_HLTXElectronTau;   //!
   TBranch        *b_HLTXMuonTau;   //!
   TBranch        *b_CandHLTHcalIsolatedTrack;   //!
   TBranch        *b_HLTMinBiasPixel;   //!
   TBranch        *b_HLTMinBias;   //!
   TBranch        *b_HLTZeroBias;   //!
   TBranch        *b_L1_SingleMu3;   //!
   TBranch        *b_L1_SingleMu5;   //!
   TBranch        *b_L1_SingleMu7;   //!
   TBranch        *b_L1_SingleMu10;   //!
   TBranch        *b_L1_SingleMu14;   //!
   TBranch        *b_L1_SingleMu20;   //!
   TBranch        *b_L1_SingleMu25;   //!
   TBranch        *b_L1_SingleIsoEG5;   //!
   TBranch        *b_L1_SingleIsoEG8;   //!
   TBranch        *b_L1_SingleIsoEG10;   //!
   TBranch        *b_L1_SingleIsoEG12;   //!
   TBranch        *b_L1_SingleIsoEG15;   //!
   TBranch        *b_L1_SingleIsoEG20;   //!
   TBranch        *b_L1_SingleIsoEG25;   //!
   TBranch        *b_L1_SingleEG5;   //!
   TBranch        *b_L1_SingleEG8;   //!
   TBranch        *b_L1_SingleEG10;   //!
   TBranch        *b_L1_SingleEG12;   //!
   TBranch        *b_L1_SingleEG15;   //!
   TBranch        *b_L1_SingleEG20;   //!
   TBranch        *b_L1_SingleEG25;   //!
   TBranch        *b_L1_SingleJet15;   //!
   TBranch        *b_L1_SingleJet20;   //!
   TBranch        *b_L1_SingleJet30;   //!
   TBranch        *b_L1_SingleJet50;   //!
   TBranch        *b_L1_SingleJet70;   //!
   TBranch        *b_L1_SingleJet100;   //!
   TBranch        *b_L1_SingleJet150;   //!
   TBranch        *b_L1_SingleJet200;   //!
   TBranch        *b_L1_SingleTauJet10;   //!
   TBranch        *b_L1_SingleTauJet20;   //!
   TBranch        *b_L1_SingleTauJet30;   //!
   TBranch        *b_L1_SingleTauJet35;   //!
   TBranch        *b_L1_SingleTauJet40;   //!
   TBranch        *b_L1_SingleTauJet60;   //!
   TBranch        *b_L1_SingleTauJet80;   //!
   TBranch        *b_L1_SingleTauJet100;   //!
   TBranch        *b_L1_HTT100;   //!
   TBranch        *b_L1_HTT200;   //!
   TBranch        *b_L1_HTT250;   //!
   TBranch        *b_L1_HTT300;   //!
   TBranch        *b_L1_HTT400;   //!
   TBranch        *b_L1_HTT500;   //!
   TBranch        *b_L1_ETM20;   //!
   TBranch        *b_L1_ETM30;   //!
   TBranch        *b_L1_ETM40;   //!
   TBranch        *b_L1_ETM50;   //!
   TBranch        *b_L1_ETM60;   //!
   TBranch        *b_L1_DoubleMu3;   //!
   TBranch        *b_L1_DoubleIsoEG8;   //!
   TBranch        *b_L1_DoubleIsoEG10;   //!
   TBranch        *b_L1_DoubleEG5;   //!
   TBranch        *b_L1_DoubleEG10;   //!
   TBranch        *b_L1_DoubleEG15;   //!
   TBranch        *b_L1_DoubleJet70;   //!
   TBranch        *b_L1_DoubleJet100;   //!
   TBranch        *b_L1_DoubleTauJet20;   //!
   TBranch        *b_L1_DoubleTauJet30;   //!
   TBranch        *b_L1_DoubleTauJet35;   //!
   TBranch        *b_L1_DoubleTauJet40;   //!
   TBranch        *b_L1_Mu3_IsoEG5;   //!
   TBranch        *b_L1_Mu5_IsoEG10;   //!
   TBranch        *b_L1_Mu3_EG12;   //!
   TBranch        *b_L1_Mu3_Jet15;   //!
   TBranch        *b_L1_Mu5_Jet15;   //!
   TBranch        *b_L1_Mu3_Jet70;   //!
   TBranch        *b_L1_Mu5_Jet20;   //!
   TBranch        *b_L1_Mu5_TauJet20;   //!
   TBranch        *b_L1_Mu5_TauJet30;   //!
   TBranch        *b_L1_IsoEG10_EG10;   //!
   TBranch        *b_L1_IsoEG10_Jet15;   //!
   TBranch        *b_L1_IsoEG10_Jet20;   //!
   TBranch        *b_L1_IsoEG10_Jet30;   //!
   TBranch        *b_L1_IsoEG10_Jet70;   //!
   TBranch        *b_L1_IsoEG10_TauJet20;   //!
   TBranch        *b_L1_IsoEG10_TauJet30;   //!
   TBranch        *b_L1_EG10_Jet15;   //!
   TBranch        *b_L1_EG12_Jet20;   //!
   TBranch        *b_L1_EG12_Jet70;   //!
   TBranch        *b_L1_EG12_TauJet40;   //!
   TBranch        *b_L1_Jet70_TauJet40;   //!
   TBranch        *b_L1_Mu3_HTT200;   //!
   TBranch        *b_L1_IsoEG10_HTT200;   //!
   TBranch        *b_L1_EG12_HTT200;   //!
   TBranch        *b_L1_Jet70_HTT200;   //!
   TBranch        *b_L1_TauJet40_HTT200;   //!
   TBranch        *b_L1_Mu3_ETM30;   //!
   TBranch        *b_L1_IsoEG10_ETM30;   //!
   TBranch        *b_L1_EG12_ETM30;   //!
   TBranch        *b_L1_Jet70_ETM40;   //!
   TBranch        *b_L1_TauJet20_ETM20;   //!
   TBranch        *b_L1_TauJet30_ETM30;   //!
   TBranch        *b_L1_TauJet30_ETM40;   //!
   TBranch        *b_L1_HTT100_ETM30;   //!
   TBranch        *b_L1_TripleMu3;   //!
   TBranch        *b_L1_TripleIsoEG5;   //!
   TBranch        *b_L1_TripleEG10;   //!
   TBranch        *b_L1_TripleJet50;   //!
   TBranch        *b_L1_TripleTauJet40;   //!
   TBranch        *b_L1_DoubleMu3_IsoEG5;   //!
   TBranch        *b_L1_DoubleMu3_EG10;   //!
   TBranch        *b_L1_DoubleIsoEG5_Mu3;   //!
   TBranch        *b_L1_DoubleEG10_Mu3;   //!
   TBranch        *b_L1_DoubleMu3_HTT200;   //!
   TBranch        *b_L1_DoubleIsoEG5_HTT200;   //!
   TBranch        *b_L1_DoubleEG10_HTT200;   //!
   TBranch        *b_L1_DoubleJet50_HTT200;   //!
   TBranch        *b_L1_DoubleTauJet40_HTT200;   //!
   TBranch        *b_L1_DoubleMu3_ETM20;   //!
   TBranch        *b_L1_DoubleIsoEG5_ETM20;   //!
   TBranch        *b_L1_DoubleEG10_ETM20;   //!
   TBranch        *b_L1_DoubleJet50_ETM20;   //!
   TBranch        *b_L1_DoubleTauJet40_ETM20;   //!
   TBranch        *b_L1_QuadJet30;   //!
   TBranch        *b_L1_ExclusiveDoubleIsoEG6;   //!
   TBranch        *b_L1_ExclusiveDoubleJet60;   //!
   TBranch        *b_L1_ExclusiveJet25_Gap_Jet25;   //!
   TBranch        *b_L1_IsoEG10_Jet20_ForJet10;   //!
   TBranch        *b_L1_MinBias_HTT10;   //!
   TBranch        *b_L1_ZeroBias;   //!

   OHltTree(TTree *tree=0,int ntrig=0);
   virtual ~OHltTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
	 virtual void			SetStandardHLTPath();
	 virtual void			SetMapBitOfStandardHLTPath();
   virtual void     Show(Long64_t entry = -1);
   void Loop(std::vector<int> *, std::vector<int> *, std::vector<int> * 
	     ,std::vector< std::vector<int> > * overlapCount
	     ,std::vector<TString> trignames, std::map<TString,int> map_TrigPrescls
	     ,int n=-1,bool doMuonCut=false,bool doElecCut=false);


 private:
   int Ntrig;
   std::vector<int> triggerBit;
   std::vector<int> triggerBitNoPrescale;
   std::vector<int> previousBitsFired;
   std::vector<int> allOtherBitsFired;
	 std::vector<int> BitOfStandardHLTPath;
	 std::map<TString,int> map_BitOfStandardHLTPath;


};

#endif

#ifdef OHltTree_cxx
OHltTree::OHltTree(TTree *tree, int ntrig)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/tmp/test.root");
      if (!f) {
         f = new TFile("/tmp/test.root");
      }
      tree = (TTree*)gDirectory->Get("HltTree");

   }
   Init(tree);

   Ntrig = ntrig;
   triggerBit.reserve(Ntrig);
   triggerBitNoPrescale.reserve(Ntrig);
   previousBitsFired.reserve(Ntrig);
   allOtherBitsFired.reserve(Ntrig);
	 BitOfStandardHLTPath.reserve(Ntrig);

   for (int it = 0; it < Ntrig; it++){
     triggerBit.push_back(false);
     triggerBitNoPrescale.push_back(false);
     previousBitsFired.push_back(false);
     allOtherBitsFired.push_back(false);
   }
}

OHltTree::~OHltTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t OHltTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t OHltTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void OHltTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("NrecoJetCal", &NrecoJetCal, &b_NrecoJetCal);
   fChain->SetBranchAddress("NrecoJetGen", &NrecoJetGen, &b_NrecoJetGen);
   fChain->SetBranchAddress("NrecoTowCal", &NrecoTowCal, &b_NrecoTowCal);
   fChain->SetBranchAddress("recoJetCalPt", recoJetCalPt, &b_recoJetCalPt);
   fChain->SetBranchAddress("recoJetCalPhi", recoJetCalPhi, &b_recoJetCalPhi);
   fChain->SetBranchAddress("recoJetCalEta", recoJetCalEta, &b_recoJetCalEta);
   fChain->SetBranchAddress("recoJetCalEt", recoJetCalEt, &b_recoJetCalEt);
   fChain->SetBranchAddress("recoJetCalE", recoJetCalE, &b_recoJetCalE);
   fChain->SetBranchAddress("recoJetGenPt", &recoJetGenPt, &b_recoJetGenPt);
   fChain->SetBranchAddress("recoJetGenPhi", &recoJetGenPhi, &b_recoJetGenPhi);
   fChain->SetBranchAddress("recoJetGenEta", &recoJetGenEta, &b_recoJetGenEta);
   fChain->SetBranchAddress("recoJetGenEt", &recoJetGenEt, &b_recoJetGenEt);
   fChain->SetBranchAddress("recoJetGenE", &recoJetGenE, &b_recoJetGenE);
   fChain->SetBranchAddress("recoTowEt", recoTowEt, &b_recoTowEt);
   fChain->SetBranchAddress("recoTowEta", recoTowEta, &b_recoTowEta);
   fChain->SetBranchAddress("recoTowPhi", recoTowPhi, &b_recoTowPhi);
   fChain->SetBranchAddress("recoTowE", recoTowE, &b_recoTowE);
   fChain->SetBranchAddress("recoTowEm", recoTowEm, &b_recoTowEm);
   fChain->SetBranchAddress("recoTowHad", recoTowHad, &b_recoTowHad);
   fChain->SetBranchAddress("recoTowOE", recoTowOE, &b_recoTowOE);
   fChain->SetBranchAddress("recoMetCal", &recoMetCal, &b_recoMetCal);
   fChain->SetBranchAddress("recoMetCalPhi", &recoMetCalPhi, &b_recoMetCalPhi);
   fChain->SetBranchAddress("recoMetCalSum", &recoMetCalSum, &b_recoMetCalSum);
   fChain->SetBranchAddress("recoMetGen", &recoMetGen, &b_recoMetGen);
   fChain->SetBranchAddress("recoMetGenPhi", &recoMetGenPhi, &b_recoMetGenPhi);
   fChain->SetBranchAddress("recoMetGenSum", &recoMetGenSum, &b_recoMetGenSum);
   fChain->SetBranchAddress("recoHTCal", &recoHTCal, &b_recoHTCal);
   fChain->SetBranchAddress("recoHTCalPhi", &recoHTCalPhi, &b_recoHTCalPhi);
   fChain->SetBranchAddress("recoHTCalSum", &recoHTCalSum, &b_recoHTCalSum);
   fChain->SetBranchAddress("NohTau2", &NohTau2, &b_NohTau2);
   fChain->SetBranchAddress("ohTau2Eiso", ohTau2Eiso, &b_ohTau2Eiso);
   fChain->SetBranchAddress("ohTau2L2Tpt", ohTau2L2Tpt, &b_ohTau2L2Tpt);
   fChain->SetBranchAddress("ohTau2L2Tiso", ohTau2L2Tiso, &b_ohTau2L2Tiso);
   fChain->SetBranchAddress("NohTau1", &NohTau1, &b_NohTau1);
   fChain->SetBranchAddress("ohTau1Eiso", &ohTau1Eiso, &b_ohTau1Eiso);
   fChain->SetBranchAddress("ohTau1L2Tpt", &ohTau1L2Tpt, &b_ohTau1L2Tpt);
   fChain->SetBranchAddress("ohTau1L2Tiso", &ohTau1L2Tiso, &b_ohTau1L2Tiso);
   fChain->SetBranchAddress("ohTau1L3Tpt", &ohTau1L3Tpt, &b_ohTau1L3Tpt);
   fChain->SetBranchAddress("ohTau1L3Tiso", &ohTau1L3Tiso, &b_ohTau1L3Tiso);
   fChain->SetBranchAddress("NrecoElec", &NrecoElec, &b_NrecoElec);
   fChain->SetBranchAddress("recoElecPt", &recoElecPt, &b_recoElecPt);
   fChain->SetBranchAddress("recoElecPhi", &recoElecPhi, &b_recoElecPhi);
   fChain->SetBranchAddress("recoElecEta", &recoElecEta, &b_recoElecEta);
   fChain->SetBranchAddress("recoElecEt", &recoElecEt, &b_recoElecEt);
   fChain->SetBranchAddress("recoElecE", &recoElecE, &b_recoElecE);
   fChain->SetBranchAddress("NrecoPhot", &NrecoPhot, &b_NrecoPhot);
   fChain->SetBranchAddress("recoPhotPt", &recoPhotPt, &b_recoPhotPt);
   fChain->SetBranchAddress("recoPhotPhi", &recoPhotPhi, &b_recoPhotPhi);
   fChain->SetBranchAddress("recoPhotEta", &recoPhotEta, &b_recoPhotEta);
   fChain->SetBranchAddress("recoPhotEt", &recoPhotEt, &b_recoPhotEt);
   fChain->SetBranchAddress("recoPhotE", &recoPhotE, &b_recoPhotE);
   fChain->SetBranchAddress("NohPhot", &NohPhot, &b_NohPhot);
   fChain->SetBranchAddress("ohPhotEt", ohPhotEt, &b_ohPhotEt);
   fChain->SetBranchAddress("ohPhotEta", ohPhotEta, &b_ohPhotEta);
   fChain->SetBranchAddress("ohPhotPhi", ohPhotPhi, &b_ohPhotPhi);
   fChain->SetBranchAddress("ohPhotEiso", ohPhotEiso, &b_ohPhotEiso);
   fChain->SetBranchAddress("ohPhotHiso", ohPhotHiso, &b_ohPhotHiso);
   fChain->SetBranchAddress("ohPhotTiso", ohPhotTiso, &b_ohPhotTiso);
   fChain->SetBranchAddress("ohPhotL1iso", ohPhotL1iso, &b_ohPhotL1iso);
   fChain->SetBranchAddress("NohEle", &NohEle, &b_NohEle);
   fChain->SetBranchAddress("ohEleEt", &ohEleEt, &b_ohEleEt);
   fChain->SetBranchAddress("ohEleEta", &ohEleEta, &b_ohEleEta);
   fChain->SetBranchAddress("ohElePhi", &ohElePhi, &b_ohElePhi);
   fChain->SetBranchAddress("ohEleE", &ohEleE, &b_ohEleE);
   fChain->SetBranchAddress("ohEleP", &ohEleP, &b_ohEleP);
   fChain->SetBranchAddress("ohEleHiso", &ohEleHiso, &b_ohEleHiso);
   fChain->SetBranchAddress("ohEleTiso", &ohEleTiso, &b_ohEleTiso);
   fChain->SetBranchAddress("ohEleL1iso", &ohEleL1iso, &b_ohEleL1iso);
   fChain->SetBranchAddress("NrecoMuon", &NrecoMuon, &b_NrecoMuon);
   fChain->SetBranchAddress("recoMuonPt", &recoMuonPt, &b_recoMuonPt);
   fChain->SetBranchAddress("recoMuonPhi", &recoMuonPhi, &b_recoMuonPhi);
   fChain->SetBranchAddress("recoMuonEta", &recoMuonEta, &b_recoMuonEta);
   fChain->SetBranchAddress("recoMuonEt", &recoMuonEt, &b_recoMuonEt);
   fChain->SetBranchAddress("recoMuonE", &recoMuonE, &b_recoMuonE);
   fChain->SetBranchAddress("NohMuL2", &NohMuL2, &b_NohMuL2);
   fChain->SetBranchAddress("ohMuL2Pt", ohMuL2Pt, &b_ohMuL2Pt);
   fChain->SetBranchAddress("ohMuL2PtErr", ohMuL2PtErr, &b_ohMuL2PtErr);
   fChain->SetBranchAddress("ohMuL2Phi", ohMuL2Phi, &b_ohMuL2Phi);
   fChain->SetBranchAddress("ohMuL2Eta", ohMuL2Eta, &b_ohMuL2Eta);
   fChain->SetBranchAddress("ohMuL2Chg", ohMuL2Chg, &b_ohMuL2Chg);
   fChain->SetBranchAddress("ohMuL2Iso", ohMuL2Iso, &b_ohMuL2Iso);
   fChain->SetBranchAddress("ohMuL2Dr", ohMuL2Dr, &b_ohMuL2Dr);
   fChain->SetBranchAddress("ohMuL2Dz", ohMuL2Dz, &b_ohMuL2Dz);
   fChain->SetBranchAddress("NohMuL3", &NohMuL3, &b_NohMuL3);
   fChain->SetBranchAddress("ohMuL3Pt", ohMuL3Pt, &b_ohMuL3Pt);
   fChain->SetBranchAddress("ohMuL3PtErr", ohMuL3PtErr, &b_ohMuL3PtErr);
   fChain->SetBranchAddress("ohMuL3Phi", ohMuL3Phi, &b_ohMuL3Phi);
   fChain->SetBranchAddress("ohMuL3Eta", ohMuL3Eta, &b_ohMuL3Eta);
   fChain->SetBranchAddress("ohMuL3Chg", ohMuL3Chg, &b_ohMuL3Chg);
   fChain->SetBranchAddress("ohMuL3Iso", ohMuL3Iso, &b_ohMuL3Iso);
   fChain->SetBranchAddress("ohMuL3Dr", ohMuL3Dr, &b_ohMuL3Dr);
   fChain->SetBranchAddress("ohMuL3Dz", ohMuL3Dz, &b_ohMuL3Dz);
   fChain->SetBranchAddress("NMCpart", &NMCpart, &b_NMCpart);
   fChain->SetBranchAddress("MCpid", MCpid, &b_MCpid);
   fChain->SetBranchAddress("MCvtxX", MCvtxX, &b_MCvtxX);
   fChain->SetBranchAddress("MCvtxY", MCvtxY, &b_MCvtxY);
   fChain->SetBranchAddress("MCvtxZ", MCvtxZ, &b_MCvtxZ);
   fChain->SetBranchAddress("MCpt", MCpt, &b_MCpt);
   fChain->SetBranchAddress("MCeta", MCeta, &b_MCeta);
   fChain->SetBranchAddress("MCphi", MCphi, &b_MCphi);
   fChain->SetBranchAddress("MCPtHat", &MCPtHat, &b_MCPtHat);
   fChain->SetBranchAddress("MCmu3", &MCmu3, &b_MCmu3);
   fChain->SetBranchAddress("MCel1", &MCel1, &b_MCel1);
   fChain->SetBranchAddress("MCbb", &MCbb, &b_MCbb);
   fChain->SetBranchAddress("MCab", &MCab, &b_MCab);
   fChain->SetBranchAddress("NL1IsolEm", &NL1IsolEm, &b_NL1IsolEm);
   fChain->SetBranchAddress("L1IsolEmEt", L1IsolEmEt, &b_L1IsolEmEt);
   fChain->SetBranchAddress("L1IsolEmE", L1IsolEmE, &b_L1IsolEmE);
   fChain->SetBranchAddress("L1IsolEmEta", L1IsolEmEta, &b_L1IsolEmEta);
   fChain->SetBranchAddress("L1IsolEmPhi", L1IsolEmPhi, &b_L1IsolEmPhi);
   fChain->SetBranchAddress("NL1NIsolEm", &NL1NIsolEm, &b_NL1NIsolEm);
   fChain->SetBranchAddress("L1NIsolEmEt", L1NIsolEmEt, &b_L1NIsolEmEt);
   fChain->SetBranchAddress("L1NIsolEmE", L1NIsolEmE, &b_L1NIsolEmE);
   fChain->SetBranchAddress("L1NIsolEmEta", L1NIsolEmEta, &b_L1NIsolEmEta);
   fChain->SetBranchAddress("L1NIsolEmPhi", L1NIsolEmPhi, &b_L1NIsolEmPhi);
   fChain->SetBranchAddress("NL1Mu", &NL1Mu, &b_NL1Mu);
   fChain->SetBranchAddress("L1MuPt", L1MuPt, &b_L1MuPt);
   fChain->SetBranchAddress("L1MuE", L1MuE, &b_L1MuE);
   fChain->SetBranchAddress("L1MuEta", L1MuEta, &b_L1MuEta);
   fChain->SetBranchAddress("L1MuPhi", L1MuPhi, &b_L1MuPhi);
   fChain->SetBranchAddress("L1MuIsol", L1MuIsol, &b_L1MuIsol);
   fChain->SetBranchAddress("L1MuMip", L1MuMip, &b_L1MuMip);
   fChain->SetBranchAddress("L1MuFor", L1MuFor, &b_L1MuFor);
   fChain->SetBranchAddress("L1MuRPC", L1MuRPC, &b_L1MuRPC);
   fChain->SetBranchAddress("L1MuQal", L1MuQal, &b_L1MuQal);
   fChain->SetBranchAddress("NL1CenJet", &NL1CenJet, &b_NL1CenJet);
   fChain->SetBranchAddress("L1CenJetEt", L1CenJetEt, &b_L1CenJetEt);
   fChain->SetBranchAddress("L1CenJetE", L1CenJetE, &b_L1CenJetE);
   fChain->SetBranchAddress("L1CenJetEta", L1CenJetEta, &b_L1CenJetEta);
   fChain->SetBranchAddress("L1CenJetPhi", L1CenJetPhi, &b_L1CenJetPhi);
   fChain->SetBranchAddress("NL1ForJet", &NL1ForJet, &b_NL1ForJet);
   fChain->SetBranchAddress("L1ForJetEt", L1ForJetEt, &b_L1ForJetEt);
   fChain->SetBranchAddress("L1ForJetE", L1ForJetE, &b_L1ForJetE);
   fChain->SetBranchAddress("L1ForJetEta", L1ForJetEta, &b_L1ForJetEta);
   fChain->SetBranchAddress("L1ForJetPhi", L1ForJetPhi, &b_L1ForJetPhi);
   fChain->SetBranchAddress("NL1Tau", &NL1Tau, &b_NL1Tau);
   fChain->SetBranchAddress("L1TauEt", L1TauEt, &b_L1TauEt);
   fChain->SetBranchAddress("L1TauE", L1TauE, &b_L1TauE);
   fChain->SetBranchAddress("L1TauEta", L1TauEta, &b_L1TauEta);
   fChain->SetBranchAddress("L1TauPhi", L1TauPhi, &b_L1TauPhi);
   fChain->SetBranchAddress("L1Met", &L1Met, &b_L1Met);
   fChain->SetBranchAddress("L1MetPhi", &L1MetPhi, &b_L1MetPhi);
   fChain->SetBranchAddress("L1MetTot", &L1MetTot, &b_L1MetTot);
   fChain->SetBranchAddress("L1MetHad", &L1MetHad, &b_L1MetHad);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("HLT1jet", &HLT1jet, &b_HLT1jet);
   fChain->SetBranchAddress("HLT2jet", &HLT2jet, &b_HLT2jet);
   fChain->SetBranchAddress("HLT3jet", &HLT3jet, &b_HLT3jet);
   fChain->SetBranchAddress("HLT4jet", &HLT4jet, &b_HLT4jet);
   fChain->SetBranchAddress("HLT1MET", &HLT1MET, &b_HLT1MET);
   fChain->SetBranchAddress("HLT2jetAco", &HLT2jetAco, &b_HLT2jetAco);
   fChain->SetBranchAddress("HLT1jet1METAco", &HLT1jet1METAco, &b_HLT1jet1METAco);
   fChain->SetBranchAddress("HLT1jet1MET", &HLT1jet1MET, &b_HLT1jet1MET);
   fChain->SetBranchAddress("HLT2jet1MET", &HLT2jet1MET, &b_HLT2jet1MET);
   fChain->SetBranchAddress("HLT3jet1MET", &HLT3jet1MET, &b_HLT3jet1MET);
   fChain->SetBranchAddress("HLT4jet1MET", &HLT4jet1MET, &b_HLT4jet1MET);
   fChain->SetBranchAddress("HLT1MET1HT", &HLT1MET1HT, &b_HLT1MET1HT);
   fChain->SetBranchAddress("CandHLT1SumET", &CandHLT1SumET, &b_CandHLT1SumET);
   fChain->SetBranchAddress("HLT1jetPE1", &HLT1jetPE1, &b_HLT1jetPE1);
   fChain->SetBranchAddress("HLT1jetPE3", &HLT1jetPE3, &b_HLT1jetPE3);
   fChain->SetBranchAddress("HLT1jetPE5", &HLT1jetPE5, &b_HLT1jetPE5);
   fChain->SetBranchAddress("CandHLT1jetPE7", &CandHLT1jetPE7, &b_CandHLT1jetPE7);
   fChain->SetBranchAddress("CandHLT1METPre1", &CandHLT1METPre1, &b_CandHLT1METPre1);
   fChain->SetBranchAddress("CandHLT1METPre2", &CandHLT1METPre2, &b_CandHLT1METPre2);
   fChain->SetBranchAddress("CandHLT1METPre3", &CandHLT1METPre3, &b_CandHLT1METPre3);
   fChain->SetBranchAddress("CandHLT2jetAve30", &CandHLT2jetAve30, &b_CandHLT2jetAve30);
   fChain->SetBranchAddress("CandHLT2jetAve60", &CandHLT2jetAve60, &b_CandHLT2jetAve60);
   fChain->SetBranchAddress("CandHLT2jetAve110", &CandHLT2jetAve110, &b_CandHLT2jetAve110);
   fChain->SetBranchAddress("CandHLT2jetAve150", &CandHLT2jetAve150, &b_CandHLT2jetAve150);
   fChain->SetBranchAddress("CandHLT2jetAve200", &CandHLT2jetAve200, &b_CandHLT2jetAve200);
   fChain->SetBranchAddress("HLT2jetvbfMET", &HLT2jetvbfMET, &b_HLT2jetvbfMET);
   fChain->SetBranchAddress("HLTS2jet1METNV", &HLTS2jet1METNV, &b_HLTS2jet1METNV);
   fChain->SetBranchAddress("HLTS2jet1METAco", &HLTS2jet1METAco, &b_HLTS2jet1METAco);
   fChain->SetBranchAddress("CandHLTSjet1MET1Aco", &CandHLTSjet1MET1Aco, &b_CandHLTSjet1MET1Aco);
   fChain->SetBranchAddress("CandHLTSjet2MET1Aco", &CandHLTSjet2MET1Aco, &b_CandHLTSjet2MET1Aco);
   fChain->SetBranchAddress("CandHLTS2jetAco", &CandHLTS2jetAco, &b_CandHLTS2jetAco);
   fChain->SetBranchAddress("CandHLTJetMETRapidityGap", &CandHLTJetMETRapidityGap, &b_CandHLTJetMETRapidityGap);
   fChain->SetBranchAddress("HLT1Electron", &HLT1Electron, &b_HLT1Electron);
   fChain->SetBranchAddress("HLT1ElectronRelaxed", &HLT1ElectronRelaxed, &b_HLT1ElectronRelaxed);
   fChain->SetBranchAddress("HLT2Electron", &HLT2Electron, &b_HLT2Electron);
   fChain->SetBranchAddress("HLT2ElectronRelaxed", &HLT2ElectronRelaxed, &b_HLT2ElectronRelaxed);
   fChain->SetBranchAddress("HLT1Photon", &HLT1Photon, &b_HLT1Photon);
   fChain->SetBranchAddress("HLT1PhotonRelaxed", &HLT1PhotonRelaxed, &b_HLT1PhotonRelaxed);
   fChain->SetBranchAddress("HLT2Photon", &HLT2Photon, &b_HLT2Photon);
   fChain->SetBranchAddress("HLT2PhotonRelaxed", &HLT2PhotonRelaxed, &b_HLT2PhotonRelaxed);
   fChain->SetBranchAddress("HLT1EMHighEt", &HLT1EMHighEt, &b_HLT1EMHighEt);
   fChain->SetBranchAddress("HLT1EMVeryHighEt", &HLT1EMVeryHighEt, &b_HLT1EMVeryHighEt);
   fChain->SetBranchAddress("CandHLT2ElectronZCounter", &CandHLT2ElectronZCounter, &b_CandHLT2ElectronZCounter);
   fChain->SetBranchAddress("CandHLT2ElectronExclusive", &CandHLT2ElectronExclusive, &b_CandHLT2ElectronExclusive);
   fChain->SetBranchAddress("CandHLT2PhotonExclusive", &CandHLT2PhotonExclusive, &b_CandHLT2PhotonExclusive);
   fChain->SetBranchAddress("CandHLT1PhotonL1Isolated", &CandHLT1PhotonL1Isolated, &b_CandHLT1PhotonL1Isolated);
   fChain->SetBranchAddress("HLT1MuonIso", &HLT1MuonIso, &b_HLT1MuonIso);
   fChain->SetBranchAddress("HLT1MuonNonIso", &HLT1MuonNonIso, &b_HLT1MuonNonIso);
   fChain->SetBranchAddress("CandHLT2MuonIso", &CandHLT2MuonIso, &b_CandHLT2MuonIso);
   fChain->SetBranchAddress("HLT2MuonNonIso", &HLT2MuonNonIso, &b_HLT2MuonNonIso);
   fChain->SetBranchAddress("HLT2MuonJPsi", &HLT2MuonJPsi, &b_HLT2MuonJPsi);
   fChain->SetBranchAddress("HLT2MuonUpsilon", &HLT2MuonUpsilon, &b_HLT2MuonUpsilon);
   fChain->SetBranchAddress("HLT2MuonZ", &HLT2MuonZ, &b_HLT2MuonZ);
   fChain->SetBranchAddress("HLTNMuonNonIso", &HLTNMuonNonIso, &b_HLTNMuonNonIso);
   fChain->SetBranchAddress("HLT2MuonSameSign", &HLT2MuonSameSign, &b_HLT2MuonSameSign);
   fChain->SetBranchAddress("CandHLT1MuonPrescalePt3", &CandHLT1MuonPrescalePt3, &b_CandHLT1MuonPrescalePt3);
   fChain->SetBranchAddress("CandHLT1MuonPrescalePt5", &CandHLT1MuonPrescalePt5, &b_CandHLT1MuonPrescalePt5);
   fChain->SetBranchAddress("CandHLT1MuonPrescalePt7x7", &CandHLT1MuonPrescalePt7x7, &b_CandHLT1MuonPrescalePt7x7);
   fChain->SetBranchAddress("CandHLT1MuonPrescalePt7x10", &CandHLT1MuonPrescalePt7x10, &b_CandHLT1MuonPrescalePt7x10);
   fChain->SetBranchAddress("CandHLT1MuonLevel1", &CandHLT1MuonLevel1, &b_CandHLT1MuonLevel1);
   fChain->SetBranchAddress("HLTB1Jet", &HLTB1Jet, &b_HLTB1Jet);
   fChain->SetBranchAddress("HLTB2Jet", &HLTB2Jet, &b_HLTB2Jet);
   fChain->SetBranchAddress("HLTB3Jet", &HLTB3Jet, &b_HLTB3Jet);
   fChain->SetBranchAddress("HLTB4Jet", &HLTB4Jet, &b_HLTB4Jet);
   fChain->SetBranchAddress("HLTBHT", &HLTBHT, &b_HLTBHT);
   fChain->SetBranchAddress("HLTB1JetMu", &HLTB1JetMu, &b_HLTB1JetMu);
   fChain->SetBranchAddress("HLTB2JetMu", &HLTB2JetMu, &b_HLTB2JetMu);
   fChain->SetBranchAddress("HLTB3JetMu", &HLTB3JetMu, &b_HLTB3JetMu);
   fChain->SetBranchAddress("HLTB4JetMu", &HLTB4JetMu, &b_HLTB4JetMu);
   fChain->SetBranchAddress("HLTBHTMu", &HLTBHTMu, &b_HLTBHTMu);
   fChain->SetBranchAddress("HLTBJPsiMuMu", &HLTBJPsiMuMu, &b_HLTBJPsiMuMu);
   fChain->SetBranchAddress("HLT1Tau", &HLT1Tau, &b_HLT1Tau);
   fChain->SetBranchAddress("HLT1Tau1MET", &HLT1Tau1MET, &b_HLT1Tau1MET);
   fChain->SetBranchAddress("HLT2TauPixel", &HLT2TauPixel, &b_HLT2TauPixel);
   fChain->SetBranchAddress("HLTXElectronBJet", &HLTXElectronBJet, &b_HLTXElectronBJet);
   fChain->SetBranchAddress("HLTXMuonBJet", &HLTXMuonBJet, &b_HLTXMuonBJet);
   fChain->SetBranchAddress("HLTXMuonBJetSoftMuon", &HLTXMuonBJetSoftMuon, &b_HLTXMuonBJetSoftMuon);
   fChain->SetBranchAddress("HLTXElectron1Jet", &HLTXElectron1Jet, &b_HLTXElectron1Jet);
   fChain->SetBranchAddress("HLTXElectron2Jet", &HLTXElectron2Jet, &b_HLTXElectron2Jet);
   fChain->SetBranchAddress("HLTXElectron3Jet", &HLTXElectron3Jet, &b_HLTXElectron3Jet);
   fChain->SetBranchAddress("HLTXElectron4Jet", &HLTXElectron4Jet, &b_HLTXElectron4Jet);
   fChain->SetBranchAddress("HLTXMuonJets", &HLTXMuonJets, &b_HLTXMuonJets);
   fChain->SetBranchAddress("HLTXElectronMuon", &HLTXElectronMuon, &b_HLTXElectronMuon);
   fChain->SetBranchAddress("HLTXElectronMuonRelaxed", &HLTXElectronMuonRelaxed, &b_HLTXElectronMuonRelaxed);
   fChain->SetBranchAddress("HLTXElectronTau", &HLTXElectronTau, &b_HLTXElectronTau);
   fChain->SetBranchAddress("HLTXMuonTau", &HLTXMuonTau, &b_HLTXMuonTau);
   fChain->SetBranchAddress("CandHLTHcalIsolatedTrack", &CandHLTHcalIsolatedTrack, &b_CandHLTHcalIsolatedTrack);
   fChain->SetBranchAddress("HLTMinBiasPixel", &HLTMinBiasPixel, &b_HLTMinBiasPixel);
   fChain->SetBranchAddress("HLTMinBias", &HLTMinBias, &b_HLTMinBias);
   fChain->SetBranchAddress("HLTZeroBias", &HLTZeroBias, &b_HLTZeroBias);
   fChain->SetBranchAddress("L1_SingleMu3", &L1_SingleMu3, &b_L1_SingleMu3);
   fChain->SetBranchAddress("L1_SingleMu5", &L1_SingleMu5, &b_L1_SingleMu5);
   fChain->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7, &b_L1_SingleMu7);
   fChain->SetBranchAddress("L1_SingleMu10", &L1_SingleMu10, &b_L1_SingleMu10);
   fChain->SetBranchAddress("L1_SingleMu14", &L1_SingleMu14, &b_L1_SingleMu14);
   fChain->SetBranchAddress("L1_SingleMu20", &L1_SingleMu20, &b_L1_SingleMu20);
   fChain->SetBranchAddress("L1_SingleMu25", &L1_SingleMu25, &b_L1_SingleMu25);
   fChain->SetBranchAddress("L1_SingleIsoEG5", &L1_SingleIsoEG5, &b_L1_SingleIsoEG5);
   fChain->SetBranchAddress("L1_SingleIsoEG8", &L1_SingleIsoEG8, &b_L1_SingleIsoEG8);
   fChain->SetBranchAddress("L1_SingleIsoEG10", &L1_SingleIsoEG10, &b_L1_SingleIsoEG10);
   fChain->SetBranchAddress("L1_SingleIsoEG12", &L1_SingleIsoEG12, &b_L1_SingleIsoEG12);
   fChain->SetBranchAddress("L1_SingleIsoEG15", &L1_SingleIsoEG15, &b_L1_SingleIsoEG15);
   fChain->SetBranchAddress("L1_SingleIsoEG20", &L1_SingleIsoEG20, &b_L1_SingleIsoEG20);
   fChain->SetBranchAddress("L1_SingleIsoEG25", &L1_SingleIsoEG25, &b_L1_SingleIsoEG25);
   fChain->SetBranchAddress("L1_SingleEG5", &L1_SingleEG5, &b_L1_SingleEG5);
   fChain->SetBranchAddress("L1_SingleEG8", &L1_SingleEG8, &b_L1_SingleEG8);
   fChain->SetBranchAddress("L1_SingleEG10", &L1_SingleEG10, &b_L1_SingleEG10);
   fChain->SetBranchAddress("L1_SingleEG12", &L1_SingleEG12, &b_L1_SingleEG12);
   fChain->SetBranchAddress("L1_SingleEG15", &L1_SingleEG15, &b_L1_SingleEG15);
   fChain->SetBranchAddress("L1_SingleEG20", &L1_SingleEG20, &b_L1_SingleEG20);
   fChain->SetBranchAddress("L1_SingleEG25", &L1_SingleEG25, &b_L1_SingleEG25);
   fChain->SetBranchAddress("L1_SingleJet15", &L1_SingleJet15, &b_L1_SingleJet15);
   fChain->SetBranchAddress("L1_SingleJet20", &L1_SingleJet20, &b_L1_SingleJet20);
   fChain->SetBranchAddress("L1_SingleJet30", &L1_SingleJet30, &b_L1_SingleJet30);
   fChain->SetBranchAddress("L1_SingleJet50", &L1_SingleJet50, &b_L1_SingleJet50);
   fChain->SetBranchAddress("L1_SingleJet70", &L1_SingleJet70, &b_L1_SingleJet70);
   fChain->SetBranchAddress("L1_SingleJet100", &L1_SingleJet100, &b_L1_SingleJet100);
   fChain->SetBranchAddress("L1_SingleJet150", &L1_SingleJet150, &b_L1_SingleJet150);
   fChain->SetBranchAddress("L1_SingleJet200", &L1_SingleJet200, &b_L1_SingleJet200);
   fChain->SetBranchAddress("L1_SingleTauJet10", &L1_SingleTauJet10, &b_L1_SingleTauJet10);
   fChain->SetBranchAddress("L1_SingleTauJet20", &L1_SingleTauJet20, &b_L1_SingleTauJet20);
   fChain->SetBranchAddress("L1_SingleTauJet30", &L1_SingleTauJet30, &b_L1_SingleTauJet30);
   fChain->SetBranchAddress("L1_SingleTauJet35", &L1_SingleTauJet35, &b_L1_SingleTauJet35);
   fChain->SetBranchAddress("L1_SingleTauJet40", &L1_SingleTauJet40, &b_L1_SingleTauJet40);
   fChain->SetBranchAddress("L1_SingleTauJet60", &L1_SingleTauJet60, &b_L1_SingleTauJet60);
   fChain->SetBranchAddress("L1_SingleTauJet80", &L1_SingleTauJet80, &b_L1_SingleTauJet80);
   fChain->SetBranchAddress("L1_SingleTauJet100", &L1_SingleTauJet100, &b_L1_SingleTauJet100);
   fChain->SetBranchAddress("L1_HTT100", &L1_HTT100, &b_L1_HTT100);
   fChain->SetBranchAddress("L1_HTT200", &L1_HTT200, &b_L1_HTT200);
   fChain->SetBranchAddress("L1_HTT250", &L1_HTT250, &b_L1_HTT250);
   fChain->SetBranchAddress("L1_HTT300", &L1_HTT300, &b_L1_HTT300);
   fChain->SetBranchAddress("L1_HTT400", &L1_HTT400, &b_L1_HTT400);
   fChain->SetBranchAddress("L1_HTT500", &L1_HTT500, &b_L1_HTT500);
   fChain->SetBranchAddress("L1_ETM20", &L1_ETM20, &b_L1_ETM20);
   fChain->SetBranchAddress("L1_ETM30", &L1_ETM30, &b_L1_ETM30);
   fChain->SetBranchAddress("L1_ETM40", &L1_ETM40, &b_L1_ETM40);
   fChain->SetBranchAddress("L1_ETM50", &L1_ETM50, &b_L1_ETM50);
   fChain->SetBranchAddress("L1_ETM60", &L1_ETM60, &b_L1_ETM60);
   fChain->SetBranchAddress("L1_DoubleMu3", &L1_DoubleMu3, &b_L1_DoubleMu3);
   fChain->SetBranchAddress("L1_DoubleIsoEG8", &L1_DoubleIsoEG8, &b_L1_DoubleIsoEG8);
   fChain->SetBranchAddress("L1_DoubleIsoEG10", &L1_DoubleIsoEG10, &b_L1_DoubleIsoEG10);
   fChain->SetBranchAddress("L1_DoubleEG5", &L1_DoubleEG5, &b_L1_DoubleEG5);
   fChain->SetBranchAddress("L1_DoubleEG10", &L1_DoubleEG10, &b_L1_DoubleEG10);
   fChain->SetBranchAddress("L1_DoubleEG15", &L1_DoubleEG15, &b_L1_DoubleEG15);
   fChain->SetBranchAddress("L1_DoubleJet70", &L1_DoubleJet70, &b_L1_DoubleJet70);
   fChain->SetBranchAddress("L1_DoubleJet100", &L1_DoubleJet100, &b_L1_DoubleJet100);
   fChain->SetBranchAddress("L1_DoubleTauJet20", &L1_DoubleTauJet20, &b_L1_DoubleTauJet20);
   fChain->SetBranchAddress("L1_DoubleTauJet30", &L1_DoubleTauJet30, &b_L1_DoubleTauJet30);
   fChain->SetBranchAddress("L1_DoubleTauJet35", &L1_DoubleTauJet35, &b_L1_DoubleTauJet35);
   fChain->SetBranchAddress("L1_DoubleTauJet40", &L1_DoubleTauJet40, &b_L1_DoubleTauJet40);
   fChain->SetBranchAddress("L1_Mu3_IsoEG5", &L1_Mu3_IsoEG5, &b_L1_Mu3_IsoEG5);
   fChain->SetBranchAddress("L1_Mu5_IsoEG10", &L1_Mu5_IsoEG10, &b_L1_Mu5_IsoEG10);
   fChain->SetBranchAddress("L1_Mu3_EG12", &L1_Mu3_EG12, &b_L1_Mu3_EG12);
   fChain->SetBranchAddress("L1_Mu3_Jet15", &L1_Mu3_Jet15, &b_L1_Mu3_Jet15);
   fChain->SetBranchAddress("L1_Mu5_Jet15", &L1_Mu5_Jet15, &b_L1_Mu5_Jet15);
   fChain->SetBranchAddress("L1_Mu3_Jet70", &L1_Mu3_Jet70, &b_L1_Mu3_Jet70);
   fChain->SetBranchAddress("L1_Mu5_Jet20", &L1_Mu5_Jet20, &b_L1_Mu5_Jet20);
   fChain->SetBranchAddress("L1_Mu5_TauJet20", &L1_Mu5_TauJet20, &b_L1_Mu5_TauJet20);
   fChain->SetBranchAddress("L1_Mu5_TauJet30", &L1_Mu5_TauJet30, &b_L1_Mu5_TauJet30);
   fChain->SetBranchAddress("L1_IsoEG10_EG10", &L1_IsoEG10_EG10, &b_L1_IsoEG10_EG10);
   fChain->SetBranchAddress("L1_IsoEG10_Jet15", &L1_IsoEG10_Jet15, &b_L1_IsoEG10_Jet15);
   fChain->SetBranchAddress("L1_IsoEG10_Jet20", &L1_IsoEG10_Jet20, &b_L1_IsoEG10_Jet20);
   fChain->SetBranchAddress("L1_IsoEG10_Jet30", &L1_IsoEG10_Jet30, &b_L1_IsoEG10_Jet30);
   fChain->SetBranchAddress("L1_IsoEG10_Jet70", &L1_IsoEG10_Jet70, &b_L1_IsoEG10_Jet70);
   fChain->SetBranchAddress("L1_IsoEG10_TauJet20", &L1_IsoEG10_TauJet20, &b_L1_IsoEG10_TauJet20);
   fChain->SetBranchAddress("L1_IsoEG10_TauJet30", &L1_IsoEG10_TauJet30, &b_L1_IsoEG10_TauJet30);
   fChain->SetBranchAddress("L1_EG10_Jet15", &L1_EG10_Jet15, &b_L1_EG10_Jet15);
   fChain->SetBranchAddress("L1_EG12_Jet20", &L1_EG12_Jet20, &b_L1_EG12_Jet20);
   fChain->SetBranchAddress("L1_EG12_Jet70", &L1_EG12_Jet70, &b_L1_EG12_Jet70);
   fChain->SetBranchAddress("L1_EG12_TauJet40", &L1_EG12_TauJet40, &b_L1_EG12_TauJet40);
   fChain->SetBranchAddress("L1_Jet70_TauJet40", &L1_Jet70_TauJet40, &b_L1_Jet70_TauJet40);
   fChain->SetBranchAddress("L1_Mu3_HTT200", &L1_Mu3_HTT200, &b_L1_Mu3_HTT200);
   fChain->SetBranchAddress("L1_IsoEG10_HTT200", &L1_IsoEG10_HTT200, &b_L1_IsoEG10_HTT200);
   fChain->SetBranchAddress("L1_EG12_HTT200", &L1_EG12_HTT200, &b_L1_EG12_HTT200);
   fChain->SetBranchAddress("L1_Jet70_HTT200", &L1_Jet70_HTT200, &b_L1_Jet70_HTT200);
   fChain->SetBranchAddress("L1_TauJet40_HTT200", &L1_TauJet40_HTT200, &b_L1_TauJet40_HTT200);
   fChain->SetBranchAddress("L1_Mu3_ETM30", &L1_Mu3_ETM30, &b_L1_Mu3_ETM30);
   fChain->SetBranchAddress("L1_IsoEG10_ETM30", &L1_IsoEG10_ETM30, &b_L1_IsoEG10_ETM30);
   fChain->SetBranchAddress("L1_EG12_ETM30", &L1_EG12_ETM30, &b_L1_EG12_ETM30);
   fChain->SetBranchAddress("L1_Jet70_ETM40", &L1_Jet70_ETM40, &b_L1_Jet70_ETM40);
   fChain->SetBranchAddress("L1_TauJet20_ETM20", &L1_TauJet20_ETM20, &b_L1_TauJet20_ETM20);
   fChain->SetBranchAddress("L1_TauJet30_ETM30", &L1_TauJet30_ETM30, &b_L1_TauJet30_ETM30);
   fChain->SetBranchAddress("L1_TauJet30_ETM40", &L1_TauJet30_ETM40, &b_L1_TauJet30_ETM40);
   fChain->SetBranchAddress("L1_HTT100_ETM30", &L1_HTT100_ETM30, &b_L1_HTT100_ETM30);
   fChain->SetBranchAddress("L1_TripleMu3", &L1_TripleMu3, &b_L1_TripleMu3);
   fChain->SetBranchAddress("L1_TripleIsoEG5", &L1_TripleIsoEG5, &b_L1_TripleIsoEG5);
   fChain->SetBranchAddress("L1_TripleEG10", &L1_TripleEG10, &b_L1_TripleEG10);
   fChain->SetBranchAddress("L1_TripleJet50", &L1_TripleJet50, &b_L1_TripleJet50);
   fChain->SetBranchAddress("L1_TripleTauJet40", &L1_TripleTauJet40, &b_L1_TripleTauJet40);
   fChain->SetBranchAddress("L1_DoubleMu3_IsoEG5", &L1_DoubleMu3_IsoEG5, &b_L1_DoubleMu3_IsoEG5);
   fChain->SetBranchAddress("L1_DoubleMu3_EG10", &L1_DoubleMu3_EG10, &b_L1_DoubleMu3_EG10);
   fChain->SetBranchAddress("L1_DoubleIsoEG5_Mu3", &L1_DoubleIsoEG5_Mu3, &b_L1_DoubleIsoEG5_Mu3);
   fChain->SetBranchAddress("L1_DoubleEG10_Mu3", &L1_DoubleEG10_Mu3, &b_L1_DoubleEG10_Mu3);
   fChain->SetBranchAddress("L1_DoubleMu3_HTT200", &L1_DoubleMu3_HTT200, &b_L1_DoubleMu3_HTT200);
   fChain->SetBranchAddress("L1_DoubleIsoEG5_HTT200", &L1_DoubleIsoEG5_HTT200, &b_L1_DoubleIsoEG5_HTT200);
   fChain->SetBranchAddress("L1_DoubleEG10_HTT200", &L1_DoubleEG10_HTT200, &b_L1_DoubleEG10_HTT200);
   fChain->SetBranchAddress("L1_DoubleJet50_HTT200", &L1_DoubleJet50_HTT200, &b_L1_DoubleJet50_HTT200);
   fChain->SetBranchAddress("L1_DoubleTauJet40_HTT200", &L1_DoubleTauJet40_HTT200, &b_L1_DoubleTauJet40_HTT200);
   fChain->SetBranchAddress("L1_DoubleMu3_ETM20", &L1_DoubleMu3_ETM20, &b_L1_DoubleMu3_ETM20);
   fChain->SetBranchAddress("L1_DoubleIsoEG5_ETM20", &L1_DoubleIsoEG5_ETM20, &b_L1_DoubleIsoEG5_ETM20);
   fChain->SetBranchAddress("L1_DoubleEG10_ETM20", &L1_DoubleEG10_ETM20, &b_L1_DoubleEG10_ETM20);
   fChain->SetBranchAddress("L1_DoubleJet50_ETM20", &L1_DoubleJet50_ETM20, &b_L1_DoubleJet50_ETM20);
   fChain->SetBranchAddress("L1_DoubleTauJet40_ETM20", &L1_DoubleTauJet40_ETM20, &b_L1_DoubleTauJet40_ETM20);
   fChain->SetBranchAddress("L1_QuadJet30", &L1_QuadJet30, &b_L1_QuadJet30);
   fChain->SetBranchAddress("L1_ExclusiveDoubleIsoEG6", &L1_ExclusiveDoubleIsoEG6, &b_L1_ExclusiveDoubleIsoEG6);
   fChain->SetBranchAddress("L1_ExclusiveDoubleJet60", &L1_ExclusiveDoubleJet60, &b_L1_ExclusiveDoubleJet60);
   fChain->SetBranchAddress("L1_ExclusiveJet25_Gap_Jet25", &L1_ExclusiveJet25_Gap_Jet25, &b_L1_ExclusiveJet25_Gap_Jet25);
   fChain->SetBranchAddress("L1_IsoEG10_Jet20_ForJet10", &L1_IsoEG10_Jet20_ForJet10, &b_L1_IsoEG10_Jet20_ForJet10);
   fChain->SetBranchAddress("L1_MinBias_HTT10", &L1_MinBias_HTT10, &b_L1_MinBias_HTT10);
   fChain->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias, &b_L1_ZeroBias);
	
 

   Notify();
}

void OHltTree::SetMapBitOfStandardHLTPath() {

  map_BitOfStandardHLTPath["HLT1jet"] =  HLT1jet;
  map_BitOfStandardHLTPath["HLT2jet"] =  HLT2jet;
  map_BitOfStandardHLTPath["HLT3jet"] =  HLT3jet;
  map_BitOfStandardHLTPath["HLT4jet"] =  HLT4jet;
  map_BitOfStandardHLTPath["HLT1MET"] =  HLT1MET;
  map_BitOfStandardHLTPath["HLT2jetAco"] =  HLT2jetAco;
  map_BitOfStandardHLTPath["HLT1jet1METAco"] =  HLT1jet1METAco;
  map_BitOfStandardHLTPath["HLT1jet1MET"] =  HLT1jet1MET;
  map_BitOfStandardHLTPath["HLT2jet1MET"] =  HLT2jet1MET;
  map_BitOfStandardHLTPath["HLT3jet1MET"] =  HLT3jet1MET;
  map_BitOfStandardHLTPath["HLT4jet1MET"] =  HLT4jet1MET;
  map_BitOfStandardHLTPath["HLT1MET1HT"] =  HLT1MET1HT;
  map_BitOfStandardHLTPath["HLT1SumET"] =  CandHLT1SumET;
  map_BitOfStandardHLTPath["HLT1jetPE1"] =  HLT1jetPE1;
  map_BitOfStandardHLTPath["HLT1jetPE3"] =  HLT1jetPE3;
  map_BitOfStandardHLTPath["HLT1jetPE5"] =  HLT1jetPE5;
  map_BitOfStandardHLTPath["HLT1jetPE7"] =  CandHLT1jetPE7;
  map_BitOfStandardHLTPath["HLT1METPre1"] =  CandHLT1METPre1;
  map_BitOfStandardHLTPath["HLT1METPre2"] =  CandHLT1METPre2;
  map_BitOfStandardHLTPath["HLT1METPre3"] =  CandHLT1METPre3;
  map_BitOfStandardHLTPath["HLT2jetAve30"] =  CandHLT2jetAve30;
  map_BitOfStandardHLTPath["HLT2jetAve60"] =  CandHLT2jetAve60;
  map_BitOfStandardHLTPath["HLT2jetAve110"] =  CandHLT2jetAve110;
  map_BitOfStandardHLTPath["HLT2jetAve150"] =  CandHLT2jetAve150;
  map_BitOfStandardHLTPath["HLT2jetAve200"] =  CandHLT2jetAve200;
  map_BitOfStandardHLTPath["HLT2jetvbfMET"] =  HLT2jetvbfMET;
  map_BitOfStandardHLTPath["HLTS2jet1METNV"] =  HLTS2jet1METNV;
  map_BitOfStandardHLTPath["HLTS2jet1METAco"] =  HLTS2jet1METAco;
  map_BitOfStandardHLTPath["HLTSjet1MET1Aco"] =  CandHLTSjet1MET1Aco;
  map_BitOfStandardHLTPath["HLTSjet2MET1Aco"] =  CandHLTSjet2MET1Aco;
  map_BitOfStandardHLTPath["HLTS2jetAco"] =  CandHLTS2jetAco;
  map_BitOfStandardHLTPath["HLTJetMETRapidityGap"] =  CandHLTJetMETRapidityGap;
  map_BitOfStandardHLTPath["HLT1Electron"] =  HLT1Electron;
  map_BitOfStandardHLTPath["HLT1ElectronRelaxed"] =  HLT1ElectronRelaxed;
  map_BitOfStandardHLTPath["HLT2Electron"] =  HLT2Electron;
  map_BitOfStandardHLTPath["HLT2ElectronRelaxed"] =  HLT2ElectronRelaxed;
  map_BitOfStandardHLTPath["HLT1Photon"] =  HLT1Photon;
  map_BitOfStandardHLTPath["HLT1PhotonRelaxed"] =  HLT1PhotonRelaxed;
  map_BitOfStandardHLTPath["HLT2Photon"] =  HLT2Photon;
  map_BitOfStandardHLTPath["HLT2PhotonRelaxed"] =  HLT2PhotonRelaxed;
  map_BitOfStandardHLTPath["HLT1EMHighEt"] =  HLT1EMHighEt;
  map_BitOfStandardHLTPath["HLT1EMVeryHighEt"] =  HLT1EMVeryHighEt;
  map_BitOfStandardHLTPath["HLT2ElectronZCounter"] =  CandHLT2ElectronZCounter;
  map_BitOfStandardHLTPath["HLT2ElectronExclusive"] =  CandHLT2ElectronExclusive;
  map_BitOfStandardHLTPath["HLT2PhotonExclusive"] =  CandHLT2PhotonExclusive;
  map_BitOfStandardHLTPath["HLT1PhotonL1Isolated"] =  CandHLT1PhotonL1Isolated;
  map_BitOfStandardHLTPath["HLT1MuonIso"] =  HLT1MuonIso;
  map_BitOfStandardHLTPath["HLT1MuonNonIso"] =  HLT1MuonNonIso;
  map_BitOfStandardHLTPath["HLT2MuonIso"] =  CandHLT2MuonIso;
  map_BitOfStandardHLTPath["HLT2MuonNonIso"] =  HLT2MuonNonIso;
  map_BitOfStandardHLTPath["HLT2MuonJPsi"] =  HLT2MuonJPsi;
  map_BitOfStandardHLTPath["HLT2MuonUpsilon"] =  HLT2MuonUpsilon;
  map_BitOfStandardHLTPath["HLT2MuonZ"] =  HLT2MuonZ;
  map_BitOfStandardHLTPath["HLTNMuonNonIso"] =  HLTNMuonNonIso;
  map_BitOfStandardHLTPath["HLT2MuonSameSign"] =  HLT2MuonSameSign;
  map_BitOfStandardHLTPath["HLT1MuonPrescalePt3"] =  CandHLT1MuonPrescalePt3;
  map_BitOfStandardHLTPath["HLT1MuonPrescalePt5"] =  CandHLT1MuonPrescalePt5;
  map_BitOfStandardHLTPath["HLT1MuonPrescalePt7x7"] =  CandHLT1MuonPrescalePt7x7;
  map_BitOfStandardHLTPath["HLT1MuonPrescalePt7x10"] =  CandHLT1MuonPrescalePt7x10;
  map_BitOfStandardHLTPath["HLT1MuonLevel1"] =  CandHLT1MuonLevel1;
  map_BitOfStandardHLTPath["HLTB1Jet"] =  HLTB1Jet;
  map_BitOfStandardHLTPath["HLTB2Jet"] =  HLTB2Jet;
  map_BitOfStandardHLTPath["HLTB3Jet"] =  HLTB3Jet;
  map_BitOfStandardHLTPath["HLTB4Jet"] =  HLTB4Jet;
  map_BitOfStandardHLTPath["HLTBHT"] =  HLTBHT;
  map_BitOfStandardHLTPath["HLTB1JetMu"] =  HLTB1JetMu;
  map_BitOfStandardHLTPath["HLTB2JetMu"] =  HLTB2JetMu;
  map_BitOfStandardHLTPath["HLTB3JetMu"] =  HLTB3JetMu;
  map_BitOfStandardHLTPath["HLTB4JetMu"] =  HLTB4JetMu;
  map_BitOfStandardHLTPath["HLTBHTMu"] =  HLTBHTMu;
  map_BitOfStandardHLTPath["HLTBJPsiMuMu"] =  HLTBJPsiMuMu;
  map_BitOfStandardHLTPath["HLT1Tau"] =  HLT1Tau;
  map_BitOfStandardHLTPath["HLT1Tau1MET"] =  HLT1Tau1MET;
  map_BitOfStandardHLTPath["HLT2TauPixel"] =  HLT2TauPixel;
  map_BitOfStandardHLTPath["HLTXElectronBJet"] =  HLTXElectronBJet;
  map_BitOfStandardHLTPath["HLTXMuonBJet"] =  HLTXMuonBJet;
  map_BitOfStandardHLTPath["HLTXMuonBJetSoftMuon"] =  HLTXMuonBJetSoftMuon;
  map_BitOfStandardHLTPath["HLTXElectron1Jet"] =  HLTXElectron1Jet;
  map_BitOfStandardHLTPath["HLTXElectron2Jet"] =  HLTXElectron2Jet;
  map_BitOfStandardHLTPath["HLTXElectron3Jet"] =  HLTXElectron3Jet;
  map_BitOfStandardHLTPath["HLTXElectron4Jet"] =  HLTXElectron4Jet;
  map_BitOfStandardHLTPath["HLTXMuonJets"] =  HLTXMuonJets;
  map_BitOfStandardHLTPath["HLTXElectronMuon"] =  HLTXElectronMuon;
  map_BitOfStandardHLTPath["HLTXElectronMuonRelaxed"] =  HLTXElectronMuonRelaxed;
  map_BitOfStandardHLTPath["HLTXElectronTau"] =  HLTXElectronTau;
  map_BitOfStandardHLTPath["HLTXMuonTau"] =  HLTXMuonTau;
  map_BitOfStandardHLTPath["HLTHcalIsolatedTrack"] =  CandHLTHcalIsolatedTrack;
  //map_BitOfStandardHLTPath["HLTHcalIsolatedTrackNoEcalIsol"] = CandHLTHcalIsolatedTrackNoEcalIsol;
  map_BitOfStandardHLTPath["HLTMinBiasPixel"] =  HLTMinBiasPixel;
  //map_BitOfStandardHLTPath["HLTMinBiasForAlignment"] = CandHLTMinBiasForAlignment;
  map_BitOfStandardHLTPath["HLTMinBias"] =  HLTMinBias;
  map_BitOfStandardHLTPath["HLTZeroBias"] =  HLTZeroBias;
  map_BitOfStandardHLTPath["HLTMinBiasPixel"] =  HLTMinBiasPixel;
  map_BitOfStandardHLTPath["HLTMinBias"] =  HLTMinBias;
  map_BitOfStandardHLTPath["HLTZeroBias"] =  HLTZeroBias;

}

void OHltTree::SetStandardHLTPath() {

  BitOfStandardHLTPath.push_back(HLT1jet); 
  BitOfStandardHLTPath.push_back(HLT2jet); 
  BitOfStandardHLTPath.push_back(HLT3jet); 
  BitOfStandardHLTPath.push_back(HLT4jet); 
  BitOfStandardHLTPath.push_back(HLT1MET); 
  BitOfStandardHLTPath.push_back(HLT2jetAco); 
  BitOfStandardHLTPath.push_back(HLT1jet1METAco); 
  BitOfStandardHLTPath.push_back(HLT1jet1MET); 
  BitOfStandardHLTPath.push_back(HLT2jet1MET); 
  BitOfStandardHLTPath.push_back(HLT3jet1MET); 
  BitOfStandardHLTPath.push_back(HLT4jet1MET); 
  BitOfStandardHLTPath.push_back(HLT1MET1HT); 
  BitOfStandardHLTPath.push_back(CandHLT1SumET); 
  BitOfStandardHLTPath.push_back(HLT1jetPE1); 
  BitOfStandardHLTPath.push_back(HLT1jetPE3); 
  BitOfStandardHLTPath.push_back(HLT1jetPE5); 
  BitOfStandardHLTPath.push_back(CandHLT1jetPE7); 
  BitOfStandardHLTPath.push_back(CandHLT1METPre1); 
  BitOfStandardHLTPath.push_back(CandHLT1METPre2); 
  BitOfStandardHLTPath.push_back(CandHLT1METPre3); 
  BitOfStandardHLTPath.push_back(CandHLT2jetAve30); 
  BitOfStandardHLTPath.push_back(CandHLT2jetAve60); 
  BitOfStandardHLTPath.push_back(CandHLT2jetAve110); 
  BitOfStandardHLTPath.push_back(CandHLT2jetAve150); 
  BitOfStandardHLTPath.push_back(CandHLT2jetAve200); 
  BitOfStandardHLTPath.push_back(HLT2jetvbfMET); 
  BitOfStandardHLTPath.push_back(HLTS2jet1METNV); 
  BitOfStandardHLTPath.push_back(HLTS2jet1METAco); 
  BitOfStandardHLTPath.push_back(CandHLTSjet1MET1Aco); 
  BitOfStandardHLTPath.push_back(CandHLTSjet2MET1Aco); 
  BitOfStandardHLTPath.push_back(CandHLTS2jetAco); 
  BitOfStandardHLTPath.push_back(CandHLTJetMETRapidityGap); 
  BitOfStandardHLTPath.push_back(HLT1Electron); 
  BitOfStandardHLTPath.push_back(HLT1ElectronRelaxed); 
  BitOfStandardHLTPath.push_back(HLT2Electron); 
  BitOfStandardHLTPath.push_back(HLT2ElectronRelaxed); 
  BitOfStandardHLTPath.push_back(HLT1Photon); 
  BitOfStandardHLTPath.push_back(HLT1PhotonRelaxed); 
  BitOfStandardHLTPath.push_back(HLT2Photon); 
  BitOfStandardHLTPath.push_back(HLT2PhotonRelaxed); 
  BitOfStandardHLTPath.push_back(HLT1EMHighEt); 
  BitOfStandardHLTPath.push_back(HLT1EMVeryHighEt); 
  BitOfStandardHLTPath.push_back(CandHLT2ElectronZCounter); 
  BitOfStandardHLTPath.push_back(CandHLT2ElectronExclusive); 
  BitOfStandardHLTPath.push_back(CandHLT2PhotonExclusive); 
  BitOfStandardHLTPath.push_back(CandHLT1PhotonL1Isolated); 
  BitOfStandardHLTPath.push_back(HLT1MuonIso); 
  BitOfStandardHLTPath.push_back(HLT1MuonNonIso); 
  BitOfStandardHLTPath.push_back(CandHLT2MuonIso); 
  BitOfStandardHLTPath.push_back(HLT2MuonNonIso); 
  BitOfStandardHLTPath.push_back(HLT2MuonJPsi); 
  BitOfStandardHLTPath.push_back(HLT2MuonUpsilon); 
  BitOfStandardHLTPath.push_back(HLT2MuonZ); 
  BitOfStandardHLTPath.push_back(HLTNMuonNonIso); 
  BitOfStandardHLTPath.push_back(HLT2MuonSameSign); 
  BitOfStandardHLTPath.push_back(CandHLT1MuonPrescalePt3); 
  BitOfStandardHLTPath.push_back(CandHLT1MuonPrescalePt5); 
  BitOfStandardHLTPath.push_back(CandHLT1MuonPrescalePt7x7); 
  BitOfStandardHLTPath.push_back(CandHLT1MuonPrescalePt7x10); 
  BitOfStandardHLTPath.push_back(CandHLT1MuonLevel1); 
  BitOfStandardHLTPath.push_back(HLTB1Jet); 
  BitOfStandardHLTPath.push_back(HLTB2Jet); 
  BitOfStandardHLTPath.push_back(HLTB3Jet); 
  BitOfStandardHLTPath.push_back(HLTB4Jet); 
  BitOfStandardHLTPath.push_back(HLTBHT); 
  BitOfStandardHLTPath.push_back(HLTB1JetMu); 
  BitOfStandardHLTPath.push_back(HLTB2JetMu); 
  BitOfStandardHLTPath.push_back(HLTB3JetMu); 
  BitOfStandardHLTPath.push_back(HLTB4JetMu); 
  BitOfStandardHLTPath.push_back(HLTBHTMu); 
  BitOfStandardHLTPath.push_back(HLTBJPsiMuMu); 
  BitOfStandardHLTPath.push_back(HLT1Tau); 
  BitOfStandardHLTPath.push_back(HLT1Tau1MET); 
  BitOfStandardHLTPath.push_back(HLT2TauPixel); 
  BitOfStandardHLTPath.push_back(HLTXElectronBJet); 
  BitOfStandardHLTPath.push_back(HLTXMuonBJet); 
  BitOfStandardHLTPath.push_back(HLTXMuonBJetSoftMuon); 
  BitOfStandardHLTPath.push_back(HLTXElectron1Jet); 
  BitOfStandardHLTPath.push_back(HLTXElectron2Jet); 
  BitOfStandardHLTPath.push_back(HLTXElectron3Jet); 
  BitOfStandardHLTPath.push_back(HLTXElectron4Jet); 
  BitOfStandardHLTPath.push_back(HLTXMuonJets); 
  BitOfStandardHLTPath.push_back(HLTXElectronMuon); 
  BitOfStandardHLTPath.push_back(HLTXElectronMuonRelaxed); 
  BitOfStandardHLTPath.push_back(HLTXElectronTau); 
  BitOfStandardHLTPath.push_back(HLTXMuonTau); 
  BitOfStandardHLTPath.push_back(CandHLTHcalIsolatedTrack); 
  //BitOfStandardHLTPath.push_back(CandHLTHcalIsolatedTrackNoEcalIsol);
  BitOfStandardHLTPath.push_back(HLTMinBiasPixel); 
  //BitOfStandardHLTPath.push_back(CandHLTMinBiasForAlignment);
  BitOfStandardHLTPath.push_back(HLTMinBias); 
  BitOfStandardHLTPath.push_back(HLTZeroBias); 
  BitOfStandardHLTPath.push_back(HLTMinBiasPixel); 
  BitOfStandardHLTPath.push_back(HLTMinBias); 
  BitOfStandardHLTPath.push_back(HLTZeroBias); 

}

Bool_t OHltTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void OHltTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t OHltTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef OHltTree_cxx