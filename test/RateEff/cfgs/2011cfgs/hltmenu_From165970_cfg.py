#------------------------------------------------------
# Configuration file for Rate & Efficiency calculations
#------------------------------------------------------
# This version is compliant with RateEff-02-XX-XX
# using logical parser for L1 seeds
#

##########################################
# General Menu & Run conditions
##########################################
run:{
    nEntries = -1;
    nPrintStatusEvery = 10000; # print out status every n events processed
    menuTag  = "HLT_Menu1E33v1May19Train";
    alcaCondition = "startup";
    versionTag  = "20110429_DS_Global"; # using column 0 of the menu 
    isRealData = true;
    doPrintAll = true;
    doDeterministicPrescale =true;
    dsList = "Datasets.list";
    readRefPrescalesFromNtuple = true;

};

########################################## 
# Run information for real data 
########################################## 
data:{ 
 # Enter the length of 1 lumi section and prescale factor of the dataset
 lumiSectionLength = 23.3;
 lumiScaleFactor = 1.35; # from 1.04(run 165970) to 1.4(May19th)

 prescaleNormalization = 1; #

 runLumiblockList = ( 
    (165970, 38, 329 ) # (runnr, minLumiBlock, maxLumiBlock)
 );

};

##########################################
# Beam conditions
##########################################
beam:{
 bunchCrossingTime = 25.0E-09; # Design: 25 ns Startup: 75 ns
 iLumi = 1E33;
 maxFilledBunches = 3564;
 nFilledBunches = 800;
 cmsEnergy = 7.; # Collision energy in TeV
};

##########################################
# Samples & Processes
##########################################
process:{
 isPhysicsSample = [0]; #Must be an int type
 names = ["minbias"];
 fnames = ["openhlt_*.root"];

#fnal
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__Commissioning_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__Cosmics_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__DoubleElectron_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__DoubleMu_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__ElectronHad_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__HT_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__Photon_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__PhotonHad_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__SingleElectron_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__Jet_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__MultiJet_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__MET_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__BTag_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__MinimumBias_Run2011A-v1__20110601_0529/"];
paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__SingleMu_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__MuOnia_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__MuEG_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__MuHad_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__Tau_Run2011A-v1__20110601_0529/"];
## paths = ["dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpctrig/Commish2011/r165970__TauPlusX_Run2011A-v1__20110601_0529/"];

#castor, cern
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__Commissioning_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__Cosmics_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__DoubleElectron_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__DoubleMu_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__ElectronHad_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__HT_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__Photon_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__PhotonHad_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__SingleElectron_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__Jet_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__MultiJet_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__MET_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__BTag_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__MinimumBias_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__SingleMu_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__MuOnia_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__MuEG_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__MuHad_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__Tau_Run2011A-v1__20110601_0529/"];
## paths = ["rfio:/castor/cern.ch/user/l/lucieg/OpenHLT/Commish2011/r165970__TauPlusX_Run2011A-v1__20110601_0529/"];



 doMuonCuts = [false];
 doElecCuts = [false];
 sigmas = [7.13E10]; # xsecs * filter efficiencies for minbias
};


##########################################
# Menu
##########################################
menu:{
 isL1Menu = false; # Default is false: is HLT Menu
 doL1preloop = true; 

  ## preFilterByBits = "";



  # (TriggerName, Prescale, EventSize)
 triggers = (
#
############# dataset Commissioning ###############
##    ("HLT_Activity_Ecal_SC7_v5", "L1_ZeroBias_Ext", 1, 0.15),
##    ("HLT_L1SingleJet16_v2", "L1_SingleJet16", 1, 0.15),
##    ("HLT_L1SingleJet36_v2", "L1_SingleJet36", 1, 0.15),
##    ("HLT_L1SingleMuOpen_v2", "L1_SingleMuOpen", 1, 0.15),
##    ("HLT_L1SingleMuOpen_DT_v2", "L1_SingleMuOpen", 1, 0.15),
##    ("HLT_Mu5_TkMu0_OST_Jpsi_Tight_B5Q7_v4", "L1_SingleMu5_Eta1p5_Q80", 1, 0.15),
##    ("HLT_L1SingleEG5_v2", "L1_SingleEG5", 1, 0.15),
##    ("HLT_L1SingleEG12_v2", "L1_SingleEG12", 1, 0.15),
##    ("HLT_BeamGas_HF_v5", "L1_BeamGas_Hf", 1, 0.15),
##    ("HLT_BeamGas_BSC_v3", "L1_BeamGas_Bsc", 1, 0.15),
##    ("HLT_L1_PreCollisions_v2", "L1_PreCollisions", 1, 0.15),
##    ("HLT_L1_Interbunch_BSC_v2", "L1_InterBunch_Bsc", 1, 0.15),
##    ("HLT_IsoTrackHE_v5", "L1_SingleJet68", 1, 0.15),
##    ("HLT_IsoTrackHB_v4", "L1_SingleJet68", 1, 0.15)##,
## ############# dataset Cosmics ###############
##    ("HLT_BeamHalo_v3", "L1_BeamHalo", 1, 0.15),
##    ("HLT_L1SingleMuOpen_AntiBPTX_v2", "L1_SingleMuOpen", 1, 0.15),
##    ("HLT_L1TrackerCosmics_v3", "L1Tech_RPC_TTU_pointing_Cosmics.v0", 1, 0.15),
##    ("HLT_RegionalCosmicTracking_v4", "L1Tech_RPC_TTU_pointing_Cosmics.v0 AND L1_SingleMuOpen", 1, 0.15)##,
## ############# dataset DoubleElectron ###############
##    ("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v5", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele8_v5", "L1_SingleEG5", 1, 0.15),
##    ("HLT_Ele8_CaloIdL_CaloIsoVL_v5", "L1_SingleEG5", 1, 0.15),
##    ("HLT_Ele8_CaloIdL_TrkIdVL_v5", "L1_SingleEG5", 1, 0.15),
##    ("HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4", "L1_SingleEG5", 1, 0.15),
##    ("HLT_Ele17_CaloIdL_CaloIsoVL_v5", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v5", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v3", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v6", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFT_v1", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_DoubleEle45_CaloIdL_v1", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v5", "L1_SingleEG5", 1, 0.15),
##    ("HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v6", "L1_TripleEG5", 1, 0.15),
##    ("HLT_TripleEle10_CaloIdL_TrkIdVL_v6", "L1_TripleEG5", 1, 0.15)##,
## ############# dataset DoubleMu ###############
##    ("HLT_L1DoubleMu0_v2", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_L2DoubleMu0_v4", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_L2DoubleMu23_NoVertex_v3", "L1_DoubleMu3", 1, 0.15),
##    ("HLT_DoubleMu3_v5", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_DoubleMu6_v3", "L1_DoubleMu3", 1, 0.15),
##    ("HLT_DoubleMu7_v3", "L1_DoubleMu3", 1, 0.15),
##    ("HLT_DoubleMu45_v1", "L1_DoubleMu3", 1, 0.15),
##    ("HLT_DoubleMu4_Acoplanarity03_v4", "L1_DoubleMu3", 1, 0.15),
##    ("HLT_DoubleMu5_Acoplanarity03_v1", "L1_DoubleMu3", 1, 0.15),
##    ("HLT_Mu13_Mu8_v2", "L1_DoubleMu3", 1, 0.15),
##    ("HLT_Mu17_Mu8_v2", "L1_DoubleMu3", 1, 0.15),
##    ("HLT_TripleMu5_v4", "L1_DoubleMu3", 1, 0.15),
##    ("HLT_Mu8_Jet40_v6", "L1_Mu3_Jet20_Central", 1, 0.15)##,
## ############# dataset ElectronHad ###############
##    ("HLT_DoubleEle8_CaloIdL_TrkIdVL_v2", "L1_DoubleEG5", 1, 0.15),
##    ("HLT_HT250_Ele5_CaloIdVL_TrkIdVL_CaloIsoVL_TrkIsoVL_PFMHT35_v5", "L1_HTT100", 1, 0.15),
##    ("HLT_HT300_Ele5_CaloIdVL_TrkIdVL_CaloIsoVL_TrkIsoVL_PFMHT40_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT350_Ele5_CaloIdVL_TrkIdVL_CaloIsoVL_TrkIsoVL_PFMHT45_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_Ele8_CaloIdT_TrkIdT_DiJet30_v2", "L1_DoubleEG5_HTT50", 1, 0.15),
##    ("HLT_Ele8_CaloIdT_TrkIdT_TriJet30_v2", "L1_DoubleEG5_HTT50", 1, 0.15),
##    ("HLT_Ele8_CaloIdT_TrkIdT_QuadJet30_v2", "L1_DoubleEG5_HTT50", 1, 0.15),
##    ("HLT_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v5", "L1_EG5_HTT75", 1, 0.15),
##    ("HLT_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT250_v5", "L1_EG5_HTT75", 1, 0.15),
##    ("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R005_MR200_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_Ele10_CaloIdL_TrkIdVL_CaloIsoVL_TrkIsoVL_R025_MR200_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_Ele10_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_R020_MR200_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_Jet35_Jet25_Deta3_v4", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_Jet35_Jet25_Deta2_v4", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele15_CaloIdVT_TrkIdT_Jet35_Jet25_Deta2_v4", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v5", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele25_CaloIdVT_TrkIdT_DiCentralJet30_v4", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele25_CaloIdVT_TrkIdT_TriCentralJet30_v4", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele25_CaloIdVT_TrkIdT_QuadCentralJet30_v1", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_v1", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_DiCentralJet30_v1", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30_v1", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_QuadCentralJet30_v1", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_BTagIP_v1", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_BTagIP_v5", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele17_CaloIdVT_TrkIdT_CentralJet30_CentralJet25_v4", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele17_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT15_v4", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT20_v4", "L1_SingleEG12", 1, 0.15),
##    ("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT150_v3", "L1_DoubleEG5_HTT50", 1, 0.15),
##    ("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT150_v3", "L1_DoubleEG5_HTT50", 1, 0.15)##,
## ############# dataset HT ###############
##    ("HLT_DiJet130_PT130_v3", "L1_SingleJet68", 1, 0.15),
##    ("HLT_DiJet160_PT160_v3", "L1_SingleJet92", 1, 0.15),
##    ("HLT_HT150_v5", "L1_HTT50", 1, 0.15),
##    ("HLT_HT150_AlphaT0p60_v4", "L1_HTT75", 1, 0.15),
##    ("HLT_HT200_v5", "L1_HTT75", 1, 0.15),
##    ("HLT_HT200_AlphaT0p53_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT200_AlphaT0p60_v4", "L1_HTT75", 1, 0.15),
##    ("HLT_HT250_v5", "L1_HTT100", 1, 0.15),
##    ("HLT_HT250_AlphaT0p53_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT250_AlphaT0p54_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT250_DoubleDisplacedJet60_v4", "L1_HTT100", 1, 0.15),
##    ("HLT_HT250_MHT60_v6", "L1_HTT100", 1, 0.15),
##    ("HLT_HT250_MHT70_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT250_MHT80_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT300_v6", "L1_HTT100", 1, 0.15),
##    ("HLT_HT300_MHT75_v7", "L1_HTT100", 1, 0.15),
##    ("HLT_HT300_PFMHT55_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT300_CentralJet30_BTagIP_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT300_CentralJet30_BTagIP_PFMHT55_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT300_CentralJet30_BTagIP_PFMHT75_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT300_AlphaT0p52_v4", "L1_HTT100", 1, 0.15),
##    ("HLT_HT300_AlphaT0p53_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT350_v5", "L1_HTT100", 1, 0.15),
##    ("HLT_HT350_AlphaT0p51_v4", "L1_HTT100", 1, 0.15),
##    ("HLT_HT350_AlphaT0p53_v4", "L1_HTT100", 1, 0.15),
##    ("HLT_HT400_v5", "L1_HTT100", 1, 0.15),
##    ("HLT_HT400_AlphaT0p51_v4", "L1_HTT100", 1, 0.15),
##    ("HLT_HT450_v5", "L1_HTT100", 1, 0.15),
##    ("HLT_HT500_v5", "L1_HTT100", 1, 0.15),
##    ("HLT_HT550_v5", "L1_HTT100", 1, 0.15),
##    ("HLT_R014_MR150_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R014_MR150_CentralJet40_BTagIP_v4", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R014_MR450_CentralJet40_BTagIP_v4", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R020_MR150_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R020_MR350_CentralJet40_BTagIP_v4", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R020_MR500_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R020_MR550_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R025_MR150_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R025_MR250_CentralJet40_BTagIP_v4", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R025_MR400_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R025_MR450_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R033_MR300_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R033_MR350_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R038_MR200_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_R038_MR250_v3", "L1_DoubleJet36_Central", 1, 0.15)##,
## ############# dataset Photon ###############
##    ("HLT_Photon20_CaloIdVL_IsoL_v4", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Photon20_R9Id_Photon18_R9Id_v5", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Photon26_Photon18_v5", "L1_SingleEG15", 1, 0.15),
##    ("HLT_Photon26_IsoVL_Photon18_v5", "L1_SingleEG15", 1, 0.15),
##    ("HLT_Photon26_IsoVL_Photon18_IsoVL_v5", "L1_SingleEG15", 1, 0.15),
##    ("HLT_Photon26_CaloIdL_IsoVL_Photon18_v5", "L1_SingleEG15", 1, 0.15),
##    ("HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v4", "L1_SingleEG15", 1, 0.15),
##    ("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v5", "L1_SingleEG15", 1, 0.15),
##    ("HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v4", "L1_SingleEG15", 1, 0.15),
##    ("HLT_Photon26_R9Id_Photon18_R9Id_v2", "L1_SingleEG15", 1, 0.15),
##    ("HLT_Photon30_CaloIdVL_v5", "L1_SingleEG15", 1, 0.15),
##    ("HLT_Photon30_CaloIdVL_IsoL_v5", "L1_SingleEG15", 1, 0.15),
##    ("HLT_Photon36_IsoVL_Photon22_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon36_CaloIdL_Photon22_CaloIdL_v4", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon36_CaloIdL_IsoVL_Photon22_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v1", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v1", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon36_CaloId_IsoVL_Photon22_R9Id_v1", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL_v1", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon36_R9Id_Photon22_R9Id_v1", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon40_CaloIdL_Photon28_CaloIdL_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon50_CaloIdVL_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon50_CaloIdVL_IsoL_v4", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon75_CaloIdVL_v5", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon75_CaloIdVL_IsoL_v5", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon90_CaloIdVL_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon90_CaloIdVL_IsoL_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon125_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon200_NoHE_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_DoublePhoton33_v5", "L1_SingleEG20", 1, 0.15),
##    ("HLT_DoublePhoton33_HEVT_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_DoublePhoton50_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_DoublePhoton60_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_DoublePhoton5_IsoVL_CEP_v4", "L1_DoubleEG2_FwdVeto", 1, 0.15),
##    ("HLT_DoubleEle33_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_DoubleEle33_CaloIdL_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_DoublePhoton40_MR150_v3", "L1_SingleEG20", 1, 0.15),
##    ("HLT_DoublePhoton40_R014_MR150_v3", "L1_SingleEG20", 1, 0.15)##,
## ############# dataset PhotonHad ###############
##    ("HLT_Photon70_CaloIdL_HT300_v6", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon70_CaloIdL_HT350_v5", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon70_CaloIdL_MHT50_v6", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon70_CaloIdL_MHT70_v5", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon40_R005_MR150_v3", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon40_R014_MR450_v3", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon40_R020_MR300_v3", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon40_R025_MR200_v3", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Photon40_R038_MR150_v3", "L1_SingleEG20", 1, 0.15)##,
## ############# dataset SingleElectron ###############
##    ("HLT_Ele25_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v1", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele25_WP80_PFMT40_v1", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele27_WP70_PFMT40_PFMHT20_v1", "L1_SingleEG15", 1, 0.15),
##    ("HLT_Ele32_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Ele42_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL_v1", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Ele52_CaloIdVT_TrkIdT_v2", "L1_SingleEG20", 1, 0.15),
##    ("HLT_Ele65_CaloIdVT_TrkIdT_v1", "L1_SingleEG20", 1, 0.15)##,
## ############# dataset Jet ###############
##    ("HLT_Jet30_v4", "L1_SingleJet16", 1, 0.15),
##    ("HLT_Jet60_v4", "L1_SingleJet36", 1, 0.15),
##    ("HLT_Jet80_v4", "L1_SingleJet52", 1, 0.15),
##    ("HLT_Jet110_v4", "L1_SingleJet68", 1, 0.15),
##    ("HLT_Jet150_v4", "L1_SingleJet92", 1, 0.15),
##    ("HLT_Jet190_v4", "L1_SingleJet92", 1, 0.15),
##    ("HLT_Jet240_v4", "L1_SingleJet92", 1, 0.15),
##    ("HLT_Jet300_v3", "L1_SingleJet92", 1, 0.15),
##    ("HLT_Jet370_v4", "L1_SingleJet92", 1, 0.15),
##    ("HLT_Jet370_NoJetID_v4", "L1_SingleJet92", 1, 0.15),
##    ("HLT_DiJetAve30_v4", "L1_SingleJet16", 1, 0.15),
##    ("HLT_DiJetAve60_v4", "L1_SingleJet36", 1, 0.15),
##    ("HLT_DiJetAve80_v4", "L1_SingleJet52", 1, 0.15),
##    ("HLT_DiJetAve110_v4", "L1_SingleJet68", 1, 0.15),
##    ("HLT_DiJetAve150_v4", "L1_SingleJet92", 1, 0.15),
##    ("HLT_DiJetAve190_v4", "L1_SingleJet92", 1, 0.15),
##    ("HLT_DiJetAve240_v4", "L1_SingleJet92", 1, 0.15),
##    ("HLT_DiJetAve300_v4", "L1_SingleJet92", 1, 0.15),
##    ("HLT_DiJetAve370_v4", "L1_SingleJet92", 1, 0.15)##,
## ############# dataset MultiJet ###############
##    ("HLT_DoubleJet30_ForwardBackward_v5", "L1_DoubleForJet32_EtaOpp", 1, 0.15),
##    ("HLT_DoubleJet60_ForwardBackward_v5", "L1_DoubleForJet32_EtaOpp", 1, 0.15),
##    ("HLT_DoubleJet70_ForwardBackward_v5", "L1_DoubleForJet32_EtaOpp", 1, 0.15),
##    ("HLT_DoubleJet80_ForwardBackward_v5", "L1_DoubleForJet44_EtaOpp", 1, 0.15),
##    ("HLT_CentralJet46_BTagIP3D_CentralJet38_BTagIP3D_v1", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_QuadJet40_v5", "L1_QuadJet20_Central", 1, 0.15),
##    ("HLT_QuadJet40_IsoPFTau40_v7", "L1_QuadJet20_Central", 1, 0.15),
##    ("HLT_QuadJet45_IsoPFTau45_v2", "L1_QuadJet20_Central", 1, 0.15),
##    ("HLT_QuadJet50_Jet40_Jet30_v1", "L1_QuadJet20_Central", 1, 0.15),
##    ("HLT_QuadJet60_v4", "L1_QuadJet20_Central", 1, 0.15),
##    ("HLT_QuadJet70_v4", "L1_QuadJet20_Central", 1, 0.15),
##    ("HLT_ExclDiJet60_HFOR_v4", "L1_SingleJet36", 1, 0.15),
##    ("HLT_ExclDiJet60_HFAND_v4", "L1_SingleJet36_FwdVeto", 1, 0.15),
##    ("HLT_L1ETM30_v2", "L1_ETM30", 1, 0.15),
##    ("HLT_L1DoubleJet36Central_v2", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_L1MultiJet_v2", "L1_HTT50 OR L1_TripleJet28_Central OR L1_QuadJet20_Central", 1, 0.15)##,
## ############# dataset MET ###############
##    ("HLT_CentralJet80_MET65_v4", "L1_ETM30", 1, 0.15),
##    ("HLT_CentralJet80_MET80HF_v3", "L1_ETM30", 1, 0.15),
##    ("HLT_CentralJet80_MET100_v4", "L1_ETM30", 1, 0.15),
##    ("HLT_CentralJet80_MET160_v4", "L1_ETM30", 1, 0.15),
##    ("HLT_DiJet60_MET45_v4", "L1_ETM20", 1, 0.15),
##    ("HLT_DiCentralJet20_MET80_v2", "L1_ETM30", 1, 0.15),
##    ("HLT_DiCentralJet20_BTagIP_MET65_v3", "L1_ETM30", 1, 0.15),
##    ("HLT_PFMHT150_v7", "L1_ETM30", 1, 0.15),
##    ("HLT_MET65_v1", "L1_ETM30", 1, 0.15),
##    ("HLT_MET65_HBHENoiseFiltered_v1", "L1_ETM30", 1, 0.15),
##    ("HLT_MET100_v4", "L1_ETM30", 1, 0.15),
##    ("HLT_MET100_HBHENoiseFiltered_v2", "L1_ETM30", 1, 0.15),
##    ("HLT_MET120_v4", "L1_ETM30", 1, 0.15),
##    ("HLT_MET120_HBHENoiseFiltered_v2", "L1_ETM30", 1, 0.15),
##    ("HLT_MET200_v4", "L1_ETM30", 1, 0.15),
##    ("HLT_MET200_HBHENoiseFiltered_v2", "L1_ETM30", 1, 0.15),
##    ("HLT_L2Mu60_1Hit_MET40_v1", "L1_SingleMu20", 1, 0.15),
##    ("HLT_L2Mu60_1Hit_MET60_v1", "L1_SingleMu20", 1, 0.15)##,
## ############# dataset BTag ###############
##    ("HLT_BTagMu_DiJet20_Mu5_v5", "L1_Mu3_Jet16_Central", 1, 0.15),
##    ("HLT_BTagMu_DiJet40_Mu5_v5", "L1_Mu3_Jet20_Central", 1, 0.15),
##    ("HLT_BTagMu_DiJet70_Mu5_v5", "L1_Mu3_Jet28_Central", 1, 0.15),
##    ("HLT_BTagMu_DiJet110_Mu5_v5", "L1_Mu3_Jet28_Central", 1, 0.15)##,
## ############# dataset MinimumBias ###############
##    ("HLT_JetE30_NoBPTX_v4", "L1_SingleJet20_NotBptxOR", 1, 0.15),
##    ("HLT_JetE30_NoBPTX_NoHalo_v6", "L1_SingleJet20_NotBptxOR_NotMuBeamHalo", 1, 0.15),
##    ("HLT_JetE30_NoBPTX3BX_NoHalo_v6", "L1_SingleJet20_NotBptxOR_NotMuBeamHalo", 1, 0.15),
##    ("HLT_JetE50_NoBPTX3BX_NoHalo_v2", "L1_SingleJet32_NotBptxOR_NotMuBeamHalo", 1, 0.15),
##    ("HLT_PixelTracks_Multiplicity80_v3", "L1_ETT220", 1, 0.15),
##    ("HLT_PixelTracks_Multiplicity100_v3", "L1_ETT220", 1, 0.15),
##    ("HLT_ZeroBias_v3", "L1_ZeroBias_Ext", 1, 0.15),
##    ("HLT_Physics_v1", "", 1, 0.15),
##    ("HLT_Random_v1", "", 1, 0.15)##,
############# dataset SingleMu ###############
   ("HLT_L1SingleMu10_v2", "L1_SingleMu10", 1, 0.15),
   ("HLT_L1SingleMu20_v2", "L1_SingleMu20", 1, 0.15),
   ("HLT_L2Mu10_v3", "L1_SingleMu10", 1, 0.15),
   ("HLT_L2Mu20_v3", "L1_SingleMu12", 1, 0.15),
   ("HLT_Mu3_v5", "L1_SingleMuOpen", 1, 0.15),
   ("HLT_Mu5_v5", "L1_SingleMu3", 1, 0.15),
   ("HLT_Mu8_v3", "L1_SingleMu3", 1, 0.15),
   ("HLT_Mu12_v3", "L1_SingleMu7", 1, 0.15),
   ("HLT_Mu15_v4", "L1_SingleMu10", 1, 0.15),
   ("HLT_Mu20_v3", "L1_SingleMu12", 1, 0.15),
   ("HLT_Mu24_v3", "L1_SingleMu12", 1, 0.15),
   ("HLT_Mu30_v3", "L1_SingleMu12", 1, 0.15),
   ("HLT_Mu40_v1", "L1_SingleMu16", 1, 0.15),
   ("HLT_Mu100_v1", "L1_SingleMu16", 1, 0.15),
   ("HLT_IsoMu12_v5", "L1_SingleMu7", 1, 0.15),
   ("HLT_IsoMu15_v9", "L1_SingleMu10", 1, 0.15),
   ("HLT_IsoMu17_v9", "L1_SingleMu10", 1, 0.15),
   ("HLT_IsoMu24_v5", "L1_SingleMu12", 1, 0.15),
   ("HLT_IsoMu30_v5", "L1_SingleMu12", 1, 0.15)##,
## ############# dataset MuOnia ###############
##    ("HLT_DoubleMu2_Bs_v3", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Dimuon0_Jpsi_v1", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Dimuon0_Upsilon_v1", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Dimuon4_Bs_Barrel_v3", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Dimuon5_Upsilon_Barrel_v1", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Dimuon6_Bs_v2", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Dimuon7_LowMass_Displaced_v2", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Dimuon7_Jpsi_Displaced_v1", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Dimuon7_Jpsi_X_Barrel_v1", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Dimuon7_PsiPrime_v1", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Dimuon10_Jpsi_Barrel_v1", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Dimuon0_Jpsi_Muon_v2", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Dimuon0_Upsilon_Muon_v2", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Mu5_L2Mu2_Jpsi_v4", "L1_DoubleMu0", 1, 0.15),
##    ("HLT_Mu5_Track2_Jpsi_v4", "L1_SingleMu3", 1, 0.15),
##    ("HLT_Mu7_Track7_Jpsi_v5", "L1_SingleMu7", 1, 0.15)##,
############# dataset MuEG ###############
##    ("HLT_Mu5_DoubleEle8_CaloIdL_TrkIdVL_v2", "L1_MuOpen_EG5", 1, 0.15),
##    ("HLT_Mu8_Ele17_CaloIdL_v5", "L1_MuOpen_EG5", 1, 0.15),
##    ("HLT_Mu8_Photon20_CaloIdVT_IsoT_v5", "L1_MuOpen_EG5", 1, 0.15),
##    ("HLT_Mu15_Photon20_CaloIdL_v6", "L1_MuOpen_EG5", 1, 0.15),
##    ("HLT_Mu15_DoublePhoton15_CaloIdL_v6", "L1_MuOpen_EG5", 1, 0.15),
##    ("HLT_Mu17_Ele8_CaloIdL_v5", "L1_MuOpen_EG5", 1, 0.15),
##    ("HLT_DoubleMu5_Ele8_v6", "L1_MuOpen_EG5", 1, 0.15),
##    ("HLT_DoubleMu5_Ele8_CaloIdL_TrkIdVL_v6", "L1_MuOpen_EG5", 1, 0.15)##,
## ############# dataset MuHad ###############
##    ("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT150_v3", "L1_Mu0_HTT50", 1, 0.15),
##    ("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT150_v3", "L1_Mu0_HTT50", 1, 0.15),
##    ("HLT_Mu8_R005_MR200_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_Mu8_R020_MR200_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_Mu8_R025_MR200_v3", "L1_DoubleJet36_Central", 1, 0.15),
##    ("HLT_HT250_Mu5_PFMHT35_v5", "L1_HTT100", 1, 0.15),
##    ("HLT_HT250_Mu15_PFMHT20_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT300_Mu5_PFMHT40_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT350_Mu5_PFMHT45_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_Mu3_DiJet30_v2", "L1_Mu3_Jet20_Central", 1, 0.15),
##    ("HLT_Mu3_TriJet30_v2", "L1_Mu3_Jet20_Central", 1, 0.15),
##    ("HLT_Mu3_QuadJet30_v2", "L1_Mu3_Jet20_Central", 1, 0.15),
##    ("HLT_Mu17_CentralJet30_v6", "L1_SingleMu10", 1, 0.15),
##    ("HLT_Mu17_DiCentralJet30_v6", "L1_SingleMu10", 1, 0.15),
##    ("HLT_Mu17_TriCentralJet30_v6", "L1_SingleMu10", 1, 0.15),
##    ("HLT_Mu17_QuadCentralJet30_v1", "L1_SingleMu10", 1, 0.15),
##    ("HLT_Mu12_DiCentralJet30_BTagIP3D_v1", "L1_SingleMu10", 1, 0.15),
##    ("HLT_Mu17_CentralJet30_BTagIP_v5", "L1_SingleMu10", 1, 0.15),
##    ("HLT_Mu15_HT200_v3", "L1_Mu0_HTT50", 1, 0.15),
##    ("HLT_Mu20_HT200_v3", "L1_Mu0_HTT50", 1, 0.15),
##    ("HLT_IsoMu17_CentralJet30_v1", "L1_SingleMu10", 1, 0.15),
##    ("HLT_IsoMu17_DiCentralJet30_v1", "L1_SingleMu10", 1, 0.15),
##    ("HLT_IsoMu17_TriCentralJet30_v1", "L1_SingleMu10", 1, 0.15),
##    ("HLT_IsoMu17_QuadCentralJet30_v1", "L1_SingleMu10", 1, 0.15),
##    ("HLT_IsoMu17_CentralJet30_BTagIP_v5", "L1_SingleMu10", 1, 0.15),
##    ("HLT_DoubleMu3_HT150_v3", "L1_Mu0_HTT50", 1, 0.15),
##    ("HLT_DoubleMu3_HT200_v6", "L1_Mu0_HTT50", 1, 0.15)##,
## ############# dataset Tau ###############
##    ("HLT_IsoPFTau35_Trk20_v2", "L1_SingleTauJet52 OR L1_SingleJet68", 1, 0.15),
##    ("HLT_IsoPFTau35_Trk20_MET60_v2", "L1_SingleTauJet52 OR L1_SingleJet68", 1, 0.15),
##    ("HLT_IsoPFTau45_Trk20_MET60_v2", "L1_SingleTauJet68 OR L1_SingleJet92", 1, 0.15),
##    ("HLT_DoubleIsoPFTau35_Trk5_eta2p1_v2", "L1_DoubleTauJet28 OR L1_DoubleJet52", 1, 0.15),
##    ("HLT_DoubleIsoPFTau40_Trk5_eta2p1_v2", "L1_DoubleTauJet36 OR L1_DoubleJet52", 1, 0.15)##,
## ############# dataset TauPlusX ###############
##    ("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5", "L1_SingleEG12", 1, 0.15),
##    ("HLT_HT250_DoubleIsoPFTau10_Trk3_PFMHT35_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT300_DoubleIsoPFTau10_Trk3_PFMHT40_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_HT350_DoubleIsoPFTau10_Trk3_PFMHT45_v3", "L1_HTT100", 1, 0.15),
##    ("HLT_Mu15_LooseIsoPFTau15_v4", "L1_SingleMu10", 1, 0.15),
##    ("HLT_IsoMu15_LooseIsoPFTau15_v4", "L1_SingleMu10", 1, 0.15),
##    ("HLT_IsoMu15_LooseIsoPFTau20_v2", "L1_SingleMu10", 1, 0.15),
##    ("HLT_IsoMu15_TightIsoPFTau20_v2", "L1_SingleMu10", 1, 0.15),
##    ("HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20_v2", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v8", "L1_SingleEG12", 1, 0.15),
##    ("HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2", "L1_SingleEG15", 1, 0.15)##,
## # 
 );

 # For L1 prescale preloop to be used in HLT mode only
 L1triggers = ( 
#
  ("L1_SingleMu7", 1),
  ("L1_SingleMu10", 1),
  ("L1_SingleMu12", 1),
  ("L1_SingleMu16", 1),
  ("L1_SingleMu20", 1),
  ("L1_SingleJet36_FwdVeto", 1),
  ("L1_BeamGas_Bsc", 1),
  ("L1_ETT220", 1),
  ("L1_SingleEG5", 1),
  ("L1_SingleEG20", 1),
  ("L1_Mu3_Jet20_Central", 1),
  ("L1_HTT50", 1),
  ("L1_TripleJet28_Central", 1),
  ("L1_QuadJet20_Central", 1),
  ("L1_DoubleForJet32_EtaOpp", 1),
  ("L1_SingleMuOpen", 1),
  ("L1_TripleEG5", 1),
  ("L1Tech_RPC_TTU_pointing_Cosmics.v0", 1),
  ("L1_HTT75", 1),
  ("L1_SingleTauJet68", 1),
  ("L1_SingleJet92", 1),
  ("L1_SingleJet52", 1),
  ("L1_ETM30", 1),
  ("L1_SingleJet20_NotBptxOR_NotMuBeamHalo", 1),
  ("L1_BeamHalo", 1),
  ("L1_DoubleJet36_Central", 1),
  ("L1_EG5_HTT75", 1),
  ("L1_ETM20", 1),
  ("L1_SingleJet36", 1),
  ("L1_SingleJet16", 1),
  ("L1_SingleJet68", 1),
  ("L1_SingleJet128", 1),
  ("L1_SingleTauJet52", 1),
  ("L1_SingleMu3", 1),
  ("L1_SingleIsoEG12", 1),
  ("L1_SingleEG12", 1),
  ("L1_SingleEG15", 1),
  ("L1_SingleEG30", 1),
  ("L1_ZeroBias_Ext", 1),
  ("L1Tech_HCAL_HO_totalOR.v0", 1),
  ("L1Tech_HCAL_HBHE_totalOR.v0", 1),
  ("L1_InterBunch_Bsc", 1),
  ("L1_Mu3_Jet16_Central", 1),
  ("L1_Mu3_Jet28_Central", 1),
  ("L1_SingleJet20_NotBptxOR", 1),
  ("L1_PreCollisions", 1),
  ("L1_MuOpen_EG5", 1),
  ("L1_SingleJet32_NotBptxOR_NotMuBeamHalo", 1),
  ("L1_DoubleMu0", 1),
  ("L1_DoubleMu3", 1),
  ("L1_SingleMu5_Eta1p5_Q80", 1),
  ("L1_HTT100", 1),
  ("L1_DoubleTauJet36", 1),
  ("L1_DoubleJet52", 1),
  ("L1_DoubleEG5_HTT50", 1),
  ("L1_DoubleEG10", 1),
  ("L1_DoubleEG2_FwdVeto", 1),
  ("L1_DoubleEG3", 1),
  ("L1_DoubleEG5", 1),
  ("L1_DoubleEG8", 1),
  ("L1_DoubleEG_12_5", 1),
  ("L1_DoubleIsoEG10", 1),
  ("L1_SingleEG12_Eta2p17", 1),
  ("L1_SingleIsoEG12_Eta2p17", 1),
  ("L1_SingleMu25", 1),
  ("L1_DoubleMu5", 1),
  ("L1_DoubleForJet44_EtaOpp", 1),
  ("L1_BeamGas_Hf", 1),
  ("L1_DoubleTauJet28", 1),
  ("L1_Mu0_HTT50", 1),
  ("L1_DoubleEG5_HTT75", 1),
  ("L1_EG5_HTT100", 1),
  ("L1_EG5_HTT125", 1),
  ("L1_SingleJet80_Central", 1),
  ("L1_TripleEG7", 1)
# 
 );

};

##########################################
#
# Only for experts:
# Select certain branches to speed up code.
# Modify only if you know what you do!
#
##########################################
branch:{
  doSelectBranches = true; #only set to true if you really know what you do!
  selectBranchL1 = true; 
  selectBranchHLT = true;
  selectBranchOpenHLT = true; 
  selectBranchReco = true;
  selectBranchL1extra = true; 
  selectBranchMC = false; 
};
### eof
