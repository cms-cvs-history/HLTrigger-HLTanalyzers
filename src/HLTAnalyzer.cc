// File: HLTAnalyzer.cc
// Description:  Example of Analysis driver originally from Jeremy Mans, 
// Date:  13-October-2006

#include <boost/foreach.hpp>

#include "HLTrigger/HLTanalyzers/interface/HLTAnalyzer.h"
#include "HLTMessages.h"

typedef std::pair<const char *, const edm::InputTag *> MissingCollectionInfo;
  
template <class T>
static inline
bool getCollection(const edm::Event & event, std::vector<MissingCollectionInfo> & missing, edm::Handle<T> & handle, const edm::InputTag & name, const char * description) 
{
  event.getByLabel(name, handle);
  bool valid = handle.isValid();
  if (not valid) {
    missing.push_back( std::make_pair(description, & name) );
    handle.clear();
  }
  return valid;
}

// Boiler-plate constructor definition of an analyzer module:
HLTAnalyzer::HLTAnalyzer(edm::ParameterSet const& conf) {

  // If your module takes parameters, here is where you would define
  // their names and types, and access them to initialize internal
  // variables. Example as follows:
  std::cout << " Beginning HLTAnalyzer Analysis " << std::endl;

  recjets_          = conf.getParameter<edm::InputTag> ("recjets");
  genjets_          = conf.getParameter<edm::InputTag> ("genjets");
  recmet_           = conf.getParameter<edm::InputTag> ("recmet");
  genmet_           = conf.getParameter<edm::InputTag> ("genmet");
  ht_               = conf.getParameter<edm::InputTag> ("ht");
  calotowers_       = conf.getParameter<edm::InputTag> ("calotowers");
  Electron_         = conf.getParameter<edm::InputTag> ("Electron");
  Photon_           = conf.getParameter<edm::InputTag> ("Photon");
  muon_             = conf.getParameter<edm::InputTag> ("muon");
  mctruth_          = conf.getParameter<edm::InputTag> ("mctruth");
  genEventScale_    = conf.getParameter<edm::InputTag> ("genEventScale");

  // keep this separate from l1extramc_ as needed by FastSim:
  //    This is purposefully done this way to allow FastSim to run with OpenHLT: 
  //    The {FastSim+OpenHLT} package runs on the head of HLTrigger/HLTanalyzers 
  //    where there is purposefully this duplication because FastSim does the 
  //    simulation of muons seperately, and needs the same collection. 
  l1extramu_        = conf.getParameter<std::string>   ("l1extramu");
  m_l1extramu       = edm::InputTag(l1extramu_, "");
  
  // read the L1Extra collection name, and add the instance names as needed
  l1extramc_        = conf.getParameter<std::string>   ("l1extramc");
  m_l1extraemi      = edm::InputTag(l1extramc_, "Isolated");
  m_l1extraemn      = edm::InputTag(l1extramc_, "NonIsolated");
  m_l1extrajetc     = edm::InputTag(l1extramc_, "Central");
  m_l1extrajetf     = edm::InputTag(l1extramc_, "Forward");
  m_l1extrataujet   = edm::InputTag(l1extramc_, "Tau");
  m_l1extramet      = edm::InputTag(l1extramc_, "");

  hltresults_       = conf.getParameter<edm::InputTag> ("hltresults");
  gtReadoutRecord_  = conf.getParameter<edm::InputTag> ("l1GtReadoutRecord");
  gtObjectMap_      = conf.getParameter<edm::InputTag> ("l1GtObjectMapRecord");

  // only keep the module label for GCT in 2.X.X - comment from Pedrame:
  //    As far as I (pragmatically) know, this is the way it works up to 2XX series; 
  //    I know that Len is working on making it work in 3XX series.
  gctCounts_        = edm::InputTag( conf.getParameter<edm::InputTag>("l1GctCounts").label(), "" );
  
  MuCandTag2_       = conf.getParameter<edm::InputTag> ("MuCandTag2");
  MuIsolTag2_       = conf.getParameter<edm::InputTag> ("MuIsolTag2");
  MuCandTag3_       = conf.getParameter<edm::InputTag> ("MuCandTag3");
  MuIsolTag3_       = conf.getParameter<edm::InputTag> ("MuIsolTag3");
  MuLinkTag_        = conf.getParameter<edm::InputTag> ("MuLinkTag");
  HLTTau_           = conf.getParameter<edm::InputTag> ("HLTTau");

  // btag OpenHLT input collections
  m_rawBJets                = conf.getParameter<edm::InputTag>("CommonBJetsL2");
  m_correctedBJets          = conf.getParameter<edm::InputTag>("CorrectedBJetsL2");
  m_lifetimeBJetsL25        = conf.getParameter<edm::InputTag>("LifetimeBJetsL25");
  m_lifetimeBJetsL3         = conf.getParameter<edm::InputTag>("LifetimeBJetsL3");
  m_lifetimeBJetsL25Relaxed = conf.getParameter<edm::InputTag>("LifetimeBJetsL25Relaxed");
  m_lifetimeBJetsL3Relaxed  = conf.getParameter<edm::InputTag>("LifetimeBJetsL3Relaxed");
  m_softmuonBJetsL25        = conf.getParameter<edm::InputTag>("SoftmuonBJetsL25");
  m_softmuonBJetsL3         = conf.getParameter<edm::InputTag>("SoftmuonBJetsL3");
  m_performanceBJetsL25     = conf.getParameter<edm::InputTag>("PerformanceBJetsL25");
  m_performanceBJetsL3      = conf.getParameter<edm::InputTag>("PerformanceBJetsL3");

  m_file = 0;   // set to null
  errCnt = 0;

  // read run parameters with a default value 
  edm::ParameterSet runParameters = conf.getParameter<edm::ParameterSet>("RunParameters");
  _HistName = runParameters.getUntrackedParameter<std::string>("HistogramFile", "test.root");
  _EtaMin   = runParameters.getUntrackedParameter<double>("EtaMin", -5.2);
  _EtaMax   = runParameters.getUntrackedParameter<double>("EtaMax",  5.2);

//   cout << "---------- Input Parameters ---------------------------" << endl;
//   cout << "  Output histograms written to: " << _HistName << std::endl;
//   cout << "  EtaMin: " << _EtaMin << endl;    
//   cout << "  EtaMax: " << _EtaMax << endl;    
//   cout << "  Monte:  " << _Monte << endl;    
//   cout << "  Debug:  " << _Debug << endl;    
//   cout << "-------------------------------------------------------" << endl;  

  // open the tree file
  m_file = new TFile(_HistName.c_str(), "RECREATE");
  if (m_file)
    m_file->cd();

  // Initialize the tree
  HltTree = new TTree("HltTree", "");

  // Setup the different analysis
  jet_analysis_.setup(conf, HltTree);
  bjet_analysis_.setup(conf, HltTree);
  elm_analysis_.setup(conf, HltTree);
  muon_analysis_.setup(conf, HltTree);
  mct_analysis_.setup(conf, HltTree);
  hlt_analysis_.setup(conf, HltTree);
  evt_header_.setup(HltTree);
}

// Boiler-plate "analyze" method declaration for an analyzer module.
void HLTAnalyzer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {

  // To get information from the event setup, you must request the "Record"
  // which contains it and then extract the object you need
  //edm::ESHandle<CaloGeometry> geometry;
  //iSetup.get<IdealGeometryRecord>().get(geometry);

  // These declarations create handles to the types of records that you want
  // to retrieve from event "iEvent".
  edm::Handle<CaloJetCollection>                    recjets;
  edm::Handle<GenJetCollection>                     genjets;
  edm::Handle<CaloTowerCollection>                  caloTowers;
  edm::Handle<CaloMETCollection>                    recmet;
  edm::Handle<GenMETCollection>                     genmet;
  edm::Handle<METCollection>                        ht;
  edm::Handle<CandidateView>                        mctruth;
  edm::Handle<double>                               genEventScale;
  edm::Handle<GsfElectronCollection>                electrons;
  edm::Handle<PhotonCollection>                     photons;
  edm::Handle<MuonCollection>                       muon;
  edm::Handle<edm::TriggerResults>                  hltresults;
  edm::Handle<l1extra::L1EmParticleCollection>      l1extemi, l1extemn;
  edm::Handle<l1extra::L1MuonParticleCollection>    l1extmu;
  edm::Handle<l1extra::L1JetParticleCollection>     l1extjetc, l1extjetf, l1exttaujet;
  edm::Handle<l1extra::L1EtMissParticleCollection>  l1extmet;
  edm::Handle<L1GlobalTriggerReadoutRecord>         l1GtRR;
  edm::Handle<L1GlobalTriggerObjectMapRecord>       l1GtOMRec;
  edm::Handle<L1GlobalTriggerObjectMap>             l1GtOM;
  edm::Handle<L1GctJetCountsCollection>             l1GctCounts;
  edm::Handle<RecoChargedCandidateCollection>       mucands2, mucands3;
  edm::Handle<edm::ValueMap<bool> >                 isoMap2,  isoMap3;
  edm::Handle<MuonTrackLinksCollection>             mulinks;
  edm::Handle<reco::HLTTauCollection>               taus;

  // btag OpenHLT input collections
  edm::Handle<edm::View<reco::Jet> >                hRawBJets;
  edm::Handle<edm::View<reco::Jet> >                hCorrectedBJets;
  edm::Handle<reco::JetTagCollection>               hLifetimeBJetsL25;
  edm::Handle<reco::JetTagCollection>               hLifetimeBJetsL3;
  edm::Handle<reco::JetTagCollection>               hLifetimeBJetsL25Relaxed;
  edm::Handle<reco::JetTagCollection>               hLifetimeBJetsL3Relaxed;
  edm::Handle<reco::JetTagCollection>               hSoftmuonBJetsL25;
  edm::Handle<reco::JetTagCollection>               hSoftmuonBJetsL3;
  edm::Handle<reco::JetTagCollection>               hPerformanceBJetsL25;
  edm::Handle<reco::JetTagCollection>               hPerformanceBJetsL3;

  // extract the collections from the event, check their validity and log which are missing
  std::vector<MissingCollectionInfo> missing;

  getCollection( iEvent, missing, recjets,         recjets_,           kRecjets );
  getCollection( iEvent, missing, genjets,         genjets_,           kGenjets );
  getCollection( iEvent, missing, recmet,          recmet_,            kRecmet );
  getCollection( iEvent, missing, genmet,          genmet_,            kGenmet );
  getCollection( iEvent, missing, caloTowers,      calotowers_,        kCaloTowers );
  getCollection( iEvent, missing, ht,              ht_,                kHt );
  getCollection( iEvent, missing, electrons,       Electron_,          kElectrons );
  getCollection( iEvent, missing, photons,         Photon_,            kPhotons );
  getCollection( iEvent, missing, muon,            muon_,              kMuon );
  getCollection( iEvent, missing, taus,            HLTTau_,            kTaus );
  getCollection( iEvent, missing, hltresults,      hltresults_,        kHltresults );
  getCollection( iEvent, missing, l1extemi,        m_l1extraemi,       kL1extemi );
  getCollection( iEvent, missing, l1extemn,        m_l1extraemn,       kL1extemn );
  getCollection( iEvent, missing, l1extmu,         m_l1extramu,        kL1extmu );
  getCollection( iEvent, missing, l1extjetc,       m_l1extrajetc,      kL1extjetc );
  getCollection( iEvent, missing, l1extjetf,       m_l1extrajetf,      kL1extjetf );
  getCollection( iEvent, missing, l1exttaujet,     m_l1extrataujet,    kL1exttaujet );
  getCollection( iEvent, missing, l1extmet,        m_l1extramet,       kL1extmet );
  getCollection( iEvent, missing, l1GtRR,          gtReadoutRecord_,   kL1GtRR );
  getCollection( iEvent, missing, l1GtOMRec,       gtObjectMap_,       kL1GtOMRec );
  getCollection( iEvent, missing, l1GctCounts,     gctCounts_,         kL1GctCounts );
  getCollection( iEvent, missing, mctruth,         mctruth_,           kMctruth );
  getCollection( iEvent, missing, genEventScale,   genEventScale_,     kGenEventScale );
  getCollection( iEvent, missing, mucands2,        MuCandTag2_,        kMucands2 );
  getCollection( iEvent, missing, mucands3,        MuCandTag3_,        kMucands3 );
  getCollection( iEvent, missing, isoMap2,         MuIsolTag2_,        kIsoMap2 );
  getCollection( iEvent, missing, isoMap3,         MuIsolTag3_,        kIsoMap3 );
  getCollection( iEvent, missing, mulinks,         MuLinkTag_,         kMulinks );
  getCollection( iEvent, missing, hRawBJets,                m_rawBJets,                 kBTagJets );
  getCollection( iEvent, missing, hCorrectedBJets,          m_correctedBJets,           kBTagCorrectedJets );
  getCollection( iEvent, missing, hLifetimeBJetsL25,        m_lifetimeBJetsL25,         kBTagLifetimeBJetsL25 );
  getCollection( iEvent, missing, hLifetimeBJetsL3,         m_lifetimeBJetsL3,          kBTagLifetimeBJetsL3 );
  getCollection( iEvent, missing, hLifetimeBJetsL25Relaxed, m_lifetimeBJetsL25Relaxed,  kBTagLifetimeBJetsL25Relaxed );
  getCollection( iEvent, missing, hLifetimeBJetsL3Relaxed,  m_lifetimeBJetsL3Relaxed,   kBTagLifetimeBJetsL3Relaxed );
  getCollection( iEvent, missing, hSoftmuonBJetsL25,        m_softmuonBJetsL25,         kBTagSoftmuonBJetsL25 );
  getCollection( iEvent, missing, hSoftmuonBJetsL3,         m_softmuonBJetsL3,          kBTagSoftmuonBJetsL3 );
  getCollection( iEvent, missing, hPerformanceBJetsL25,     m_performanceBJetsL25,      kBTagPerformanceBJetsL25 );
  getCollection( iEvent, missing, hPerformanceBJetsL3,      m_performanceBJetsL3,       kBTagPerformanceBJetsL3 );

  // print missing collections
  if (not missing.empty() and (errCnt < errMax())) {
    errCnt++;
    std::stringstream out;       
    out <<  "OpenHLT analyser - missing collections:";
    BOOST_FOREACH(const MissingCollectionInfo & entry, missing)
      out << "\n\t" << entry.first << ": " << entry.second->encode();
    edm::LogPrint("OpenHLT") << out.str() << std::endl; 
    if (errCnt == errMax())
      edm::LogWarning("OpenHLT") << "Maximum error count reached -- No more messages will be printed.";
  }

  // run the analysis, passing required event fragments
  jet_analysis_.analyze(
    recjets.product(),
    genjets.product(),
    recmet.product(),
    genmet.product(),
    ht.product(),
    taus.product(),
    caloTowers.product(),
    HltTree);
  
  muon_analysis_.analyze(
    muon.product(),
    mucands2.product(),
    isoMap2.product(),
    mucands3.product(),
    isoMap3.product(),
    mulinks.product(),
    HltTree);
  
  elm_analysis_.analyze(iEvent, iSetup, 
    electrons.product(),
    photons.product(),
    HltTree);
  
  mct_analysis_.analyze(
    mctruth.product(),
    genEventScale.product(),
    HltTree);
  
  hlt_analysis_.analyze(
    hltresults.product(),
    l1extemi.product(),
    l1extemn.product(),
    l1extmu.product(),
    l1extjetc.product(),
    l1extjetf.product(),
    l1exttaujet.product(),
    l1extmet.product(),
    l1GtRR.product(),
    l1GtOMRec.product(),
    l1GctCounts.product(),
    HltTree);
  
  bjet_analysis_.analyze(
    hRawBJets.product(), 
    hCorrectedBJets.product(),
    hLifetimeBJetsL25.product(),
    hLifetimeBJetsL3.product(),
    hLifetimeBJetsL25Relaxed.product(),
    hLifetimeBJetsL3Relaxed.product(),
    hSoftmuonBJetsL25.product(),
    hSoftmuonBJetsL3.product(),
    hPerformanceBJetsL25.product(),
    hPerformanceBJetsL3.product(),
    HltTree);

  evt_header_.analyze(iEvent, HltTree);

  // std::cout << " Ending Event Analysis" << std::endl;
  // After analysis, fill the variables tree
  if (m_file)
    m_file->cd();
  HltTree->Fill();
}

// "endJob" is an inherited method that you may implement to do post-EOF processing and produce final output.
void HLTAnalyzer::endJob() {

  if (m_file)
    m_file->cd();

  HltTree->Write();
  delete HltTree;
  HltTree = 0;

  if (m_file) {         // if there was a tree file...
    m_file->Write();    // write out the branches
    delete m_file;      // close and delete the file
    m_file = 0;         // set to zero to clean up
  }

}
