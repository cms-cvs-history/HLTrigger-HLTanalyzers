#include <cmath>

#include <TTree.h>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h" 
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"                   // this is needed here due to IMO a bug somewhere else
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "HLTrigger/HLTanalyzers/interface/HLTBJet.h"

const int kMaxBJets  = 4;
const int kMaxTracks = 100;

HLTBJet::HLTBJet() 
{
  initLifetime();
  initSoftmuon();
  initPerformance();
}

void HLTBJet::initLifetime() 
{
  m_lifetimeBJets = 0;
  m_lifetimeBJetL2Energy         = new float[kMaxBJets];
  m_lifetimeBJetL2ET             = new float[kMaxBJets];
  m_lifetimeBJetL2Eta            = new float[kMaxBJets];
  m_lifetimeBJetL2Phi            = new float[kMaxBJets];
  m_lifetimeBJetL25Discriminator = new float[kMaxBJets];
  m_lifetimeBJetL3Discriminator  = new float[kMaxBJets];
  m_pixelTracks = 0;
  m_pixelTrackPt                 = new float[kMaxTracks];
  m_pixelTrackEta                = new float[kMaxTracks];
  m_pixelTrackPhi                = new float[kMaxTracks];
  m_pixelTrackChi2               = new float[kMaxTracks];
  m_regionalTracks = 0;
  m_regionalTrackPt              = new float[kMaxTracks];
  m_regionalTrackEta             = new float[kMaxTracks];
  m_regionalTrackPhi             = new float[kMaxTracks];
  m_regionalTrackChi2            = new float[kMaxTracks];
  m_regionalSeedPt               = new float[kMaxTracks];
  m_regionalSeedEta              = new float[kMaxTracks];
  m_regionalSeedPhi              = new float[kMaxTracks];
}

void HLTBJet::initSoftmuon() 
{
  m_softmuonBJets = 0;
  m_softmuonBJetL2Energy         = new float[kMaxBJets];
  m_softmuonBJetL2ET             = new float[kMaxBJets];
  m_softmuonBJetL2Eta            = new float[kMaxBJets];
  m_softmuonBJetL2Phi            = new float[kMaxBJets];
  m_softmuonBJetL25Discriminator = new float[kMaxBJets];
  m_softmuonBJetL3Discriminator  = new float[kMaxBJets];
}

void HLTBJet::initPerformance() 
{
  m_performanceBJets = 0;
  m_performanceBJetL2Energy         = new float[kMaxBJets];
  m_performanceBJetL2ET             = new float[kMaxBJets];
  m_performanceBJetL2Eta            = new float[kMaxBJets];
  m_performanceBJetL2Phi            = new float[kMaxBJets];
  m_performanceBJetL25Discriminator = new float[kMaxBJets];
  m_performanceBJetL3Discriminator  = new float[kMaxBJets];
}

HLTBJet::~HLTBJet() 
{ }

void HLTBJet::setup(const edm::ParameterSet & conf, TTree * tree)
{
  setupLifetime(conf, tree);
  setupSoftmuon(conf, tree);
  setupPerformance(conf, tree);
}

void HLTBJet::setupLifetime(const edm::ParameterSet & conf, TTree * tree)
{ 
  tree->Branch("NohBJetLife",                       & m_lifetimeBJets,                  "NohBJetLife/I");
  tree->Branch("ohBJetLifeL2E",                     m_lifetimeBJetL2Energy,             "ohBJetLifeL2E[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL2ET",                    m_lifetimeBJetL2ET,                 "ohBJetLifeL2ET[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL2Eta",                   m_lifetimeBJetL2Eta,                "ohBJetLifeL2Eta[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL2Phi",                   m_lifetimeBJetL2Phi,                "ohBJetLifeL2Phi[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL25Discriminator",        m_lifetimeBJetL25Discriminator,     "ohBJetLifeL25Discriminator[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL3Discriminator",         m_lifetimeBJetL3Discriminator,      "ohBJetLifeL3Discriminator[NohBJetLife]/F");
  tree->Branch("NohBJetPixelTracks",                & m_pixelTracks,                    "NohBJetPixelTracks/I");
  tree->Branch("ohBJetLifePixelTrackPt",            m_pixelTrackPt,                     "ohBJetLifePixelTrackPt[NohBJetPixelTracks]/F");
  tree->Branch("ohBJetLifePixelTrackEta",           m_pixelTrackEta,                    "ohBJetLifePixelTrackEta[NohBJetPixelTracks]/F");
  tree->Branch("ohBJetLifePixelTrackPhi",           m_pixelTrackPhi,                    "ohBJetLifePixelTrackPhi[NohBJetPixelTracks]/F");
  tree->Branch("ohBJetLifePixelTrackChi2",          m_pixelTrackChi2,                   "ohBJetLifePixelTrackChi2[NohBJetPixelTracks]/F");
  tree->Branch("NohBJetRegionalTracks",             & m_regionalTracks,                 "NohBJetRegionalTracks/I");
  tree->Branch("ohBJetLifeRegionalTrackPt",         m_regionalTrackPt,                  "ohBJetLifeRegionalTrackPt[NohBJetRegionalTracks]/F");
  tree->Branch("ohBJetLifeRegionalTrackEta",        m_regionalTrackEta,                 "ohBJetLifeRegionalTrackEta[NohBJetRegionalTracks]/F");
  tree->Branch("ohBJetLifeRegionalTrackPhi",        m_regionalTrackPhi,                 "ohBJetLifeRegionalTrackPhi[NohBJetRegionalTracks]/F");
  tree->Branch("ohBJetLifeRegionalTrackChi2",       m_regionalTrackChi2,                "ohBJetLifeRegionalTrackChi2[NohBJetRegionalTracks]/F");
  tree->Branch("ohBJetLifeRegionalSeedPt",          m_regionalSeedPt,                   "ohBJetLifeRegionalSeedPt[NohBJetRegionalTracks]/F");
  tree->Branch("ohBJetLifeRegionalSeedEta",         m_regionalSeedEta,                  "ohBJetLifeRegionalSeedEta[NohBJetRegionalTracks]/F");
  tree->Branch("ohBJetLifeRegionalSeedPhi",         m_regionalSeedPhi,                  "ohBJetLifeRegionalSeedPhi[NohBJetRegionalTracks]/F");
}

void HLTBJet::setupSoftmuon(const edm::ParameterSet & conf, TTree * tree)
{ 
  tree->Branch("NohBJetSoftm",                      & m_softmuonBJets,                  "NohBJetSoftm/I");
  tree->Branch("ohBJetSoftmL2E",                    m_softmuonBJetL2Energy,             "ohBJetSoftmL2E[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL2ET",                   m_softmuonBJetL2ET,                 "ohBJetSoftmL2ET[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL2Eta",                  m_softmuonBJetL2Eta,                "ohBJetSoftmL2Eta[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL2Phi",                  m_softmuonBJetL2Phi,                "ohBJetSoftmL2Phi[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL25Discriminator",       m_softmuonBJetL25Discriminator,     "ohBJetSoftmL25Discriminator[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL3Discriminator",        m_softmuonBJetL3Discriminator,      "ohBJetSoftmL3Discriminator[NohBJetSoftm]/F");
}

void HLTBJet::setupPerformance(const edm::ParameterSet & conf, TTree * tree)
{ 
  tree->Branch("NohBJetPerf",                       & m_performanceBJets,               "NohBJetPerf/I");
  tree->Branch("ohBJetPerfL2E",                     m_performanceBJetL2Energy,          "ohBJetPerfL2E[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL2ET",                    m_performanceBJetL2ET,              "ohBJetPerfL2ET[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL2Eta",                   m_performanceBJetL2Eta,             "ohBJetPerfL2Eta[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL2Phi",                   m_performanceBJetL2Phi,             "ohBJetPerfL2Phi[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL25Discriminator",        m_performanceBJetL25Discriminator,  "ohBJetPerfL25Discriminator[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL3Discriminator",         m_performanceBJetL3Discriminator,   "ohBJetPerfL3Discriminator[NohBJetPerf]/F");
}

void HLTBJet::update(const edm::EventSetup & setup) 
{
  edm::ESHandle<GlobalTrackingGeometry> h_geometry;
  setup.get<GlobalTrackingGeometryRecord>().get(h_geometry);
  m_geometry = h_geometry.product();

  edm::ESHandle<MagneticField> h_field;
  setup.get<IdealMagneticFieldRecord>().get(h_field);
  m_field = h_field.product();
}

void HLTBJet::analyzeLifetime(
    const edm::View<reco::Jet>   & lifetimeBjetL2, 
    const reco::JetTagCollection & lifetimeBjetL25, 
    const reco::JetTagCollection & lifetimeBjetL3, 
    const reco::TrackCollection  & lifetimePixelTracks, 
    const reco::TrackCollection  & lifetimeRegionalTracks)
{
  m_lifetimeBJets = lifetimeBjetL2.size();
  if (m_lifetimeBJets > kMaxBJets) m_lifetimeBJets = kMaxBJets;
  // no filter is applied, so all collections *should* have the same number of elements
  for (int i = 0; i < m_lifetimeBJets; i++) {
    m_lifetimeBJetL2Energy[i]         = lifetimeBjetL2[i].energy();
    m_lifetimeBJetL2ET[i]             = lifetimeBjetL2[i].et();
    m_lifetimeBJetL2Eta[i]            = lifetimeBjetL2[i].eta();
    m_lifetimeBJetL2Phi[i]            = lifetimeBjetL2[i].phi();
    m_lifetimeBJetL25Discriminator[i] = lifetimeBjetL25[i].second;
    m_lifetimeBJetL3Discriminator[i]  = lifetimeBjetL3[i].second;
  }
  m_pixelTracks = lifetimePixelTracks.size();
  if (m_pixelTracks > kMaxTracks) m_pixelTracks = kMaxTracks;
  for (int i = 0; i < m_pixelTracks; ++i) {
    m_pixelTrackPt[i]   = lifetimePixelTracks[i].pt();
    m_pixelTrackEta[i]  = lifetimePixelTracks[i].eta();
    m_pixelTrackPhi[i]  = lifetimePixelTracks[i].phi();
    m_pixelTrackChi2[i] = lifetimePixelTracks[i].chi2();
  }
  m_regionalTracks = lifetimeRegionalTracks.size();
  if (m_regionalTracks > kMaxTracks) m_regionalTracks = kMaxTracks;
  for (int i = 0; i < m_regionalTracks; ++i) {
    GlobalVector momentum = seedMomentum(lifetimeRegionalTracks[i]);    // (0, 0, 0) if TrackExtra or the seed ar not available
    m_regionalTrackPt[i]   = lifetimeRegionalTracks[i].pt();
    m_regionalTrackEta[i]  = lifetimeRegionalTracks[i].eta();
    m_regionalTrackPhi[i]  = lifetimeRegionalTracks[i].phi();
    m_regionalTrackChi2[i] = lifetimeRegionalTracks[i].chi2();
    m_regionalSeedPt[i]    = momentum.perp();
    m_regionalSeedEta[i]   = momentum.eta();
    m_regionalSeedPhi[i]   = momentum.phi();
  }

}

void HLTBJet::analyzeSoftmuon(
    const edm::View<reco::Jet>   & softmuonBjetL2, 
    const reco::JetTagCollection & softmuonBjetL25, 
    const reco::JetTagCollection & softmuonBjetL3)
{
  m_softmuonBJets = softmuonBjetL2.size();
  if (m_softmuonBJets > kMaxBJets) m_softmuonBJets = kMaxBJets;
  for (int i = 0; i < m_softmuonBJets; i++) {
    m_softmuonBJetL2Energy[i]         = softmuonBjetL2[i].energy();
    m_softmuonBJetL2ET[i]             = softmuonBjetL2[i].et();
    m_softmuonBJetL2Eta[i]            = softmuonBjetL2[i].eta();
    m_softmuonBJetL2Phi[i]            = softmuonBjetL2[i].phi();
    m_softmuonBJetL25Discriminator[i] = softmuonBjetL25[i].second;
    m_softmuonBJetL3Discriminator[i]  = softmuonBjetL3[i].second;
  }
}

void HLTBJet::analyzePerformance(
    const edm::View<reco::Jet>   & performanceBjetL2, 
    const reco::JetTagCollection & performanceBjetL25, 
    const reco::JetTagCollection & performanceBjetL3)
{
  m_performanceBJets = performanceBjetL2.size();
  if (m_performanceBJets > kMaxBJets) m_performanceBJets = kMaxBJets;
  for (int i = 0; i < m_performanceBJets; i++) {
    m_performanceBJetL2Energy[i]         = performanceBjetL2[i].energy();
    m_performanceBJetL2ET[i]             = performanceBjetL2[i].et();
    m_performanceBJetL2Eta[i]            = performanceBjetL2[i].eta();
    m_performanceBJetL2Phi[i]            = performanceBjetL2[i].phi();
    m_performanceBJetL25Discriminator[i] = performanceBjetL25[i].second;
    m_performanceBJetL3Discriminator[i]  = performanceBjetL3[i].second;
  }
}

GlobalVector HLTBJet::seedMomentum(const reco::Track & track) const
{
  static TrajectoryStateTransform transform;
  
  if (track.extra().isNull() or track.extra()->seedRef().isNull())
    return GlobalVector();
  else {
    const PTrajectoryStateOnDet & tsod = track.seedRef()->startingState();
    const Surface * surface = & m_geometry->idToDet( tsod.detId() )->surface();
    return transform.transientState(tsod, surface, m_field).globalMomentum();
  }
}
  
