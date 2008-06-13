#ifndef HLTrigger_HLTanalyzers_HLTBJet_h
#define HLTrigger_HLTanalyzers_HLTBJet_h

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

class TTree;
class edm::ParameterSet;
class edm::EventSetup;
class MagneticField;
class GlobalTrackingGeometry;
class TrajectorySeed;

class HLTBJet {
public:
  HLTBJet();
  ~HLTBJet();
  
  void setup(const edm::ParameterSet & conf, TTree * tree);
  void update(const edm::EventSetup & setup);

  void analyzeLifetime(
      const edm::View<reco::Jet>   & lifetimeBjetL2, 
      const reco::JetTagCollection & lifetimeBjetL25, 
      const reco::JetTagCollection & lifetimeBjetL3, 
      const reco::TrackCollection  & lifetimePixelTracks, 
      const reco::TrackCollection  & lifetimeRegionalTracks);

  void analyzeSoftmuon(
      const edm::View<reco::Jet>   & softmuonBjetL2, 
      const reco::JetTagCollection & softmuonBjetL25, 
      const reco::JetTagCollection & softmuonBjetL3);

  void analyzePerformance(
      const edm::View<reco::Jet>   & performanceBjetL2, 
      const reco::JetTagCollection & performanceBjetL25, 
      const reco::JetTagCollection & performanceBjetL3);

private:
  void initLifetime();
  void initSoftmuon();
  void initPerformance();

  void setupLifetime(const edm::ParameterSet & conf, TTree * tree);
  void setupSoftmuon(const edm::ParameterSet & conf, TTree * tree);
  void setupPerformance(const edm::ParameterSet & conf, TTree * tree);
  
  GlobalVector seedMomentum(const reco::Track & track) const;

  // MagneticField and GlobalTrackingGeometry are needed to convert from local to global coordinates
  const MagneticField *          m_field;
  const GlobalTrackingGeometry * m_geometry;
  
  // set of variables for lifetime-based b-tag
  int m_lifetimeBJets;
  float * m_lifetimeBJetL2Energy;
  float * m_lifetimeBJetL2ET;
  float * m_lifetimeBJetL2Eta;
  float * m_lifetimeBJetL2Phi;
  float * m_lifetimeBJetL25Discriminator;
  float * m_lifetimeBJetL3Discriminator;
  int m_pixelTracks;
  float * m_pixelTrackPt;
  float * m_pixelTrackEta;
  float * m_pixelTrackPhi;
  float * m_pixelTrackChi2;
  int m_regionalTracks;
  float * m_regionalTrackPt;
  float * m_regionalTrackEta;
  float * m_regionalTrackPhi;
  float * m_regionalTrackChi2;
  float * m_regionalSeedPt;
  float * m_regionalSeedEta;
  float * m_regionalSeedPhi;

  // set of variables for soft-muon-based b-tag
  int m_softmuonBJets;
  float * m_softmuonBJetL2Energy;
  float * m_softmuonBJetL2ET;
  float * m_softmuonBJetL2Eta;
  float * m_softmuonBJetL2Phi;
  float * m_softmuonBJetL25Discriminator;       // this is actually a boolean value - do not optimize
  float * m_softmuonBJetL3Discriminator;
  
  // set of variables for b-tagging performance measurements
  int m_performanceBJets;
  float * m_performanceBJetL2Energy;
  float * m_performanceBJetL2ET;
  float * m_performanceBJetL2Eta;
  float * m_performanceBJetL2Phi;
  float * m_performanceBJetL25Discriminator;    // this is actually a boolean value - do not optimize 
  float * m_performanceBJetL3Discriminator;     // this is actually a boolean value - do not optimize
};

#endif // HLTrigger_HLTanalyzers_HLTBJet_h
