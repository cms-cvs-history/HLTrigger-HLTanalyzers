#ifndef HLTrigger_HLTanalyzers_HLTBJet_h
#define HLTrigger_HLTanalyzers_HLTBJet_h

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

class edm::ParameterSet;
class TTree;

class HLTBJet {
public:
  HLTBJet();
  ~HLTBJet();
  
  void setup(const edm::ParameterSet & conf, TTree * tree);
  void analyze(
      const edm::View<reco::Jet> & lifetimeBjetL2, const reco::JetTagCollection & lifetimeBjetL25, const reco::JetTagCollection & lifetimeBjetL3, 
      const edm::View<reco::Jet> & softmuonBjetL2, const reco::JetTagCollection & softmuonBjetL25, const reco::JetTagCollection & softmuonBjetL3, 
      const edm::View<reco::Jet> & performanceBjetL2, const reco::JetTagCollection & performanceBjetL25, const reco::JetTagCollection & performanceBjetL3, 
      const reco::METCollection & ht, 
      TTree * tree);

private:
  // set of variables for lifetime-based b-tag
  int m_lifetimeBJets;
  float * m_lifetimeBJetL2Energy;
  float * m_lifetimeBJetL2ET;
  float * m_lifetimeBJetL2Eta;
  float * m_lifetimeBJetL2Phi;
  float * m_lifetimeBJetL25Energy;
  float * m_lifetimeBJetL25ET;
  float * m_lifetimeBJetL25Eta;
  float * m_lifetimeBJetL25Phi;
  float * m_lifetimeBJetL25Discriminator;
  float * m_lifetimeBJetL3Energy;
  float * m_lifetimeBJetL3ET;
  float * m_lifetimeBJetL3Eta;
  float * m_lifetimeBJetL3Phi;
  float * m_lifetimeBJetL3Discriminator;
  
  // set of variables for soft-muon-based b-tag
  int m_softmuonBJets;
  float * m_softmuonBJetL2Energy;
  float * m_softmuonBJetL2ET;
  float * m_softmuonBJetL2Eta;
  float * m_softmuonBJetL2Phi;
  float * m_softmuonBJetL25Energy;
  float * m_softmuonBJetL25ET;
  float * m_softmuonBJetL25Eta;
  float * m_softmuonBJetL25Phi;
  float * m_softmuonBJetL25Discriminator;       // this is actually a boolean value - do not optimize
  float * m_softmuonBJetL3Energy;
  float * m_softmuonBJetL3ET;
  float * m_softmuonBJetL3Eta;
  float * m_softmuonBJetL3Phi;
  float * m_softmuonBJetL3Discriminator;
  
  // set of variables for b-tagging performance measurements
  int m_performanceBJets;
  float * m_performanceBJetL2Energy;
  float * m_performanceBJetL2ET;
  float * m_performanceBJetL2Eta;
  float * m_performanceBJetL2Phi;
  float * m_performanceBJetL25Energy;
  float * m_performanceBJetL25ET;
  float * m_performanceBJetL25Eta;
  float * m_performanceBJetL25Phi;
  float * m_performanceBJetL25Discriminator;   // this is actually a boolean value - do not optimize 
  float * m_performanceBJetL3Energy;
  float * m_performanceBJetL3ET;
  float * m_performanceBJetL3Eta;
  float * m_performanceBJetL3Phi;
  float * m_performanceBJetL3Discriminator;     // this is actually a boolean value - do not optimize
};

#endif // HLTrigger_HLTanalyzers_HLTBJet_h
