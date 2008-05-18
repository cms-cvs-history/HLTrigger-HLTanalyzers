#include <cmath>

#include <TTree.h>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "HLTrigger/HLTanalyzers/interface/HLTBJet.h"

const int kMaxBJets = 4;
  
HLTBJet::HLTBJet() 
{
  m_lifetimeBJets = 0;
  m_lifetimeBJetL2Energy         = new float[kMaxBJets];
  m_lifetimeBJetL2ET             = new float[kMaxBJets];
  m_lifetimeBJetL2Eta            = new float[kMaxBJets];
  m_lifetimeBJetL2Phi            = new float[kMaxBJets];
  m_lifetimeBJetL25Energy        = new float[kMaxBJets];
  m_lifetimeBJetL25ET            = new float[kMaxBJets];
  m_lifetimeBJetL25Eta           = new float[kMaxBJets];
  m_lifetimeBJetL25Phi           = new float[kMaxBJets];
  m_lifetimeBJetL25Discriminator = new float[kMaxBJets];
  m_lifetimeBJetL3Energy         = new float[kMaxBJets];
  m_lifetimeBJetL3ET             = new float[kMaxBJets];
  m_lifetimeBJetL3Eta            = new float[kMaxBJets];
  m_lifetimeBJetL3Phi            = new float[kMaxBJets];
  m_lifetimeBJetL3Discriminator  = new float[kMaxBJets];
  m_softmuonBJets = 0;
  m_softmuonBJetL2Energy         = new float[kMaxBJets];
  m_softmuonBJetL2ET             = new float[kMaxBJets];
  m_softmuonBJetL2Eta            = new float[kMaxBJets];
  m_softmuonBJetL2Phi            = new float[kMaxBJets];
  m_softmuonBJetL25Energy        = new float[kMaxBJets];
  m_softmuonBJetL25ET            = new float[kMaxBJets];
  m_softmuonBJetL25Eta           = new float[kMaxBJets];
  m_softmuonBJetL25Phi           = new float[kMaxBJets];
  m_softmuonBJetL25Discriminator = new float[kMaxBJets];
  m_softmuonBJetL3Energy         = new float[kMaxBJets];
  m_softmuonBJetL3ET             = new float[kMaxBJets];
  m_softmuonBJetL3Eta            = new float[kMaxBJets];
  m_softmuonBJetL3Phi            = new float[kMaxBJets];
  m_softmuonBJetL3Discriminator  = new float[kMaxBJets];
  m_performanceBJets = 0;
  m_performanceBJetL2Energy         = new float[kMaxBJets];
  m_performanceBJetL2ET             = new float[kMaxBJets];
  m_performanceBJetL2Eta            = new float[kMaxBJets];
  m_performanceBJetL2Phi            = new float[kMaxBJets];
  m_performanceBJetL25Energy        = new float[kMaxBJets];
  m_performanceBJetL25ET            = new float[kMaxBJets];
  m_performanceBJetL25Eta           = new float[kMaxBJets];
  m_performanceBJetL25Phi           = new float[kMaxBJets];
  m_performanceBJetL25Discriminator = new float[kMaxBJets];
  m_performanceBJetL3Energy         = new float[kMaxBJets];
  m_performanceBJetL3ET             = new float[kMaxBJets];
  m_performanceBJetL3Eta            = new float[kMaxBJets];
  m_performanceBJetL3Phi            = new float[kMaxBJets];
  m_performanceBJetL3Discriminator  = new float[kMaxBJets];
}

HLTBJet::~HLTBJet() 
{ }

void HLTBJet::setup(const edm::ParameterSet & conf, TTree * tree)
{ 
  tree->Branch("NohBJetLife",                       & m_lifetimeBJets,                  "NohBJetLife/I");
  tree->Branch("ohBJetLifeL2E",                     m_lifetimeBJetL2Energy,             "ohBJetLifeL2E[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL2ET",                    m_lifetimeBJetL2ET,                 "ohBJetLifeL2ET[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL2Eta",                   m_lifetimeBJetL2Eta,                "ohBJetLifeL2Eta[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL2Phi",                   m_lifetimeBJetL2Phi,                "ohBJetLifeL2Phi[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL25E",                    m_lifetimeBJetL25Energy,            "ohBJetLifeL25E[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL25ET",                   m_lifetimeBJetL25ET,                "ohBJetLifeL25ET[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL25Eta",                  m_lifetimeBJetL25Eta,               "ohBJetLifeL25Eta[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL25Phi",                  m_lifetimeBJetL25Phi,               "ohBJetLifeL25Phi[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL25Discriminator",        m_lifetimeBJetL25Discriminator,     "ohBJetLifeL25Discriminator[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL3E",                     m_lifetimeBJetL3Energy,             "ohBJetLifeL3E[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL3ET",                    m_lifetimeBJetL3ET,                 "ohBJetLifeL3ET[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL3Eta",                   m_lifetimeBJetL3Eta,                "ohBJetLifeL3Eta[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL3Phi",                   m_lifetimeBJetL3Phi,                "ohBJetLifeL3Phi[NohBJetLife]/F");
  tree->Branch("ohBJetLifeL3Discriminator",         m_lifetimeBJetL3Discriminator,      "ohBJetLifeL3Discriminator[NohBJetLife]/F");
  tree->Branch("NohBJetSoftm",                      & m_softmuonBJets,                  "NohBJetSoftm/I");
  tree->Branch("ohBJetSoftmL2E",                    m_softmuonBJetL2Energy,             "ohBJetSoftmL2E[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL2ET",                   m_softmuonBJetL2ET,                 "ohBJetSoftmL2ET[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL2Eta",                  m_softmuonBJetL2Eta,                "ohBJetSoftmL2Eta[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL2Phi",                  m_softmuonBJetL2Phi,                "ohBJetSoftmL2Phi[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL25E",                   m_softmuonBJetL25Energy,            "ohBJetSoftmL25E[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL25ET",                  m_softmuonBJetL25ET,                "ohBJetSoftmL25ET[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL25Eta",                 m_softmuonBJetL25Eta,               "ohBJetSoftmL25Eta[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL25Phi",                 m_softmuonBJetL25Phi,               "ohBJetSoftmL25Phi[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL25Discriminator",       m_softmuonBJetL25Discriminator,     "ohBJetSoftmL25Discriminator[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL3E",                    m_softmuonBJetL3Energy,             "ohBJetSoftmL3E[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL3ET",                   m_softmuonBJetL3ET,                 "ohBJetSoftmL3ET[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL3Eta",                  m_softmuonBJetL3Eta,                "ohBJetSoftmL3Eta[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL3Phi",                  m_softmuonBJetL3Phi,                "ohBJetSoftmL3Phi[NohBJetSoftm]/F");
  tree->Branch("ohBJetSoftmL3Discriminator",        m_softmuonBJetL3Discriminator,      "ohBJetSoftmL3Discriminator[NohBJetSoftm]/F");
  tree->Branch("NohBJetPerf",                       & m_performanceBJets,               "NohBJetPerf/I");
  tree->Branch("ohBJetPerfL2E",                     m_performanceBJetL2Energy,          "ohBJetPerfL2E[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL2ET",                    m_performanceBJetL2ET,              "ohBJetPerfL2ET[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL2Eta",                   m_performanceBJetL2Eta,             "ohBJetPerfL2Eta[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL2Phi",                   m_performanceBJetL2Phi,             "ohBJetPerfL2Phi[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL25E",                    m_performanceBJetL25Energy,         "ohBJetPerfL25E[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL25ET",                   m_performanceBJetL25ET,             "ohBJetPerfL25ET[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL25Eta",                  m_performanceBJetL25Eta,            "ohBJetPerfL25Eta[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL25Phi",                  m_performanceBJetL25Phi,            "ohBJetPerfL25Phi[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL25Discriminator",        m_performanceBJetL25Discriminator,  "ohBJetPerfL25Discriminator[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL3E",                     m_performanceBJetL3Energy,          "ohBJetPerfL3E[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL3ET",                    m_performanceBJetL3ET,              "ohBJetPerfL3ET[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL3Eta",                   m_performanceBJetL3Eta,             "ohBJetPerfL3Eta[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL3Phi",                   m_performanceBJetL3Phi,             "ohBJetPerfL3Phi[NohBJetPerf]/F");
  tree->Branch("ohBJetPerfL3Discriminator",         m_performanceBJetL3Discriminator,   "ohBJetPerfL3Discriminator[NohBJetPerf]/F");
}

void HLTBJet::analyze(
    const edm::View<reco::Jet> & lifetimeBjetL2, const reco::JetTagCollection & lifetimeBjetL25, const reco::JetTagCollection & lifetimeBjetL3, 
    const edm::View<reco::Jet> & softmuonBjetL2, const reco::JetTagCollection & softmuonBjetL25, const reco::JetTagCollection & softmuonBjetL3, 
    const edm::View<reco::Jet> & performanceBjetL2, const reco::JetTagCollection & performanceBjetL25, const reco::JetTagCollection & performanceBjetL3, 
    const reco::METCollection & ht, 
    TTree * tree) 
{
  m_lifetimeBJets = lifetimeBjetL2.size();
  if (m_lifetimeBJets > kMaxBJets) m_lifetimeBJets = kMaxBJets;
  // no filter is applied, so all collections *should* have the same number of elements
  for (int i = 0; i < m_lifetimeBJets; i++) {
    std::cerr << '\t' << i << " L2.5 jet: " << lifetimeBjetL25[i].first.isNonnull() << " L3 jet: " << lifetimeBjetL3[i].first.isNonnull() << std::endl;
    m_lifetimeBJetL2Energy[i]         = lifetimeBjetL2[i].energy();
    m_lifetimeBJetL2ET[i]             = lifetimeBjetL2[i].et();
    m_lifetimeBJetL2Eta[i]            = lifetimeBjetL2[i].eta();
    m_lifetimeBJetL2Phi[i]            = lifetimeBjetL2[i].phi();
    m_lifetimeBJetL25Energy[i]        = lifetimeBjetL25[i].first->energy();
    m_lifetimeBJetL25ET[i]            = lifetimeBjetL25[i].first->et();
    m_lifetimeBJetL25Eta[i]           = lifetimeBjetL25[i].first->eta();
    m_lifetimeBJetL25Phi[i]           = lifetimeBjetL25[i].first->phi();
    m_lifetimeBJetL25Discriminator[i] = lifetimeBjetL25[i].second;
    m_lifetimeBJetL3Energy[i]         = lifetimeBjetL3[i].first->energy();
    m_lifetimeBJetL3ET[i]             = lifetimeBjetL3[i].first->et();
    m_lifetimeBJetL3Eta[i]            = lifetimeBjetL3[i].first->eta();
    m_lifetimeBJetL3Phi[i]            = lifetimeBjetL3[i].first->phi();
    m_lifetimeBJetL3Discriminator[i]  = lifetimeBjetL3[i].second;
  }

  m_softmuonBJets = softmuonBjetL2.size();
  if (m_softmuonBJets > kMaxBJets) m_softmuonBJets = kMaxBJets;
  // no filter is applied, so all collections *should* have the same number of elements
  for (int i = 0; i < m_softmuonBJets; i++) {
    std::cerr << '\t' << i << " L2.5 jet: " << softmuonBjetL25[i].first.isNonnull() << " L3 jet: " << softmuonBjetL3[i].first.isNonnull() << std::endl;
    m_softmuonBJetL2Energy[i]         = softmuonBjetL2[i].energy();
    m_softmuonBJetL2ET[i]             = softmuonBjetL2[i].et();
    m_softmuonBJetL2Eta[i]            = softmuonBjetL2[i].eta();
    m_softmuonBJetL2Phi[i]            = softmuonBjetL2[i].phi();
    m_softmuonBJetL25Energy[i]        = softmuonBjetL25[i].first->energy();
    m_softmuonBJetL25ET[i]            = softmuonBjetL25[i].first->et();
    m_softmuonBJetL25Eta[i]           = softmuonBjetL25[i].first->eta();
    m_softmuonBJetL25Phi[i]           = softmuonBjetL25[i].first->phi();
    m_softmuonBJetL25Discriminator[i] = softmuonBjetL25[i].second;
    m_softmuonBJetL3Energy[i]         = softmuonBjetL3[i].first->energy();
    m_softmuonBJetL3ET[i]             = softmuonBjetL3[i].first->et();
    m_softmuonBJetL3Eta[i]            = softmuonBjetL3[i].first->eta();
    m_softmuonBJetL3Phi[i]            = softmuonBjetL3[i].first->phi();
    m_softmuonBJetL3Discriminator[i]  = softmuonBjetL3[i].second;
  }

  m_performanceBJets = performanceBjetL2.size();
  if (m_performanceBJets > kMaxBJets) m_performanceBJets = kMaxBJets;
  // no filter is applied, so all collections *should* have the same number of elements
  for (int i = 0; i < m_performanceBJets; i++) {
    std::cerr << '\t' << i << " L2.5 jet: " << performanceBjetL25[i].first.isNonnull() << " L3 jet: " << performanceBjetL3[i].first.isNonnull() << std::endl;
    m_performanceBJetL2Energy[i]         = performanceBjetL2[i].energy();
    m_performanceBJetL2ET[i]             = performanceBjetL2[i].et();
    m_performanceBJetL2Eta[i]            = performanceBjetL2[i].eta();
    m_performanceBJetL2Phi[i]            = performanceBjetL2[i].phi();
    m_performanceBJetL25Energy[i]        = performanceBjetL25[i].first->energy();
    m_performanceBJetL25ET[i]            = performanceBjetL25[i].first->et();
    m_performanceBJetL25Eta[i]           = performanceBjetL25[i].first->eta();
    m_performanceBJetL25Phi[i]           = performanceBjetL25[i].first->phi();
    m_performanceBJetL25Discriminator[i] = performanceBjetL25[i].second;
    m_performanceBJetL3Energy[i]         = performanceBjetL3[i].first->energy();
    m_performanceBJetL3ET[i]             = performanceBjetL3[i].first->et();
    m_performanceBJetL3Eta[i]            = performanceBjetL3[i].first->eta();
    m_performanceBJetL3Phi[i]            = performanceBjetL3[i].first->phi();
    m_performanceBJetL3Discriminator[i]  = performanceBjetL3[i].second;
  }
}
