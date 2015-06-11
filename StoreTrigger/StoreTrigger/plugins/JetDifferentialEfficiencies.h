
// user include files
// #include "JetDifferentialEfficiencies.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

std::vector<float> StoreJetPt(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<pat::Jet>> jets_);
