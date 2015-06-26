// system include files
#include <memory>
#include <vector>

#include "JetDifferentialEfficiencies.h"
std::vector<float> StoreJetPt(const edm::Event& iEvent, edm::EDGetTokenT<std::vector<pat::Jet>> jets_){

  edm::Handle < std::vector<pat::Jet> > jets;
  iEvent.getByToken(jets_, jets);
  std::vector<float> pT;
  pT.clear();

  for( auto jet = jets->begin(); jet != jets->end(); ++jet ){ 
    if( jet->pt()<30. || std::abs(jet->eta())>3.0 ){
      continue; // skip jets with low pT or outside the tracker acceptance
    }
    pT.push_back(jet->pt());
  }

  return pT;
}