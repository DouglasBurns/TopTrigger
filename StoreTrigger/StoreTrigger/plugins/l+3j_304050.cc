// -*- C++ -*-
//
// Package:    StoreTrigger/StoreTrigger
// Class:      TTBarJet304050
// 
/**\class TTBarJet304050 TTBarJet304050.cc StoreTrigger/StoreTrigger/plugins/TTBarJet304050.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Douglas Burns
//         Created:  Thu, 28 May 2015 15:07:20 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <fstream>

// user include files
#include "l+3j_304050.h"


TTBarJet304050::TTBarJet304050(const edm::ParameterSet& iConfig) :
    hltInputTag_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTInputTag"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("HLTriggerObjects"))),
    singleleptontrigger_(iConfig.getParameter <std::string> ("SingleLeptonTriggerInput")),
    ttbarjet304050trigger_(iConfig.getParameter <std::string> ("TTBarJet304050TriggerInput")){
   //now do what ever initialization is needed

}


TTBarJet304050::~TTBarJet304050(){
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TTBarJet304050::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;

  // std::cout << "In analyze" << std::endl;
  edm::Handle < edm::TriggerResults > triggerResults;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  iEvent.getByToken(hltInputTag_, triggerResults);
  iEvent.getByToken(triggerObjects_, triggerObjects);

  const edm::TriggerNames &TrigNames = iEvent.triggerNames(*triggerResults); 

  for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) {
      if ( TrigNames.triggerName(i).find(singleleptontrigger_) != std::string::npos ) {
        singleleptonIndex = i;
        singleleptontrigger = TrigNames.triggerName(i);
      }
      if ( TrigNames.triggerName(i).find(ttbarjet304050trigger_) != std::string::npos ) {
        ttbarjet304050Index = i;
        ttbarjet304050trigger = TrigNames.triggerName(i);
      }
  }


  if ( singleleptonIndex < triggerResults->size() ) {
    SingleLeptonTrigDecision = triggerResults->accept(singleleptonIndex);
    SingleLeptonHist->Fill(SingleLeptonTrigDecision);
  }
  else {
    std::cout << "Looking for : " << singleleptontrigger_ << " but failed" << std::endl;
  }


  if ( ttbarjet304050Index < triggerResults->size() ) {
    TTBarJet304050TrigDecision = triggerResults->accept(ttbarjet304050Index);
    // if (TTBarJet304050TrigDecision == true){
    TTBarJet304050Hist->Fill(TTBarJet304050TrigDecision);
  // }
  }
  else {
    std::cout << "Looking for : " << ttbarjet304050trigger_ << " but failed" << std::endl;
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
TTBarJet304050::beginJob(){
  SingleLeptonHist = fileService->make<TH1D>("Single Lepton Trigger Decision", "Single Lepton Trigger Decision", 2, -0.5, 1.5);
  TTBarJet304050Hist = fileService->make<TH1D>("TTBarJet304050 Trigger Decision", "TTBarJet304050 Trigger Decision", 2, -0.5, 1.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TTBarJet304050::endJob(){

}

// ------------ method called when starting to processes a run  ------------
/*
void 
TTBarJet304050::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
TTBarJet304050::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TTBarJet304050::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TTBarJet304050::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTBarJet304050::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TTBarJet304050);
