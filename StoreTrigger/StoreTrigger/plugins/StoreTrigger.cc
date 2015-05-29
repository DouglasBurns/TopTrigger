// -*- C++ -*-
//
// Package:    StoreTrigger/StoreTrigger
// Class:      StoreTrigger
// 
/**\class StoreTrigger StoreTrigger.cc StoreTrigger/StoreTrigger/plugins/StoreTrigger.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emyr Clement
//         Created:  Thu, 28 May 2015 15:07:20 GMT
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "StoreTrigger.h"


StoreTrigger::StoreTrigger(const edm::ParameterSet& iConfig) :
    hltInputTag_(iConfig.getParameter < edm::InputTag > ("HLTInputTag"))
{
   //now do what ever initialization is needed

}


StoreTrigger::~StoreTrigger()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
StoreTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  std::cout << "In analyze" << std::endl;
  edm::Handle < edm::TriggerResults > triggerResults;
  iEvent.getByLabel(hltInputTag_, triggerResults);
  const edm::TriggerNames& TrigNames = iEvent.triggerNames(*triggerResults); 

  pathNames.clear();
  passTrig.clear();

  trigger = "HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1";//WPLoose
  pathNames.push_back(trigger);
  trigger = "HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v1";//WPLoose
  pathNames.push_back(trigger);
  trigger = "HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1";//WPLoose
  pathNames.push_back(trigger);
  trigger = "HLT_Ele27_eta2p1_WP75_Gsf_v1";//WPLoose
  pathNames.push_back(trigger);
  // trigger = "HLT_Ele27_eta2p1_WPTight_Gsf_v1";//WPTight
  // pathNames.push_back(trigger);
  trigger = "HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1";//WPLoose
  pathNames.push_back(trigger);
  trigger = "HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet30_v1";//WPLoose
  pathNames.push_back(trigger);
  trigger = "HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1";//WPLoose
  pathNames.push_back(trigger);
  trigger = "HLT_Ele32_eta2p1_WP75_Gsf_v1";//WPLoose
  pathNames.push_back(trigger);
  // trigger = "HLT_Ele32_eta2p1_WPTight_Gsf_v1";//WPTight
  // pathNames.push_back(trigger);
  trigger = "HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v1";//v2
  pathNames.push_back(trigger);
  trigger = "HLT_IsoMu20_eta2p1_TriCentralPFJet30_v1";//v2
  pathNames.push_back(trigger);
  trigger = "HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v1";//v2
  pathNames.push_back(trigger);
  trigger = "HLT_IsoMu20_eta2p1_v1";//v2
  pathNames.push_back(trigger);
  // trigger = "HLT_IsoMu20_v1";//v2
  // pathNames.push_back(trigger);
  trigger = "HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1";//v2
  pathNames.push_back(trigger);
  trigger = "HLT_IsoMu24_eta2p1_TriCentralPFJet30_v1";//v2
  pathNames.push_back(trigger);
  trigger = "HLT_IsoMu24_eta2p1_TriCentralPFJet50_40_30_v1";//v2
  pathNames.push_back(trigger);
  trigger = "HLT_IsoMu24_eta2p1_v1";//v2
  pathNames.push_back(trigger);

  std::cout << "Number of triggers interested in : " << pathNames.size() << std::endl;
  std::cout << "Number of trigger results : " << triggerResults->size() << std::endl;
  // Note : if the path doesn't exist, the index returned will be triggerResults->size()
  // This will then cause a seg fault when you do triggerResults->accept()

  // Print out names of all triggers in this event
  // for (unsigned int i = 0; i < triggerResults->size(); ++i) {
  //   std::cout << "Path " << i << " name : " << TrigNames.triggerName( i ) << std::endl;
  // }

  for (unsigned int j = 0; j < pathNames.size(); ++j) {

    std::cout << j << "th Trig index : " << TrigNames.triggerIndex(pathNames.at(j)) << std::endl;

    // Check if it passes
    bool TrigDecision=triggerResults->accept(TrigNames.triggerIndex(pathNames.at(j)));
    passTrig.push_back(TrigDecision);
    std::cout << "pass trigger : " << passTrig.at(j) << std::endl;

  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
StoreTrigger::beginJob()
{
  outTree = fileService->make<TTree>("Triggers", "TriggerDecisions");
  outTree->Branch("Name_of_Trigger", &pathNames);
  outTree->Branch("Trigger_Decision", &passTrig);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
StoreTrigger::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
StoreTrigger::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
StoreTrigger::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
StoreTrigger::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
StoreTrigger::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
StoreTrigger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(StoreTrigger);
