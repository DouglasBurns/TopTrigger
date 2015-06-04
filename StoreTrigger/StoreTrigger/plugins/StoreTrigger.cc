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
// Original Author:  Douglas Burns
//         Created:  Thu, 28 May 2015 15:07:20 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <fstream>

// user include files
#include "StoreTrigger.h"


StoreTrigger::StoreTrigger(const edm::ParameterSet& iConfig) :
    hltInputTag_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTInputTag"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects")))

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

  // std::cout << "In analyze" << std::endl;
  edm::Handle < edm::TriggerResults > triggerResults;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  iEvent.getByToken(hltInputTag_, triggerResults);
  iEvent.getByToken(triggerObjects_, triggerObjects);

  const edm::TriggerNames &TrigNames = iEvent.triggerNames(*triggerResults); 

  pathNames.clear();
  passTrig.clear();
  pathIndices.clear();

  std::ifstream ListofTriggers("ListofTriggers.txt");
  if ( ListofTriggers.is_open() ){
    while ( !ListofTriggers.eof() ){
      ListofTriggers >> trigger;
      std::cout << trigger << std::endl;
      pathNames.push_back(trigger);
    }
    ListofTriggers.close();
  }

  // trigger = "HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1";
  // pathNames.push_back(trigger);
  // trigger = "HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v1";
  // pathNames.push_back(trigger);
  // trigger = "HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1";
  // pathNames.push_back(trigger);
  // trigger = "HLT_Ele27_eta2p1_WP75_Gsf_v1";
  // pathNames.push_back(trigger);

  // trigger = "HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1";
  // pathNames.push_back(trigger);
  // trigger = "HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet30_v1";
  // pathNames.push_back(trigger);
  // trigger = "HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1";
  // pathNames.push_back(trigger);
  // trigger = "HLT_Ele32_eta2p1_WP75_Gsf_v1";
  // pathNames.push_back(trigger);

  // trigger = "HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v1";
  // pathNames.push_back(trigger);
  // trigger = "HLT_IsoMu20_eta2p1_TriCentralPFJet30_v1";
  // pathNames.push_back(trigger);
  // trigger = "HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v1";
  // pathNames.push_back(trigger);
  // trigger = "HLT_IsoMu20_eta2p1_v1";
  // pathNames.push_back(trigger);

  // trigger = "HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1";
  // pathNames.push_back(trigger);
  // trigger = "HLT_IsoMu24_eta2p1_TriCentralPFJet30_v1";
  // pathNames.push_back(trigger);
  // trigger = "HLT_IsoMu24_eta2p1_TriCentralPFJet50_40_30_v1";
  // pathNames.push_back(trigger);
  // trigger = "HLT_IsoMu24_eta2p1_v1";
  // pathNames.push_back(trigger);

  std::cout << "Number of triggers interested in : " << pathNames.size() << std::endl;
  std::cout << "Number of trigger results : " << triggerResults->size() << std::endl;
  // Note : if the path doesn't exist, the index returned will be triggerResults->size()
  // This will then cause a seg fault when you do triggerResults->accept()

  // Print out names of all triggers in this event
  // for (unsigned int i = 0; i < triggerResults->size(); ++i) {
  //   std::cout << "Path " << i << " name : " << TrigNames.triggerName( i ) << std::endl;
  // }



  for (unsigned int j = 0; j < pathNames.size(); ++j) {

    // Check if it passes
    TrigDecision=triggerResults->accept(TrigNames.triggerIndex(pathNames.at(j)));
    passTrig.push_back(TrigDecision);
    pathIndices.push_back(j);

  }


  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
    obj.unpackPathNames(TrigNames);
    // std::cout << "\t Trigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
    // Print trigger object collection and type
    // std::cout << "\t Collection: " << obj.collection() << std::endl;
  }






   outTree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
StoreTrigger::beginJob()
{
  outTree = fileService->make<TTree>("Triggers", "TriggerDecisions");

  outTree->Branch("Trigger_Decision", "std::vector<int>", &passTrig);
  outTree->Branch("Name_of_Trigger", "std::vector<std::string>", &pathNames);
  outTree->Branch("Trigger_Index", "std::vector<int>", &pathIndices);

  // Normal Ouput = outTree->Branch("new_v", &new_v, "new_v/F"); Look at https://root.cern.ch/root/html/TTree.html
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
