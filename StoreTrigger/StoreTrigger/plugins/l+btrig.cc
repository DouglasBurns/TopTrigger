// -*- C++ -*-
//
// Package:    StoreTrigger/StoreTrigger
// Class:      SingleTop
// 
/**\class SingleTop SingleTop.cc StoreTrigger/StoreTrigger/plugins/SingleTop.cc

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
#include "l+btrig.h"


SingleTop::SingleTop(const edm::ParameterSet& iConfig) :
    hltInputTag_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTInputTag"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("HLTriggerObjects"))),
    singleleptontrigger_(iConfig.getParameter <std::string> ("SingleLeptonTriggerInput")),
    singletoptrigger_(iConfig.getParameter <std::string> ("SingleTopTriggerInput")),
    singletopfilter_(iConfig.getParameter <std::string> ("SingleTopFilterInput")),
    bjetfilter_(iConfig.getParameter <std::string> ("BJetFilterInput")){
   //now do what ever initialization is needed
}


SingleTop::~SingleTop(){
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SingleTop::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

    using namespace edm;

    // std::cout << "In analyze" << std::endl;
    edm::Handle < edm::TriggerResults > triggerResults;
    edm::Handle < pat::TriggerObjectStandAloneCollection > triggerObjects;

    iEvent.getByToken(hltInputTag_, triggerResults);
    iEvent.getByToken(triggerObjects_, triggerObjects);

    const edm::TriggerNames &TrigNames = iEvent.triggerNames(*triggerResults); 

    for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) {
        if ( TrigNames.triggerName(i).find(singleleptontrigger_) != std::string::npos ) {
          singleleptonIndex = i;
        }
        if ( TrigNames.triggerName(i).find(singletoptrigger_) != std::string::npos ) {
          singletopIndex = i;
        }
    }


    if ( singleleptonIndex < triggerResults->size() ) {
      SingleLeptonTrigDecision = triggerResults->accept(singleleptonIndex);
      SingleLeptonHist->Fill(SingleLeptonTrigDecision);
    }
    else {
      std::cout << "Looking for : " << singleleptontrigger_ << " but failed" << std::endl;
    }


    if ( singletopIndex < triggerResults->size() ) {
      SingleTopTrigDecision = triggerResults->accept(singletopIndex);
      SingleTopHist->Fill(SingleTopTrigDecision);
    }
    else {
      std::cout << "Looking for : " << singletoptrigger_ << " but failed" << std::endl;
    }

    // std::cout << "\n === TRIGGER OBJECTS === " << std::endl;

    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
        // std::cout << "Trigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;

        for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
          if ( obj.filterLabels()[i].find(singletopfilter_) != std::string::npos ) {
            SingleTopFilterHist_Pt->Fill(obj.pt());
            SingleTopFilterHist_Eta->Fill(obj.eta());
            SingleTopFilterHist_Phi->Fill(obj.phi());
          }
        }

        for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
          if ( obj.filterLabels()[i].find(bjetfilter_) != std::string::npos ) {
            BJetFilterHist_Pt->Fill(obj.pt());
            BJetFilterHist_Eta->Fill(obj.eta());
            BJetFilterHist_Phi->Fill(obj.phi());
          }
        }


        // // Print associated trigger filters
        // std::cout << "Filters:    " << std::endl;
        // for (unsigned h = 0; h < obj.filterLabels().size(); ++h){

        //     std::cout << obj.filterLabels()[h] << std::endl;
        // } 
    }
}

// ------------ method called once each job just before starting event loop  ------------
void 
SingleTop::beginJob(){
  subDir_TrigDec = fileService->mkdir( "Trigger Decision" );
  SingleTopHist = subDir_TrigDec.make<TH1D>("Single Top Trigger Decision", singletoptrigger_.c_str(), 2, -0.5, 1.5);
  SingleLeptonHist = subDir_TrigDec.make<TH1D>("Single Lepton Trigger Decision", singleleptontrigger_.c_str(), 2, -0.5, 1.5);


  subDir_SingleTopFilter = fileService->mkdir( "SingleTopFilter" );
  SingleTopFilterHist_Pt = subDir_SingleTopFilter.make<TH1D>("Pt", "Pt", 100, 0, 300);
  SingleTopFilterHist_Eta = subDir_SingleTopFilter.make<TH1D>("Eta", "Eta", 100, -5, 5);
  SingleTopFilterHist_Phi = subDir_SingleTopFilter.make<TH1D>("Phi", "Phi", 100, -3.5, 3);

  subDir_BJetFilter = fileService->mkdir( "BJetFilter" );
  BJetFilterHist_Pt = subDir_BJetFilter.make<TH1D>("Pt", "Pt", 100, 0, 300);
  BJetFilterHist_Eta = subDir_BJetFilter.make<TH1D>("Eta", "Eta", 100, -5, 5);
  BJetFilterHist_Phi = subDir_BJetFilter.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);


}

// ------------ method called once each job just after ending the event loop  ------------
void 
SingleTop::endJob(){

}

// ------------ method called when starting to processes a run  ------------
/*
void 
SingleTop::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
SingleTop::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
SingleTop::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
SingleTop::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SingleTop::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SingleTop);
