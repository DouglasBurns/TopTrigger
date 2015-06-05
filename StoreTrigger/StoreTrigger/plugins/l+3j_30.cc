// -*- C++ -*-
//
// Package:    StoreTrigger/StoreTrigger
// Class:      TTBarJet30
// 
/**\class TTBarJet30 TTBarJet30.cc StoreTrigger/StoreTrigger/plugins/TTBarJet30.cc

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
#include "l+3j_30.h"


TTBarJet30::TTBarJet30(const edm::ParameterSet& iConfig) :
    hltInputTag_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTInputTag"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("HLTriggerObjects"))),
    singleleptontrigger_(iConfig.getParameter <std::string> ("SingleLeptonTriggerInput")),
    ttbarjet30trigger_(iConfig.getParameter <std::string> ("TTBarJet30TriggerInput")),
    symmetricjetfilter_(iConfig.getParameter <std::string> ("SymmetricJetFilterInput")){
   //now do what ever initialization is needed

}


TTBarJet30::~TTBarJet30(){
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TTBarJet30::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

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
      if ( TrigNames.triggerName(i).find(ttbarjet30trigger_) != std::string::npos ) {
          ttbarjet30Index = i;
      }
  }


  if ( singleleptonIndex < triggerResults->size() ) {
    SingleLeptonTrigDecision = triggerResults->accept(singleleptonIndex);
    SingleLeptonHist->Fill(SingleLeptonTrigDecision);
  }
  else {
    std::cout << "Looking for : " << singleleptontrigger_ << " but failed" << std::endl;
  }


  if ( ttbarjet30Index < triggerResults->size() ) {
    TTBarJet30TrigDecision = triggerResults->accept(ttbarjet30Index);
    TTBarJet30Hist->Fill(TTBarJet30TrigDecision);
  }
  else {
    std::cout << "Looking for : " << ttbarjet30trigger_ << " but failed" << std::endl;
  }



  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
      if ( obj.filterLabels()[i].find(symmetricjetfilter_) != std::string::npos ) {
        TTBarJet30Hist_Pt->Fill(obj.pt());
        TTBarJet30Hist_Eta->Fill(obj.eta());
        TTBarJet30Hist_Phi->Fill(obj.phi());
      }
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
TTBarJet30::beginJob(){

  subDir_TrigDec = fileService->mkdir( "Trigger Decision" );
  TTBarJet30Hist = subDir_TrigDec.make<TH1D>("TTBarJet30 Trigger Decision", ttbarjet30trigger_.c_str(), 2, -0.5, 1.5);
  SingleLeptonHist = subDir_TrigDec.make<TH1D>("Single Lepton Trigger Decision", singleleptontrigger_.c_str(), 2, -0.5, 1.5);

  subDir_SymmetricJetFilter = fileService->mkdir( "SymmetricJetFilter" );
  TTBarJet30Hist_Pt = subDir_SymmetricJetFilter.make<TH1D>("Pt", "Pt", 100, 0, 300);
  TTBarJet30Hist_Eta = subDir_SymmetricJetFilter.make<TH1D>("Eta", "Eta", 100, -5, 5);
  TTBarJet30Hist_Phi = subDir_SymmetricJetFilter.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TTBarJet30::endJob(){

}

// ------------ method called when starting to processes a run  ------------
/*
void 
TTBarJet30::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
TTBarJet30::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TTBarJet30::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TTBarJet30::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTBarJet30::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TTBarJet30);
