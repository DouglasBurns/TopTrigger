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
    ttbarjet304050trigger_(iConfig.getParameter <std::string> ("TTBarJet304050TriggerInput")),
    asymmetricjet30filter_(iConfig.getParameter <std::string> ("AsymmetricJet30FilterInput")),
    asymmetricjet40filter_(iConfig.getParameter <std::string> ("AsymmetricJet40FilterInput")),
    asymmetricjet50filter_(iConfig.getParameter <std::string> ("AsymmetricJet50FilterInput")){
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
      }
      if ( TrigNames.triggerName(i).find(ttbarjet304050trigger_) != std::string::npos ) {
        ttbarjet304050Index = i;
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


  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
    for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
      if ( obj.filterLabels()[i].find(asymmetricjet30filter_) != std::string::npos ) {
        TTBarJet30Hist_Pt->Fill(obj.pt());
        TTBarJet30Hist_Eta->Fill(obj.eta());
        TTBarJet30Hist_Phi->Fill(obj.phi());
      }
    }
      for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
      if ( obj.filterLabels()[i].find(asymmetricjet40filter_) != std::string::npos ) {
        TTBarJet40Hist_Pt->Fill(obj.pt());
        TTBarJet40Hist_Eta->Fill(obj.eta());
        TTBarJet40Hist_Phi->Fill(obj.phi());
      }
    }    for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
      if ( obj.filterLabels()[i].find(asymmetricjet50filter_) != std::string::npos ) {
        TTBarJet50Hist_Pt->Fill(obj.pt());
        TTBarJet50Hist_Eta->Fill(obj.eta());
        TTBarJet50Hist_Phi->Fill(obj.phi());
      }
    }
  }
}
// ------------ method called once each job just before starting event loop  ------------
void 
TTBarJet304050::beginJob(){

  subDir_TrigDec = fileService->mkdir( "Trigger Decision" );
  TTBarJet304050Hist = subDir_TrigDec.make<TH1D>("TTBarJet304050 Trigger Decision", ttbarjet304050trigger_.c_str(), 2, -0.5, 1.5);
  SingleLeptonHist = subDir_TrigDec.make<TH1D>("Single Lepton Trigger Decision", singleleptontrigger_.c_str(), 2, -0.5, 1.5);

  subDir_AsymmetricJet30Filter = fileService->mkdir( "Asymmetric30JetFilter" );
  TTBarJet30Hist_Pt = subDir_AsymmetricJet30Filter.make<TH1D>("Pt", "Pt", 100, 0, 300);
  TTBarJet30Hist_Eta = subDir_AsymmetricJet30Filter.make<TH1D>("Eta", "Eta", 100, -5, 5);
  TTBarJet30Hist_Phi = subDir_AsymmetricJet30Filter.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);

  subDir_AsymmetricJet40Filter = fileService->mkdir( "Asymmetric40JetFilter" );
  TTBarJet40Hist_Pt = subDir_AsymmetricJet40Filter.make<TH1D>("Pt", "Pt", 100, 0, 300);
  TTBarJet40Hist_Eta = subDir_AsymmetricJet40Filter.make<TH1D>("Eta", "Eta", 100, -5, 5);
  TTBarJet40Hist_Phi = subDir_AsymmetricJet40Filter.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);

  subDir_AsymmetricJet50Filter = fileService->mkdir( "Asymmetric50JetFilter" );
  TTBarJet50Hist_Pt = subDir_AsymmetricJet50Filter.make<TH1D>("Pt", "Pt", 100, 0, 300);
  TTBarJet50Hist_Eta = subDir_AsymmetricJet50Filter.make<TH1D>("Eta", "Eta", 100, -5, 5);
  TTBarJet50Hist_Phi = subDir_AsymmetricJet50Filter.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);
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
