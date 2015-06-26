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
    jets_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    singleleptontrigger_(iConfig.getParameter <std::string> ("SingleLeptonTriggerInput")),
    singletoptrigger_(iConfig.getParameter <std::string> ("SingleTopTriggerInput")),
    singletopfilter_(iConfig.getParameter <std::string> ("SingleTopFilterInput")),
    bjetfilter_(iConfig.getParameter <std::string> ("BJetFilterInput")){
}

SingleTop::~SingleTop(){
}

void
SingleTop::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

   using namespace edm;

    // Used to create output string for hadronic output | single lepton
    CombinedTrigger.append(singleleptontrigger_.c_str());
    CombinedTrigger.append(" and ");
    CombinedTrigger.append(singletoptrigger_.c_str());

    // Handle to the trigger decisions and variables 
    edm::Handle < edm::TriggerResults > triggerResults;
    edm::Handle < pat::TriggerObjectStandAloneCollection > triggerObjects;

    // Get particular information required
    iEvent.getByToken(hltInputTag_, triggerResults);
    iEvent.getByToken(triggerObjects_, triggerObjects);

    // List of trigger names
    const edm::TriggerNames &TrigNames = iEvent.triggerNames(*triggerResults); 

    // See if trigger name you inputted exists, if so note its index
    for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) {
        if ( TrigNames.triggerName(i).find(singleleptontrigger_) != std::string::npos ) {
          singleleptonIndex = i;
        }
        if ( TrigNames.triggerName(i).find(singletoptrigger_) != std::string::npos ) {
          singletopIndex = i;
        }
    }

    // If single lepton trigger name exists, find outcome and store in histogram
    if ( singleleptonIndex < triggerResults->size() ) {
      SingleLeptonTrigDecision = triggerResults->accept(singleleptonIndex);
      SingleLeptonHist->Fill(SingleLeptonTrigDecision);

      // If single lepton trigger outcome is true, see if hadronic trigger exists (for added eff)
      if (SingleLeptonTrigDecision == true){

        // If hadronic trigger name exists, find outcome and store in histogram
        if ( singletopIndex < triggerResults->size() ) {
          SingleTopCombinedTrigDecision = triggerResults->accept(singletopIndex);
          SingleTopCombinedHist->Fill(SingleTopCombinedTrigDecision);
        }

        // Cant find hadronic trigger
        else {
          std::cout << "Looking for : " << singletoptrigger_ << " but failed" << std::endl;
        }

      }

    }

    // Cant find single lepton trigger
    else {
      std::cout << "Looking for : " << singleleptontrigger_ << " but failed" << std::endl;
    }

    // If hadronic trigger name exists, find outcome and store in histogram
    if ( singletopIndex < triggerResults->size() ) {
      SingleTopTrigDecision = triggerResults->accept(singletopIndex);
      SingleTopHist->Fill(SingleTopTrigDecision);
    }

    // Cant find hadronic trigger    
    else {
      std::cout << "Looking for : " << singletoptrigger_ << " but failed" << std::endl;
    }

    StoreJets(iEvent, SingleTopTrigDecision);


    // std::cout << "\n === TRIGGER OBJECTS === " << std::endl;

    // Get list of trigger objects: jets, leptons etc
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
        // std::cout << "Trigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;

        // Look for specific filter, if it exists store trigger object variables in histogram
        for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
          if ( obj.filterLabels()[i].find(singletopfilter_) != std::string::npos ) {
            SingleTopFilterHist_Pt->Fill(obj.pt());
            SingleTopFilterHist_Eta->Fill(obj.eta());
            SingleTopFilterHist_Phi->Fill(obj.phi());
          }
        }

        // Look for specific filter, if it exists store trigger object variables in histogram
        for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
          if ( obj.filterLabels()[i].find(bjetfilter_) != std::string::npos ) {
            BJetFilterHist_Pt->Fill(obj.pt());
            BJetFilterHist_Eta->Fill(obj.eta());
            BJetFilterHist_Phi->Fill(obj.phi());
          }
        }


        // // Print trigger filters associated with each triggerObject
        // std::cout << "Filters:    " << std::endl;
        // for (unsigned h = 0; h < obj.filterLabels().size(); ++h){

        //     std::cout << obj.filterLabels()[h] << std::endl;
        // } 
    }
}

// ------------ method called once each job just before starting event loop  ------------
void 
SingleTop::beginJob(){
  // create subdirectories with histograms
  subDir_TrigDec = fileService->mkdir( "Trigger Decision" );
  SingleTopHist = subDir_TrigDec.make<TH1D>("SingleTop Trigger Decision", singletoptrigger_.c_str(), 2, -0.5, 1.5);
  SingleLeptonHist = subDir_TrigDec.make<TH1D>("Single Lepton Trigger Decision", singleleptontrigger_.c_str(), 2, -0.5, 1.5);
  SingleTopCombinedHist = subDir_TrigDec.make<TH1D>("Added SingleTop Trigger Decision", CombinedTrigger.c_str(), 2, -0.5, 1.5);

  subDir_SingleTopFilter = fileService->mkdir( "SingleTopFilter" );
  SingleTopFilterHist_Pt = subDir_SingleTopFilter.make<TH1D>("Pt", "Pt", 100, 0, 300);
  SingleTopFilterHist_Eta = subDir_SingleTopFilter.make<TH1D>("Eta", "Eta", 100, -5, 5);
  SingleTopFilterHist_Phi = subDir_SingleTopFilter.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);

  subDir_BJetFilter = fileService->mkdir( "BJetFilter" );
  BJetFilterHist_Pt = subDir_BJetFilter.make<TH1D>("Pt", "Pt", 100, 0, 300);
  BJetFilterHist_Eta = subDir_BJetFilter.make<TH1D>("Eta", "Eta", 100, -5, 5);
  BJetFilterHist_Phi = subDir_BJetFilter.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);

  subDir_TrigDiffEff = fileService->mkdir( "TrigDiffErr" );
  SingleTop_Pass_JetPtHist = subDir_TrigDiffEff.make<TH1D>("SingleTop_Pass_JetPt", "TriggerPass_Pt", 100, 0, 300);
  SingleTop_Pass_JetEtaHist = subDir_TrigDiffEff.make<TH1D>("SingleTop_Pass_JetEta", "TriggerPass_Eta", 100, -3, 3);
  SingleTop_Total_JetPtHist = subDir_TrigDiffEff.make<TH1D>("SingleTop_Total_JetPt", "Total_Pt", 100, 0, 300);
  SingleTop_Total_JetEtaHist = subDir_TrigDiffEff.make<TH1D>("SingleTop_Total_JetEta", "Total_Eta", 100, -3, 3);
}


void
SingleTop::StoreJets(const edm::Event& iEvent, bool TrigDecision){
  edm::Handle < std::vector<pat::Jet> > jets;
  iEvent.getByToken(jets_, jets);

  for( auto jet = jets->begin(); jet != jets->end(); ++jet ){ 
    if( jet->pt()<30. || std::abs(jet->eta())>3.0 ){
      continue; // skip jets with low pT or outside the tracker acceptance
    }

    pT = jet->pt();
    eta = jet->eta();

    if(TrigDecision==true){
      SingleTop_Pass_JetPtHist->Fill(pT);
      SingleTop_Pass_JetEtaHist->Fill(eta);
    }
    SingleTop_Total_JetPtHist->Fill(pT);
    SingleTop_Total_JetEtaHist->Fill(eta);
  }
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
