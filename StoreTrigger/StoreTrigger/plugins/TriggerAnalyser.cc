// -*- C++ -*-
//
// Package:    StoreTrigger/StoreTrigger
// Class:      TriggerAnalyser
// 
/**\class TriggerAnalyser TriggerAnalyser.cc StoreTrigger/StoreTrigger/plugins/TriggerAnalyser.cc

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
#include "TriggerAnalyser.h"


TriggerAnalyser::TriggerAnalyser(const edm::ParameterSet& iConfig) :
    hltInputTag_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTInputTag"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("HLTriggerObjects"))),
    jets_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    mets_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
    singleleptontrigger_(iConfig.getParameter <std::string> ("SingleLeptonTriggerInput")),
    crosstrigger_(iConfig.getParameter <std::string> ("CrossTriggerInput")),
    filter1_(iConfig.getParameter <std::string> ("FilterInput1")),
    filter2_(iConfig.getParameter <std::string> ("FilterInput2")),
    filter3_(iConfig.getParameter <std::string> ("FilterInput3")),
    hadronicleg_(iConfig.getParameter <std::string> ("HadronicLeg")),
    leptonicleg_(iConfig.getParameter <std::string> ("LeptonicLeg")){
   //now do what ever initialization is needed

}


TriggerAnalyser::~TriggerAnalyser(){
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

      CombinedTrigger.append(singleleptontrigger_.c_str());
      CombinedTrigger.append(" and ");
      CombinedTrigger.append(crosstrigger_.c_str());

      // std::cout << "Lepton Leg : " << leptonicleg_ << std::endl;

      using namespace edm;

      // std::cout << "In analyze" << std::endl;
      edm::Handle < edm::TriggerResults > triggerResults;
      edm::Handle < pat::TriggerObjectStandAloneCollection > triggerObjects;
      edm::Handle < std::vector<pat::Jet> > jets;
      edm::Handle < std::vector<pat::MET> > mets;

      iEvent.getByToken(hltInputTag_, triggerResults);
      iEvent.getByToken(triggerObjects_, triggerObjects);
      iEvent.getByToken(jets_, jets);
      iEvent.getByToken(mets_, mets);

      const edm::TriggerNames &TrigNames = iEvent.triggerNames(*triggerResults); 

      for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) {
            if ( TrigNames.triggerName(i).find(singleleptontrigger_) != std::string::npos ) {
                  singleleptonIndex = i;
            }
            if ( TrigNames.triggerName(i).find(crosstrigger_) != std::string::npos ) {
                  crossIndex = i;
            }
      }


      if ( singleleptonIndex < triggerResults->size() ) {
            SingleLeptonTrigDecision = triggerResults->accept(singleleptonIndex);
            SingleLeptonHist->Fill(SingleLeptonTrigDecision);

            if (SingleLeptonTrigDecision == true){

                  if ( crossIndex < triggerResults->size() ) {
                        CrossTriggerCombinedTrigDecision = triggerResults->accept(crossIndex);
                        CrossTriggerCombinedHist->Fill(CrossTriggerCombinedTrigDecision);
                  }

                  else {
                        std::cout << "Looking for : " << crosstrigger_ << " but failed" << std::endl;
                  }
            }
      }

      else {
            std::cout << "Looking for : " << singleleptontrigger_ << " but failed" << std::endl;
      }


      if ( crossIndex < triggerResults->size() ) {
            CrossTriggerTrigDecision = triggerResults->accept(crossIndex);
            // if (CrossTriggerTrigDecision == true){
            CrossTriggerHist->Fill(CrossTriggerTrigDecision);
            // }
      }
      else {
            std::cout << "Looking for : " << crosstrigger_ << " but failed" << std::endl;
      }



      // Select jets and store distributions used for differential efficiencies
      for( auto jet = jets->begin(); jet != jets->end(); ++jet ){ 
            // skip jets with low pT or outside the tracker acceptance
            if( jet->pt()<30. || std::abs(jet->eta())>3.0 ){
                  continue; 
            }

            jetPt = jet->pt();
            jetEta = jet->eta();

            if(CrossTriggerTrigDecision==true){
                  CrossTrigger_Pass_JetPtHist->Fill(jetPt);
                  CrossTrigger_Pass_JetEtaHist->Fill(jetEta);
            }
            CrossTrigger_Total_JetPtHist->Fill(jetPt);
            CrossTrigger_Total_JetEtaHist->Fill(jetEta);
      }








      for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 

            if ( hadronicleg_ == "SingleTop" ){
                  for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
                        if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {
                                Filter1_Pt->Fill(obj.pt());
                                Filter1_Eta->Fill(obj.eta());
                                Filter1_Phi->Fill(obj.phi());
                        }
                  }
                  for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
                        if ( obj.filterLabels()[i].find(filter2_) != std::string::npos ) {
                                Filter2_Pt->Fill(obj.pt());
                                Filter2_Eta->Fill(obj.eta());
                                Filter2_Phi->Fill(obj.phi());
                        }
                  }  
            }
            
            if ( hadronicleg_ == "TTBarJet30" ){
                  for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
                        if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {
                                Filter1_Pt->Fill(obj.pt());
                                Filter1_Eta->Fill(obj.eta());
                                Filter1_Phi->Fill(obj.phi());
                        }
                  }
            }

            if ( hadronicleg_ == "TTBarJet304050" ){
                  for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
                        if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {
                                Filter1_Pt->Fill(obj.pt());
                                Filter1_Eta->Fill(obj.eta());
                                Filter1_Phi->Fill(obj.phi());
                        }
                  }
                  for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
                        if ( obj.filterLabels()[i].find(filter2_) != std::string::npos ) {
                                Filter2_Pt->Fill(obj.pt());
                                Filter2_Eta->Fill(obj.eta());
                                Filter2_Phi->Fill(obj.phi());
                        }
                  }    
                  for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {
                        if ( obj.filterLabels()[i].find(filter3_) != std::string::npos ) {
                                Filter3_Pt->Fill(obj.pt());
                                Filter3_Eta->Fill(obj.eta());
                                Filter3_Phi->Fill(obj.phi());
                        }
                  }
            }
      }
}










// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyser::beginJob(){

      subDir_TrigDec = fileService->mkdir( "Trigger Decision" );
      CrossTriggerHist = subDir_TrigDec.make<TH1D>("CrossTrigger Trigger Decision", crosstrigger_.c_str(), 2, -0.5, 1.5);
      SingleLeptonHist = subDir_TrigDec.make<TH1D>("Single Lepton Trigger Decision", singleleptontrigger_.c_str(), 2, -0.5, 1.5);
      CrossTriggerCombinedHist = subDir_TrigDec.make<TH1D>("Added CrossTrigger Trigger Decision", CombinedTrigger.c_str(), 2, -0.5, 1.5);

      subDir_TrigDiffEff = fileService->mkdir( "Trigger_Distributions" );


      CrossTrigger_Pass_JetPtHist = subDir_TrigDiffEff.make<TH1D>("CrossTrigger_Pass_JetPt", "CrossTriggerPass_Pt", 100, 0, 300);
      CrossTrigger_Total_JetPtHist = subDir_TrigDiffEff.make<TH1D>("CrossTrigger_Total_JetPt", "Total_Pt", 100, 0, 300);
      CrossTrigger_Pass_JetEtaHist = subDir_TrigDiffEff.make<TH1D>("CrossTrigger_Pass_JetEta", "CrossTriggerPass_Eta", 100, -3, 3);
      CrossTrigger_Total_JetEtaHist = subDir_TrigDiffEff.make<TH1D>("CrossTrigger_Total_JetEta", "Total_Eta", 100, -3, 3);

      if ( hadronicleg_ == "SingleTop" ){
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );
            Filter1_Pt = subDir_Filter1.make<TH1D>("Pt", "Pt", 100, 0, 300);
            Filter1_Eta = subDir_Filter1.make<TH1D>("Eta", "Eta", 100, -5, 5);
            Filter1_Phi = subDir_Filter1.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);

            subDir_Filter2 = fileService->mkdir( filter2_.c_str() );
            Filter2_Pt = subDir_Filter2.make<TH1D>("Pt", "Pt", 100, 0, 300);
            Filter2_Eta = subDir_Filter2.make<TH1D>("Eta", "Eta", 100, -5, 5);
            Filter2_Phi = subDir_Filter2.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);
      }
      
      if ( hadronicleg_ == "TTBarJet30" ){
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );
            Filter1_Pt = subDir_Filter1.make<TH1D>("Pt", "Pt", 100, 0, 300);
            Filter1_Eta = subDir_Filter1.make<TH1D>("Eta", "Eta", 100, -5, 5);
            Filter1_Phi = subDir_Filter1.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);
      }
      
      if ( hadronicleg_ == "TTBarJet304050" ){
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );
            Filter1_Pt = subDir_Filter1.make<TH1D>("Pt", "Pt", 100, 0, 300);
            Filter1_Eta = subDir_Filter1.make<TH1D>("Eta", "Eta", 100, -5, 5);
            Filter1_Phi = subDir_Filter1.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);

            subDir_Filter2 = fileService->mkdir( filter2_.c_str() );
            Filter2_Pt = subDir_Filter2.make<TH1D>("Pt", "Pt", 100, 0, 300);
            Filter2_Eta = subDir_Filter2.make<TH1D>("Eta", "Eta", 100, -5, 5);
            Filter2_Phi = subDir_Filter2.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);

            subDir_Filter3 = fileService->mkdir( filter3_.c_str() );
            Filter3_Pt = subDir_Filter3.make<TH1D>("Pt", "Pt", 100, 0, 300);
            Filter3_Eta = subDir_Filter3.make<TH1D>("Eta", "Eta", 100, -5, 5);
            Filter3_Phi = subDir_Filter3.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);
      }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyser::endJob(){

}

// ------------ method called when starting to processes a run  ------------
/*
void 
TriggerAnalyser::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
TriggerAnalyser::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TriggerAnalyser::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TriggerAnalyser::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyser);
