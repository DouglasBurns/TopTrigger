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
    triggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTriggerResults"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("HLTriggerObjects"))),
    jets_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    mets_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
    electrons_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
    muons_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
    vertices_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
    singleleptontrigger_(iConfig.getParameter <std::string> ("SingleLeptonTriggerInput")),
    crosstrigger_(iConfig.getParameter <std::string> ("CrossTriggerInput")),
    filter1_(iConfig.getParameter <std::string> ("FilterInput1")),
    filter2_(iConfig.getParameter <std::string> ("FilterInput2")),
    filter3_(iConfig.getParameter <std::string> ("FilterInput3")),
    btagger_(iConfig.getParameter <std::string> ("BTagger")),
    electronID_(iConfig.getParameter <std::string> ("ElectronID")),
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

      using namespace edm;

      // std::cout << "In analyze" << std::endl;
      edm::Handle < edm::TriggerResults > triggerResults;
      edm::Handle < pat::TriggerObjectStandAloneCollection > triggerObjects;
      edm::Handle < std::vector<pat::Jet> > jets;
      edm::Handle < std::vector<pat::MET> > mets;
      edm::Handle < std::vector<pat::Electron> > electrons;
      edm::Handle < std::vector<pat::Muon> > muons;
      edm::Handle < std::vector<reco::Vertex> > vertices;

      iEvent.getByToken(triggerResults_, triggerResults);
      iEvent.getByToken(triggerObjects_, triggerObjects);
      iEvent.getByToken(jets_, jets);
      iEvent.getByToken(mets_, mets);
      iEvent.getByToken(electrons_, electrons);
      iEvent.getByToken(muons_, muons);
      iEvent.getByToken(vertices_, vertices);

      const edm::TriggerNames &TrigNames = iEvent.triggerNames(*triggerResults); 


      // MOCK EVENT SELECTION ----------------------------------------------------------------------- //

      // Perform Event Selection
      no_Jets = no_BJets = no_Lepton = 0;
      passMockEventSelection = false;
      for( auto jet = jets->begin(); jet != jets->end(); ++jet ){
            if (jet->pt() > 20){
                  ++no_Jets;
                  if (jet->bDiscriminator(btagger_) > 0.890 ){//find some bjet value (https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X50ns#Supported_Algorithms_and_Operati)
                        ++no_BJets;
                  }
            }
      }

      if (leptonicleg_ == "Ele"){
            for  (auto lepton = electrons->begin(); lepton != electrons->end(); ++lepton){
                  if (lepton->electronID(electronID_)){// && lepton->pt() > 35
                        ++no_Lepton;
                  }
            }
      }
      if (leptonicleg_ == "Mu"){
            for  (auto lepton = muons->begin(); lepton != muons->end(); ++lepton){
                  for( auto vertex = vertices->begin(); vertex != vertices->end(); ++vertex ){ 
                        if (!vertex->isValid()) continue;
                        if (lepton->isTightMuon(*vertex)){//  && lepton->pt() > 25 lepton->
                              ++no_Lepton;
                        }
                  }
            }
      }

      if (hadronicleg_ == "SingleTop"){
            if (no_Jets >= 2 && no_BJets >= 1 && no_Lepton == 1) passMockEventSelection = true;
      }
       if (hadronicleg_ == "TTBarJet30" || hadronicleg_ == "TTBarJet304050"){
            if (no_Jets >= 4 && no_BJets >= 2 && no_Lepton == 1) passMockEventSelection = true;
      }


      // TRIGGER NAMES ------------------------------------------------------------------------------ //

      // Find and store indices of the input triggers
      for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) {
            // std::cout << "Trigger " << i << " in Menu : " << TrigNames.triggerName(i) << std::endl;
            if ( TrigNames.triggerName(i).find(singleleptontrigger_) != std::string::npos ) {
                  singleleptonIndex = i;
            }
            if ( TrigNames.triggerName(i).find(crosstrigger_) != std::string::npos ) {
                  crossIndex = i;
            }
      }


      // TRIGGER RESULTS ---------------------------------------------------------------------------- //

      // Find trigger result for single lepton trigger and store in histogram
      if ( singleleptonIndex < triggerResults->size() ) {
            SingleLeptonTrigDecision = triggerResults->accept(singleleptonIndex);
            histContainer_["SingleLeptonHist"]->Fill(SingleLeptonTrigDecision);
            if (passMockEventSelection){
                  histContainer_["SingleLeptonHist_PassSelection"]->Fill(SingleLeptonTrigDecision);
            }

            if (SingleLeptonTrigDecision == true){
                  // If passes single lepton trigger find and store trigger results for the cross trigger in histogram
                  if ( crossIndex < triggerResults->size() ) {
                        CrossTriggerCombinedTrigDecision = triggerResults->accept(crossIndex);
                        histContainer_["CrossTriggerCombinedHist"]->Fill(CrossTriggerCombinedTrigDecision);

                        if (passMockEventSelection){
                              histContainer_["CrossTriggerCombinedHist_PassSelection"]->Fill(CrossTriggerCombinedTrigDecision);
                        }
                  }

                  else {
                        std::cout << "Looking for : " << crosstrigger_ << " but failed" << std::endl;
                  }
            }
      }
      else {
            std::cout << "Looking for : " << singleleptontrigger_ << " but failed" << std::endl;
      }

      // Find trigger result for cross trigger and store in histogram
      if ( crossIndex < triggerResults->size() ) {
            CrossTriggerTrigDecision = triggerResults->accept(crossIndex);
            histContainer_["CrossTriggerHist"]->Fill(CrossTriggerTrigDecision);
            if (passMockEventSelection){
                  histContainer_["CrossTriggerHist_PassSelection"]->Fill(CrossTriggerTrigDecision);
            }
      }
      else {
            std::cout << "Looking for : " << crosstrigger_ << " but failed" << std::endl;
      }


      if (passMockEventSelection){

            // VERTICES ----------------------------------------------------------------------------------- //

            vertexMultiplicity = 0;
            for( auto vertex = vertices->begin(); vertex != vertices->end(); ++vertex ){ 
                  if (!vertex->isValid()) continue;
                  vertexMultiplicity += 1;
            }

            if (CrossTriggerTrigDecision==true){
                  histContainer_["CrossTrigger_Pass_VertexMultiplicityHist"]->Fill(vertexMultiplicity);
            }

            histContainer_["CrossTrigger_Total_VertexMultiplicityHist"]->Fill(vertexMultiplicity);
            

            // JETS --------------------------------------------------------------------------------------- //

            // Select jets and store distributions
            jetMultiplicity_20 = jetMultiplicity_30 = jetMultiplicity_40 = jetMultiplicity_50 = 0;
            HT = 0;

            jetCSV = 0;
            forwardjeteta = 0;
            for( auto jet = jets->begin(); jet != jets->end(); ++jet ){ 
                  // skip jets with low pT or outside the tracker acceptance
                  if( jet->pt()<20. || std::abs(jet->eta())>2.6 ) continue;

                  // Number of jets in an event and HT depending on jet pt cut
                  if ( jet->pt()>20){
                        ++jetMultiplicity_20;
                  } 
                  if ( jet->pt()>30){
                        ++jetMultiplicity_30;
                  }
                  if ( jet->pt()>40){
                        ++jetMultiplicity_40;
                  } 
                  if ( jet->pt()>50){
                        ++jetMultiplicity_50;
                  } 

                  // HT is the sum of central jets over 10GeV (No MET or Leptons) for events that pass selection
                  HT += jet->pt();

                  // Most forward jet has highest eta
                  if (std::abs(jet->eta()) >= forwardjeteta) forwardjeteta = std::abs(jet->eta());

                  // Jet with greatest btag value
                  if (jet->bDiscriminator(btagger_) >= jetCSV) jetCSV = jet->bDiscriminator(btagger_);

                  // Store Jet Pt, Eta, Higest CSV, Global_HT and Multiplicity in histograms 
                  // if (passMockEventSelection == true){
                  if (CrossTriggerTrigDecision==true){
                        histContainer_["CrossTrigger_Pass_JetPtHist"]->Fill(jet->pt());
                        histContainer_["CrossTrigger_Pass_JetEtaHist"]->Fill(jet->eta());
                        histContainer_["CrossTrigger_Pass_JetPhiHist"]->Fill(jet->phi());
                        histContainer_["CrossTrigger_Pass_JetCSVHist"]->Fill(jet->bDiscriminator(btagger_));
                  }
                  histContainer_["CrossTrigger_Total_JetPtHist"]->Fill(jet->pt());
                  histContainer_["CrossTrigger_Total_JetEtaHist"]->Fill(jet->eta());
                  histContainer_["CrossTrigger_Total_JetPhiHist"]->Fill(jet->phi());
                  histContainer_["CrossTrigger_Total_JetCSVHist"]->Fill(jet->bDiscriminator(btagger_));
            
            }

            if(CrossTriggerTrigDecision==true){
                  histContainer_["CrossTrigger_Pass_JetMultiplicity_20"]->Fill(jetMultiplicity_20);
                  histContainer_["CrossTrigger_Pass_JetMultiplicity_30"]->Fill(jetMultiplicity_30);
                  histContainer_["CrossTrigger_Pass_JetMultiplicity_40"]->Fill(jetMultiplicity_40);
                  histContainer_["CrossTrigger_Pass_JetMultiplicity_50"]->Fill(jetMultiplicity_50);
                  histContainer_["CrossTrigger_Pass_greatestBtag"]->Fill(jetCSV);
                  histContainer_["CrossTrigger_Pass_forwardJetEta"]->Fill(forwardjeteta);
                  histContainer_["CrossTrigger_Pass_HT"]->Fill(HT);
            }


            histContainer_["CrossTrigger_Total_JetMultiplicity_20"]->Fill(jetMultiplicity_20);
            histContainer_["CrossTrigger_Total_JetMultiplicity_30"]->Fill(jetMultiplicity_30);
            histContainer_["CrossTrigger_Total_JetMultiplicity_40"]->Fill(jetMultiplicity_40);
            histContainer_["CrossTrigger_Total_JetMultiplicity_50"]->Fill(jetMultiplicity_50);
            histContainer_["CrossTrigger_Total_greatestBtag"]->Fill(jetCSV);
            histContainer_["CrossTrigger_Total_forwardJetEta"]->Fill(forwardjeteta);
            histContainer_["CrossTrigger_Total_HT"]->Fill(HT);  


            // MET ---------------------------------------------------------------------------------------- //

            // Store MET of event
            for( auto met = mets->begin(); met != mets->end(); ++met ){ 

                  metEnergy = met->energy();//met pt = met energy

                  if(CrossTriggerTrigDecision==true){
                        histContainer_["CrossTrigger_Pass_METHist"]->Fill(metEnergy);
                  }
            histContainer_["CrossTrigger_Total_METHist"]->Fill(metEnergy);
            }


            // LEADING LEPTONS ---------------------------------------------------------------------------- //
            // Store the leading lepton pt, eta and energy distributions in histograms
            leadingLeptonIndex.clear();
            leptonIndex = 0;
            Index = 0;
            if (leptonicleg_ == "Ele"){
                  for  (auto lepton = electrons->begin(); lepton != electrons->end(); ++lepton){
                        if (lepton->electronID(electronID_)){
                              leptonIndex = Index;
                              leadingLeptonIndex.push_back(Index);
                              if(CrossTriggerTrigDecision==true){
                                    histContainer_["CrossTrigger_Pass_LeptonPtHist"]->Fill(lepton->pt());
                                    histContainer_["CrossTrigger_Pass_LeptonEtaHist"]->Fill(lepton->eta());
                                    histContainer_["CrossTrigger_Pass_LeptonPhiHist"]->Fill(lepton->phi());
                                    histContainer_["CrossTrigger_Pass_LeptonEnergyHist"]->Fill(lepton->energy());
                              }
                              histContainer_["CrossTrigger_Total_LeptonPtHist"]->Fill(lepton->pt());
                              histContainer_["CrossTrigger_Total_LeptonEtaHist"]->Fill(lepton->eta());
                              histContainer_["CrossTrigger_Total_LeptonPhiHist"]->Fill(lepton->phi());
                              histContainer_["CrossTrigger_Total_LeptonEnergyHist"]->Fill(lepton->energy());
                        }
                  ++Index;
                  }
            }
            else if ( leptonicleg_ == "Mu" ){
                  for  (auto lepton = muons->begin(); lepton != muons->end(); ++lepton){
                        for( auto vertex = vertices->begin(); vertex != vertices->end(); ++vertex ){ 
                              if (!vertex->isValid()) continue;

                              if (lepton->isTightMuon(*vertex)){
                                    leadingLeptonIndex.push_back(Index);
                                    leptonIndex = Index;

                                    if(CrossTriggerTrigDecision==true){
                                          histContainer_["CrossTrigger_Pass_LeptonPtHist"]->Fill(lepton->pt());
                                          histContainer_["CrossTrigger_Pass_LeptonEtaHist"]->Fill(lepton->eta());
                                          histContainer_["CrossTrigger_Pass_LeptonPhiHist"]->Fill(lepton->phi());
                                          histContainer_["CrossTrigger_Pass_LeptonEnergyHist"]->Fill(lepton->energy());
                                    }
                                    histContainer_["CrossTrigger_Total_LeptonPtHist"]->Fill(lepton->pt());
                                    histContainer_["CrossTrigger_Total_LeptonEtaHist"]->Fill(lepton->eta());
                                    histContainer_["CrossTrigger_Total_LeptonPhiHist"]->Fill(lepton->phi());
                                    histContainer_["CrossTrigger_Total_LeptonEnergyHist"]->Fill(lepton->energy());
                              }
                        }
                  ++Index;
                  }
            }

            // std::cout << "Size : " << leadingLeptonIndex.size() << "  Indices : ";
            // for (uint x = 0; x<leadingLeptonIndex.size(); ++x) std::cout << " " << leadingLeptonIndex[x];
            // std::cout << std::endl;
      

      }
      // FILTERS: TRIGGER OBJECTS ------------------------------------------------------------------- //
      // std::cout << "++++++++++++++++++" << std::endl;

      if (SingleLeptonTrigDecision){

            for ( auto jet = jets->begin(); jet != jets->end(); ++jet ) {

                  // Get rid of non central and soft jets
                  if ( jet->pt()<20. || std::abs(jet->eta())>2.6 ) continue; 
                  
                  // isMatchedToLepton = false;
                  // if ( leadingLeptonIndex.size() != 0 ){
                  //       for (uint x = 0; x<leadingLeptonIndex.size(); ++x){
                  //             if (leptonicleg_ == "Ele"){
                  //                   if (electrons->size() == 0) continue;
                  //                   if ( leadingLeptonIndex.size() <= electrons->size()) continue;
                  //                   auto lepton = electrons->at(leadingLeptonIndex[x]);                              
                  //                   double const dR2 = reco::deltaR2(lepton, *jet);
                  //                   if (dR2 < 0.3 * 0.3) isMatchedToLepton = true;
                  //             }
                  //             if (leptonicleg_ == "Mu"){
                  //                   if (muons->size() == 0) continue;
                  //                   if ( leadingLeptonIndex.size() <= muons->size()) continue;

                  //                   auto lepton = muons->at(leadingLeptonIndex[x]);
                  //                   double const dR2 = reco::deltaR2(lepton, *jet);
                  //                   if (dR2 < 0.3 * 0.3) isMatchedToLepton = true;
                  //             }

                  //             if (isMatchedToLepton) break;
                  //       }  
                  // }
                  // if (isMatchedToLepton) continue;

                  // Fill clean central jets into histograms
                  histContainer_["Total_RECO_JetPt"]->Fill(jet->pt());
                  histContainer_["Total_RECO_JetBTag"]->Fill(-log(1-jet->bDiscriminator(btagger_)));

                  // Find matching filter object,
                  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 

                        // By looking through the filter labels (which filters are associated with filter object)
                        for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {

                              if ( hadronicleg_ == "SingleTop" ){

                                    // Keeping the filter objects that match the filters we are interested in
                                    if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {

                                          // Finally the matching takes place and histograms are filled
                                          double const dR2 = reco::deltaR(obj, *jet);
                                          histContainer_["DeltaR2 Matching"]->Fill(dR2);
                                          if (dR2 < 0.3){
                                                // std::cout << "Woooo! Matched with deltaR : " << dR2 << ", with obj pt of : " << obj.pt() << " to reco pt of : " << jet->pt() << std::endl;
                                                histContainer_["Filter1_Pt"]->Fill(obj.pt());
                                                histContainer_["Filter1_Eta"]->Fill(obj.eta());
                                                histContainer_["Filter1_Phi"]->Fill(obj.phi());

                                                histContainer_["Filter1_matchedJetPt"]->Fill(jet->pt());
                                                histContainer_["Filter1_matchedJetEta"]->Fill(jet->eta());
                                                histContainer_["Filter1_matchedJetPhi"]->Fill(jet->phi());

                                                distContainer_["Filter1_TriggerObject_Reco_Pt"]->Fill(obj.pt(),jet->pt());
                                          }
                                    }
                              

                                    if ( obj.filterLabels()[i].find(filter2_) != std::string::npos ) {
                                          double const dR2 = reco::deltaR2(obj, *jet);
                                          if (dR2 < 0.3 * 0.3){
                                                histContainer_["Filter2_Pt"]->Fill(obj.pt());
                                                histContainer_["Filter2_Eta"]->Fill(obj.eta());
                                                histContainer_["Filter2_Phi"]->Fill(obj.phi());

                                                histContainer_["Filter2_matchedJetPt"]->Fill(jet->pt());
                                                histContainer_["Filter2_matchedJetEta"]->Fill(jet->eta());
                                                histContainer_["Filter2_matchedJetPhi"]->Fill(jet->phi());

                                                histContainer_["Filter2_matchedJetBTag"]->Fill(-log(1-jet->bDiscriminator(btagger_)));
                                          }
                                    }
                              }
                  

                              if ( hadronicleg_ == "TTBarJet30" ){

                                    if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {

                                          double const dR2 = reco::deltaR2(obj, *jet);
                                          if (dR2 < 0.3 * 0.3){
                                                histContainer_["Filter1_Pt"]->Fill(obj.pt());
                                                histContainer_["Filter1_Eta"]->Fill(obj.eta());
                                                histContainer_["Filter1_Phi"]->Fill(obj.phi());

                                                histContainer_["Filter1_matchedJetPt"]->Fill(jet->pt());
                                                histContainer_["Filter1_matchedJetEta"]->Fill(jet->eta());
                                                histContainer_["Filter1_matchedJetPhi"]->Fill(jet->phi());
                                          }
                                    }
                              }


                              if ( hadronicleg_ == "TTBarJet304050" ){
                                    
                                    if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {

                                          double const dR2 = reco::deltaR2(obj, *jet);
                                          if (dR2 < 0.3 * 0.3){
                                                histContainer_["Filter1_Pt"]->Fill(obj.pt());
                                                histContainer_["Filter1_Eta"]->Fill(obj.eta());
                                                histContainer_["Filter1_Phi"]->Fill(obj.phi());

                                                histContainer_["Filter1_matchedJetPt"]->Fill(jet->pt());
                                                histContainer_["Filter1_matchedJetEta"]->Fill(jet->eta());
                                                histContainer_["Filter1_matchedJetPhi"]->Fill(jet->phi());
                                          }
                                    }

                                    if ( obj.filterLabels()[i].find(filter2_) != std::string::npos ) {

                                          double const dR2 = reco::deltaR2(obj, *jet);
                                          if (dR2 < 0.3 * 0.3){
                                                histContainer_["Filter2_Pt"]->Fill(obj.pt());
                                                histContainer_["Filter2_Eta"]->Fill(obj.eta());
                                                histContainer_["Filter2_Phi"]->Fill(obj.phi());

                                                histContainer_["Filter2_matchedJetPt"]->Fill(jet->pt());
                                                histContainer_["Filter2_matchedJetEta"]->Fill(jet->eta());
                                                histContainer_["Filter2_matchedJetPhi"]->Fill(jet->phi());
                                          }
                                    }

                                    if ( obj.filterLabels()[i].find(filter3_) != std::string::npos ) {

                                          double const dR2 = reco::deltaR2(obj, *jet);
                                          if (dR2 < 0.3 * 0.3){
                                                histContainer_["Filter3_Pt"]->Fill(obj.pt());
                                                histContainer_["Filter3_Eta"]->Fill(obj.eta());
                                                histContainer_["Filter3_Phi"]->Fill(obj.phi());

                                                histContainer_["Filter3_matchedJetPt"]->Fill(jet->pt());
                                                histContainer_["Filter3_matchedJetEta"]->Fill(jet->eta());
                                                histContainer_["Filter3_matchedJetPhi"]->Fill(jet->phi());
                                          }
                                    }
                              }
                        }
                  }
            }
      }
}

// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyser::beginJob(){





      // INITIALISE DIRECTORIES AND HISTOGRAMS ------------------------------------------------------ //

      subDir_TrigDec = fileService->mkdir( "Trigger Decision" );
      histContainer_["CrossTriggerHist"] = subDir_TrigDec.make<TH1F>("CrossTrigger Trigger Decision", (crosstrigger_ + "; Fail/Pass; Number of Events").c_str(), 2, -0.5, 1.5);
      histContainer_["SingleLeptonHist"] = subDir_TrigDec.make<TH1F>("Single Lepton Trigger Decision", (singleleptontrigger_ + "; Fail/Pass; Number of Events").c_str(), 2, -0.5, 1.5);
      histContainer_["CrossTriggerCombinedHist"] = subDir_TrigDec.make<TH1F>("Added CrossTrigger Trigger Decision", (crosstrigger_ + " and " + singleleptontrigger_ + "; Fail/Pass; Number of Events").c_str(), 2, -0.5, 1.5);

      histContainer_["CrossTriggerHist_PassSelection"] = subDir_TrigDec.make<TH1F>("CrossTrigger_PassSelection Trigger Decision", (crosstrigger_ + "_PassSelection; Fail/Pass; Number of Events").c_str(), 2, -0.5, 1.5);
      histContainer_["SingleLeptonHist_PassSelection"] = subDir_TrigDec.make<TH1F>("Single Lepton_PassSelection Trigger Decision", (singleleptontrigger_ + "_PassSelection; Fail/Pass; Number of Events").c_str(), 2, -0.5, 1.5);
      histContainer_["CrossTriggerCombinedHist_PassSelection"] = subDir_TrigDec.make<TH1F>("Added CrossTrigger_PassSelection Trigger Decision", (crosstrigger_ + " and " + singleleptontrigger_ + "_PassSelection; Fail/Pass; Number of Events").c_str(), 2, -0.5, 1.5);

      subDir_TrigDec_TurnOnCurves = subDir_TrigDec.mkdir( "Trigger Turn On Curves" );

      subDir_Observables = fileService->mkdir( "Trigger Observables" );
      subDir_Observables_Jet = subDir_Observables.mkdir( "Jets" );
      histContainer_["CrossTrigger_Pass_JetPtHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetPt", "Pass_Pt; RECO Jet Pt (GeV); Number of Events", 30, 0, 120);
      histContainer_["CrossTrigger_Total_JetPtHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetPt", "Total_Pt; RECO Jet Pt (GeV); Number of Events", 30, 0, 120);
      histContainer_["CrossTrigger_Pass_JetEtaHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetEta", "Pass_Eta; RECO Jet Eta; Number of Events", 100, -3, 3);
      histContainer_["CrossTrigger_Total_JetEtaHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetEta", "Total_Eta; RECO Jet Eta; Number of Events", 100, -3, 3);
      histContainer_["CrossTrigger_Pass_JetPhiHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetPhi", "Pass_Eta; Jet Eta; Number of Events", 100, -3.5, 3.5);
      histContainer_["CrossTrigger_Total_JetPhiHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetPhi", "Total_Eta; Jet Eta; Number of Events", 100, -3.5, 3.5);
      histContainer_["CrossTrigger_Pass_JetCSVHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetCSV", "Pass_JetCSV; RECO Jet pfCSVv2; Number of Events", 100, 0, 1);
      histContainer_["CrossTrigger_Total_JetCSVHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetCSV", "Total_JetCSV; RECO Jet pfCSVv2; Number of Events", 100, 0, 1);       
      histContainer_["CrossTrigger_Pass_JetMultiplicity_20"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetMultiplicity_20", "Pass_JetMultiplicity_20; RECO Jet Multiplicity; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Total_JetMultiplicity_20"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetMultiplicity_20", "Total_JetMultiplicity_20; RECO Jet Multiplicity; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Pass_JetMultiplicity_30"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetMultiplicity_30", "Pass_JetMultiplicity_30; RECO Jet Multiplicity; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Total_JetMultiplicity_30"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetMultiplicity_30", "Total_JetMultiplicity_30; RECO Jet Multiplicity; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Pass_JetMultiplicity_40"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetMultiplicity_40", "Pass_JetMultiplicity_40; RECO Jet Multiplicity; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Total_JetMultiplicity_40"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetMultiplicity_40", "Total_JetMultiplicity_40; RECO Jet Multiplicity; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Pass_JetMultiplicity_50"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetMultiplicity_50", "Pass_JetMultiplicity_50; RECO Jet Multiplicity; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Total_JetMultiplicity_50"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetMultiplicity_50", "Total_JetMultiplicity_50; RECO Jet Multiplicity; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Pass_HT"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_HT", "Pass_HT; RECO Jet HT (GeV); Number of Events", 50, 0, 700);
      histContainer_["CrossTrigger_Total_HT"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_HT", "Total_HT; RECO Jet HT (GeV); Number of Events", 50, 0, 700);
      histContainer_["CrossTrigger_Pass_greatestBtag"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_GreatestBtag", "Pass_GreatestBtag; Greatest RECO Jet pfCSVv2; Number of Events", 100, 0, 1);
      histContainer_["CrossTrigger_Total_greatestBtag"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_GreatestBtag", "Total_GreatestBtag; Greatest RECO Jet pfCSVv2; Number of Events;", 100, 0, 12);     
      histContainer_["CrossTrigger_Pass_forwardJetEta"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_forwardJetEta", "Pass_ForwardJetEta; Forward RECO Jet Eta; Number of Events", 100, 0, 3);
      histContainer_["CrossTrigger_Total_forwardJetEta"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_forwardJetEta", "Total_ForwardJetEta; Forward RECO Jet Eta; Number of Events", 100, 0, 3);
      histContainer_["Total_RECO_JetPt"] = subDir_Observables_Jet.make<TH1F>("Total_RECO_Jet_Pt", "RECO Jet Pt; RECO Jet Pt (GeV); Number of Events", 50, 0, 80);
      histContainer_["Total_RECO_JetBTag"] = subDir_Observables_Jet.make<TH1F>("Total_RECO_Jet_BTag", "RECO Jet BTag; RECO Jet BTag pfCSVv2; Number of Events", 50, 0, 10);

      subDir_Observables_Vertices = subDir_Observables.mkdir( "Vertices" );
      histContainer_["CrossTrigger_Pass_VertexMultiplicityHist"] = subDir_Observables_Vertices.make<TH1F>("CrossTrigger_Pass_VertexMultiplicityHist", "Pass_Vertex; Vertex Multiplicity; Number of Events", 30, 0, 30);
      histContainer_["CrossTrigger_Total_VertexMultiplicityHist"] = subDir_Observables_Vertices.make<TH1F>("CrossTrigger_Total_VertexMultiplicityHist", "Total_Vertex; Vertex Multiplicity; Number of Events", 30, 0, 30);

      subDir_Observables_MET = subDir_Observables.mkdir( "MET" );
      histContainer_["CrossTrigger_Pass_METHist"] = subDir_Observables_MET.make<TH1F>("CrossTrigger_Pass_MET", "Pass_MET; MET (GeV); Number of Events", 100, 0, 200);
      histContainer_["CrossTrigger_Total_METHist"] = subDir_Observables_MET.make<TH1F>("CrossTrigger_Total_MET", "Total_MET; MET (GeV); Number of Events", 100, 0, 200);

      subDir_Observables_Lepton = subDir_Observables.mkdir( "Leading_Lepton" );
      histContainer_["CrossTrigger_Pass_LeptonPtHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Pass_LeptonPt", ("Pass_Pt; RECO_" + leptonicleg_ + " Pt (GeV); Number of Events").c_str(), 30, 0, 120);
      histContainer_["CrossTrigger_Total_LeptonPtHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Total_LeptonPt", ("Total_Pt; RECO_" + leptonicleg_ + " Pt (GeV); Number of Events").c_str(), 30, 0, 120);
      histContainer_["CrossTrigger_Pass_LeptonEtaHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Pass_LeptonEta", ("Pass_Eta; RECO_" + leptonicleg_ + " Eta; Number of Events").c_str(), 100, -3, 3);
      histContainer_["CrossTrigger_Total_LeptonEtaHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Total_LeptonEta", ("Total_Eta; RECO_" + leptonicleg_ + " Eta; Number of Events").c_str(), 100, -3, 3);
      histContainer_["CrossTrigger_Pass_LeptonPhiHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Pass_LeptonPhi", ("Pass_Eta; RECO_" + leptonicleg_ + " Eta; Number of Events").c_str(), 100, -3.5, 3.5);
      histContainer_["CrossTrigger_Total_LeptonPhiHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Total_LeptonPhi", ("Total_Eta; RECO_" + leptonicleg_ + " Eta; Number of Events").c_str(), 100, -3.5, 3.5);
      histContainer_["CrossTrigger_Pass_LeptonEnergyHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Pass_LeptonEnergy", ("Pass_Energy; RECO_" + leptonicleg_ + " Energy (GeV); Number of Events").c_str(), 100, 0, 200);
      histContainer_["CrossTrigger_Total_LeptonEnergyHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Total_LeptonEnergy", ("Total_Energy; RECO_" + leptonicleg_ + " Energy (GeV); Number of Events").c_str(), 100, 0, 200);

      
      if ( hadronicleg_ == "SingleTop" ){

            // Single Lepton Filter
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Trigger Object Observables" );
            histContainer_["Filter1_Pt"] = subDir_Filter1_Observables.make<TH1F>("Pt", "Pt; Filter Object Pt; Number of Events", 50, 0, 100);
            histContainer_["Filter1_Eta"] = subDir_Filter1_Observables.make<TH1F>("Eta", "Eta; Filter Object Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter1_Phi"] = subDir_Filter1_Observables.make<TH1F>("Phi", "Phi; Filter Object Phi; Number of Events", 50, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter1_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 50, 0, 80);
            histContainer_["Filter1_matchedJetEta"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Eta", "Filter 1 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter1_matchedJetPhi"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Phi", "Filter 1 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);
            
            histContainer_["DeltaR2 Matching"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("DeltaR2","DeltaR2 Distribution; DeltaR2; Number of Events", 150, 0, 15);
            distContainer_["Filter1_TriggerObject_Reco_Pt"] = subDir_Filter1_MatchedJetObservables.make<TH2F>("Matching", "Filter Object, Matched RECO Jet Pt; Filter Obj Jet Pt (GeV); RECO Jet Pt (GeV)", 50, 0, 200, 50, 0, 200);//cndjschdjkshcjkdnsk

            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );

            // CSV Filter
            subDir_Filter2 = fileService->mkdir( filter2_.c_str() );

            subDir_Filter2_Observables = subDir_Filter2.mkdir( "Trigger Object Observables" );
            histContainer_["Filter2_Pt"] = subDir_Filter2_Observables.make<TH1F>("Pt", "Pt", 50, 0, 100);
            histContainer_["Filter2_Eta"] = subDir_Filter2_Observables.make<TH1F>("Eta", "Eta", 50, -3, 3);
            histContainer_["Filter2_Phi"] = subDir_Filter2_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter2_MatchedJetObservables = subDir_Filter2.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter2_matchedJetBTag"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet CSV", "Filter 2 matched Jet CSV; RECO Jet CSV (pfCSVv2); Number of Events", 50, 0, 10);
            histContainer_["Filter2_matchedJetPt"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Pt", "Filter 2 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 50, 0, 100);
            histContainer_["Filter2_matchedJetEta"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Eta", "Filter 2 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter2_matchedJetPhi"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Phi", "Filter 2 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);
            subDir_Filter2_TurnOnCurves = subDir_Filter2.mkdir( "Turn On Curves" );
      }
      
      if ( hadronicleg_ == "TTBarJet30" ){

            // // 30GeV Jet Filter
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Trigger Object Observables" );
            histContainer_["Filter1_Pt"] = subDir_Filter1_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter1_Eta"] = subDir_Filter1_Observables.make<TH1F>("Eta", "Eta", 50, -3, 3);
            histContainer_["Filter1_Phi"] = subDir_Filter1_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter1_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 50, 0, 80);
            histContainer_["Filter1_matchedJetEta"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Eta", "Filter 1 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter1_matchedJetPhi"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Phi", "Filter 1 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);

            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );
      }
      
      if ( hadronicleg_ == "TTBarJet304050" ){

            // 30GeV Jet Filter
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Filter Observables" );
            histContainer_["Filter1_Pt"] = subDir_Filter1_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter1_Eta"] = subDir_Filter1_Observables.make<TH1F>("Eta", "Eta", 50, -3, 3);
            histContainer_["Filter1_Phi"] = subDir_Filter1_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter1_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 50, 0, 80);
            histContainer_["Filter1_matchedJetEta"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Eta", "Filter 1 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter1_matchedJetPhi"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Phi", "Filter 1 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);

            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );

            // 40GeV Jet Filter
            subDir_Filter2 = fileService->mkdir( filter2_.c_str() );

            subDir_Filter2_Observables = subDir_Filter2.mkdir( "Filter Observables" );
            histContainer_["Filter2_Pt"] = subDir_Filter2_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter2_Eta"] = subDir_Filter2_Observables.make<TH1F>("Eta", "Eta", 50, -3, 3);
            histContainer_["Filter2_Phi"] = subDir_Filter2_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter2_MatchedJetObservables = subDir_Filter2.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter2_matchedJetPt"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Pt", "Filter 2 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 50, 0, 80);
            histContainer_["Filter2_matchedJetEta"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Eta", "Filter 2 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter2_matchedJetPhi"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Phi", "Filter 2 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);

            subDir_Filter2_TurnOnCurves = subDir_Filter2.mkdir( "Turn On Curves" );

            // 50GeV Jet Filter
            subDir_Filter3 = fileService->mkdir( filter3_.c_str() );

            subDir_Filter3_Observables = subDir_Filter3.mkdir( "Filter Observables" );
            histContainer_["Filter3_Pt"] = subDir_Filter3_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter3_Eta"] = subDir_Filter3_Observables.make<TH1F>("Eta", "Eta", 50, -3, 3);
            histContainer_["Filter3_Phi"] = subDir_Filter3_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter3_MatchedJetObservables = subDir_Filter3.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter3_matchedJetPt"] = subDir_Filter3_MatchedJetObservables.make<TH1F>("Filter 3 matched Jet Pt", "Filter 3 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 50, 0, 80);
            histContainer_["Filter3_matchedJetEta"] = subDir_Filter3_MatchedJetObservables.make<TH1F>("Filter 3 matched Jet Eta", "Filter 3 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter3_matchedJetPhi"] = subDir_Filter3_MatchedJetObservables.make<TH1F>("Filter 3 matched Jet Phi", "Filter 3 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);

            subDir_Filter3_TurnOnCurves = subDir_Filter3.mkdir( "Turn On Curves" );
      }
}

// // ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyser::endJob(){
//       //reco vs trigger object




//       // CREATE TRIGGER TURN ON CURVES -------------------------------------------------------------- //
      TCanvas *c1 = new TCanvas((crosstrigger_ + "_JetMult_20").c_str(),"c1",600,400);
      c1->SetGrid();
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_20"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMultiplicity_20"],histContainer_["CrossTrigger_Total_JetMultiplicity_20"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_20"]->SetTitle("Efficiency Plot Jet Multiplicity > 20 GeV; Jet Multiplicity; Efficiency");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_20"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_20"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_20"]->SetMaximum(1);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_20"]->Draw("ALP");
      gPad->Update();
      c1->Update();
      c1->Write();

      TCanvas *c2 = new TCanvas((crosstrigger_ + "_JetMult_30").c_str(),"c2",600,400);
      c2->SetGrid();
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_30"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMultiplicity_30"],histContainer_["CrossTrigger_Total_JetMultiplicity_30"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_30"]->SetTitle("Efficiency Plot Jet Multiplicity > 30 GeV; Jet Multiplicity; Efficiency");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_30"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_30"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_30"]->SetMaximum(1);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_30"]->Draw("ALP");
      gPad->Update();
      c2->Update();
      c2->Write();

      TCanvas *c3 = new TCanvas((crosstrigger_ + "_JetMult_40").c_str(),"c2",600,400);
      c3->SetGrid();
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_40"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMultiplicity_40"],histContainer_["CrossTrigger_Total_JetMultiplicity_40"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_40"]->SetTitle("Efficiency Plot Jet Multiplicity > 40 GeV; Jet Multiplicity; Efficiency");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_40"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_40"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_40"]->SetMaximum(1);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_40"]->Draw("ALP");
      gPad->Update();
      c3->Update();
      c3->Write();

      TCanvas *c4 = new TCanvas((crosstrigger_ + "_JetMult_50").c_str(),"c4",600,400);
      c4->SetGrid();
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_50"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMultiplicity_50"],histContainer_["CrossTrigger_Total_JetMultiplicity_50"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_50"]->SetTitle("Efficiency Plot Jet Multiplicity > 50 GeV; Jet Multiplicity; Efficiency");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_50"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_50"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_50"]->SetMaximum(1);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_50"]->Draw("ALP");
      gPad->Update();
      c4->Update();
      c4->Write();

      TCanvas *c5 = new TCanvas((crosstrigger_ + "_Global_HT_20").c_str(),"c5",600,400);
      c5->SetGrid();
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Global_HT_20"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_HT"],histContainer_["CrossTrigger_Total_HT"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Global_HT_20"]->SetTitle("Efficiency Plot HT (Sum Jets > 20 GeV); Global_HT_20 (GeV); Efficiency");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Global_HT_20"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Global_HT_20"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Global_HT_20"]->SetMaximum(1);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Global_HT_20"]->Draw("ALP");
      gPad->Update();
      c5->Update();
      c5->Write();

      TCanvas *c6 = new TCanvas((crosstrigger_ + "_MET").c_str(),"c6",600,400);
      c6->SetGrid();
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_METHist"],histContainer_["CrossTrigger_Total_METHist"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->SetTitle("Efficiency Plot MET; MET (GeV); Efficiency");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->SetMaximum(1);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->Draw("ALP");
      gPad->Update();
      c6->Update();
      c6->Write();

      TCanvas *c7 = new TCanvas((crosstrigger_ + "_VertMult").c_str(),"c7",600,400);
      c7->SetGrid();
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_VertexMultiplicityHist"],histContainer_["CrossTrigger_Total_VertexMultiplicityHist"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->SetTitle("Efficiency Plot Vertex Multiplicity; Vertex Multiplicity; Efficiency");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->SetMaximum(1);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->Draw("ALP");
      gPad->Update();
      c7->Update();
      c7->Write();

      TCanvas *c8 = new TCanvas((crosstrigger_ + "_LeadingLepton_Pt").c_str(),"c8",600,400);
      c8->SetGrid();
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_LeptonPt"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_LeptonPtHist"],histContainer_["CrossTrigger_Total_LeptonPtHist"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_LeptonPt"]->SetTitle(("Efficiency Plot " + leptonicleg_ + " Pt; " + leptonicleg_ + " Pt (GeV); Efficiency").c_str());
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_LeptonPt"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_LeptonPt"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_LeptonPt"]->SetMaximum(1);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_LeptonPt"]->Draw("ALP");
      gPad->Update();
      c8->Update();
      c8->Write();


      // CREATE FILTER TURN ON CURVES --------------------------------------------------------------- //

      if ( hadronicleg_ == "SingleTop" ){

            TCanvas *c9 = new TCanvas((crosstrigger_ + "_FilterTurnOn_Pt30").c_str(),"c9",600,400);
            c9->SetGrid();
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"] = subDir_Filter1_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter1_matchedJetPt"],histContainer_["Total_RECO_JetPt"]);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetTitle( (filter1_ + " Turn On Jet Pt; RECO Jet Pt (GeV); Efficiency").c_str() );
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerColor(4);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerStyle(21);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMaximum(1);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Draw("ALP");
            gPad->Update();
            c9->Update();
            c9->Write();

            TCanvas *c10 = new TCanvas((crosstrigger_ + "_FilterTurnOn_CSV").c_str(),"c10",600,400);
            c10->SetGrid();
            turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"] = subDir_Filter2_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter2_matchedJetBTag"],histContainer_["Total_RECO_JetBTag"]);
            turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetTitle( (filter2_ + " Turn On Jet CSV; RECO Jet CSV; Efficiency").c_str() );
            turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetMarkerColor(4);
            turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetMarkerStyle(21);
            turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetMaximum(1);
            turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->Draw("ALP");
            gPad->Update();
            c10->Update();
            c10->Write();
      }

      if ( hadronicleg_ == "TTBarJet30" ){

            TCanvas *c9 = new TCanvas((crosstrigger_ + "_FilterTurnOn_Pt30").c_str(),"c9",600,400);
            c9->SetGrid();
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"] = subDir_Filter1_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter1_matchedJetPt"],histContainer_["Total_RECO_JetPt"]);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetTitle( (filter1_ + " Turn On Jet Pt; RECO Jet Pt (GeV); Efficiency").c_str() );
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerColor(4);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerStyle(21);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMaximum(1);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Draw("ALP");
            gPad->Update();
            c9->Update();
            c9->Write();
      }

      if ( hadronicleg_ == "TTBarJet304050" ){

            TCanvas *c9 = new TCanvas((crosstrigger_ + "_FilterTurnOn_Pt30").c_str(),"c9",600,400);
            c9->SetGrid();
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"] = subDir_Filter1_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter1_matchedJetPt"],histContainer_["Total_RECO_JetPt"]);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetTitle( (filter1_ + " Turn On Jet Pt; RECO Jet Pt (GeV); Efficiency").c_str() );
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerColor(4);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerStyle(21);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMaximum(1);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Draw("ALP");
            gPad->Update();
            c9->Update();
            c9->Write();

            TCanvas *c10 = new TCanvas((crosstrigger_ + "_FilterTurnOn_Pt40").c_str(),"c10",600,400);
            c10->SetGrid();
            turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"] = subDir_Filter2_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter2_matchedJetPt"],histContainer_["Total_RECO_JetPt"]);
            turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->SetTitle( (filter2_ + " Turn On Jet Pt; RECO Jet Pt (GeV); Efficiency").c_str() );
            turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->SetMarkerColor(4);
            turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->SetMarkerStyle(21);
            turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->SetMaximum(1);
            turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->Draw("ALP");
            gPad->Update();
            c10->Update();
            c10->Write();

            TCanvas *c11 = new TCanvas((crosstrigger_ + "_FilterTurnOn_Pt50").c_str(),"c11",600,400);
            c11->SetGrid();
            turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"] = subDir_Filter3_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter3_matchedJetPt"],histContainer_["Total_RECO_JetPt"]);
            turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->SetTitle( (filter3_ + " Turn On Jet Pt; RECO Jet Pt (GeV); Efficiency").c_str() );
            turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->SetMarkerColor(4);
            turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->SetMarkerStyle(21);
            turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->SetMaximum(1);
            turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->Draw("ALP");
            gPad->Update();
            c11->Update();
            c11->Write();
      }
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
