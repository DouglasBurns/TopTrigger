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
    singleleptontrigger_(iConfig.getParameter <std::string> ("SingleLeptonTriggerInput")),
    crosstrigger_(iConfig.getParameter <std::string> ("CrossTriggerInput")),
    filter1_(iConfig.getParameter <std::string> ("FilterInput1")),
    filter2_(iConfig.getParameter <std::string> ("FilterInput2")),
    filter3_(iConfig.getParameter <std::string> ("FilterInput3")),
    btagger_(iConfig.getParameter <std::string> ("BTagger")),
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
      edm::Handle < std::vector<pat::Electron> > electrons;
      edm::Handle < std::vector<pat::Muon> > muons;

      iEvent.getByToken(triggerResults_, triggerResults);
      iEvent.getByToken(triggerObjects_, triggerObjects);
      iEvent.getByToken(jets_, jets);
      iEvent.getByToken(mets_, mets);
      iEvent.getByToken(electrons_, electrons);
      iEvent.getByToken(muons_, muons);

      const edm::TriggerNames &TrigNames = iEvent.triggerNames(*triggerResults); 


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
            SingleLeptonHist->Fill(SingleLeptonTrigDecision);

            if (SingleLeptonTrigDecision == true){
                  // If passes single lepton trigger find and store trigger results for the cross trigger in histogram
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

      // Find trigger result for cross trigger and store in histogram
      if ( crossIndex < triggerResults->size() ) {
            CrossTriggerTrigDecision = triggerResults->accept(crossIndex);
            CrossTriggerHist->Fill(CrossTriggerTrigDecision);
      }
      else {
            std::cout << "Looking for : " << crosstrigger_ << " but failed" << std::endl;
      }


      // JETS --------------------------------------------------------------------------------------- //

      jetMultiplicity = 0;
      hltHT = 0;
      jetCSV = 0;

      // Select jets and store distributions
      for( auto jet = jets->begin(); jet != jets->end(); ++jet ){ 
            // skip jets with low pT or outside the tracker acceptance
            if( jet->pt()<30. || std::abs(jet->eta())>3.0 ) continue; 
 
            // std::cout << "Btag Value : " << jet->bDiscriminator(btagger_) << std::endl;
            ++jetMultiplicity;
            hltHT += jet->pt();

            // Greatest CSV value
            if (jet->bDiscriminator(btagger_) >= jetCSV) jetCSV = jet->bDiscriminator(btagger_);

            // Store Jet Pt, Eta, Higest CSV, HT and Multiplicity in histograms 
            if (CrossTriggerTrigDecision==true){
                  CrossTrigger_Pass_JetPtHist->Fill(jet->pt());
                  CrossTrigger_Pass_JetEtaHist->Fill(jet->eta());
            }
            CrossTrigger_Total_JetPtHist->Fill(jet->pt());
            CrossTrigger_Total_JetEtaHist->Fill(jet->eta());
      }

      if(CrossTriggerTrigDecision==true){
            CrossTrigger_Pass_JetMultiplicity->Fill(jetMultiplicity);
            CrossTrigger_Pass_hltHT->Fill(hltHT);
            CrossTrigger_Pass_greatestBtag->Fill(-log(1-jetCSV));
      }

      CrossTrigger_Total_JetMultiplicity->Fill(jetMultiplicity);
      CrossTrigger_Total_hltHT->Fill(hltHT);
      CrossTrigger_Total_greatestBtag->Fill(-log(1-jetCSV));


      // MET ---------------------------------------------------------------------------------------- //

      // Store MET of event
      for( auto met = mets->begin(); met != mets->end(); ++met ){ 

            metEnergy = met->energy();//met pt = met energy

            if(CrossTriggerTrigDecision==true){
                  CrossTrigger_Pass_METEnergyHist->Fill(metEnergy);
            }
            CrossTrigger_Total_METEnergyHist->Fill(metEnergy);
      }


      // LEADING LEPTONS ---------------------------------------------------------------------------- //

      // Store the leading lepton pt, eta and energy distributions in histograms
      if ( leptonicleg_ == "Ele" ){

            //Leading lepton at (0)
            if (electrons->size() != 0){ 
                  auto lepton = electrons->at(0);

                  leptonPt = lepton.pt();
                  leptonEta = lepton.eta();
                  leptonEnergy = lepton.energy();

                  if(CrossTriggerTrigDecision==true){
                        CrossTrigger_Pass_LeptonPtHist->Fill(leptonPt);
                        CrossTrigger_Pass_LeptonEtaHist->Fill(leptonEta);
                        CrossTrigger_Pass_LeptonEnergyHist->Fill(leptonEnergy);
                  }
                  CrossTrigger_Total_LeptonPtHist->Fill(leptonPt);
                  CrossTrigger_Total_LeptonEtaHist->Fill(leptonEta);
                  CrossTrigger_Total_LeptonEnergyHist->Fill(leptonEnergy);
            }
      }
      else if ( leptonicleg_ == "Mu" ){
            if (muons->size() != 0){
                  auto lepton = muons->at(0);//leading lepton

                  leptonPt = lepton.pt();
                  leptonEta = lepton.eta();
                  leptonEnergy = lepton.energy();

                  if(CrossTriggerTrigDecision==true){
                        CrossTrigger_Pass_LeptonPtHist->Fill(leptonPt);
                        CrossTrigger_Pass_LeptonEtaHist->Fill(leptonEta);
                        CrossTrigger_Pass_LeptonEnergyHist->Fill(leptonEnergy);
                  }
                  CrossTrigger_Total_LeptonPtHist->Fill(leptonPt);
                  CrossTrigger_Total_LeptonEtaHist->Fill(leptonEta);
                  CrossTrigger_Total_LeptonEnergyHist->Fill(leptonEnergy);
            }
      }


      // FILTERS: TRIGGER OBJECTS ------------------------------------------------------------------- //

      // Get all filter objects in an event
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 

            // remove non-jet trigger objects
            isJet = false; 
            isBJet = false;           
            for (unsigned h = 0; h < obj.filterIds().size(); ++h){
                  if (obj.filterIds()[h] == 85 || obj.filterIds()[h] == 86 ) isJet = true;
                  if (obj.filterIds()[h] == 86 ) isBJet = true;
                  // std::cout << obj.filterIds()[h] << " ";
            }
            // std::cout << std::endl;
            if (!isJet) continue;
            // if (obj.pt() < 30.) continue;

            // match trigger object to RECO jet and continue if no matching is found
            double minDR2 = 9999;
            int JetIndex = 0;
            isMatched = false;
            for( auto jet = jets->begin(); jet != jets->end(); ++jet ){

                  //lose softer jets
                  // if (jet->pt() < 8.) continue;
                      
                  double const dR2 = reco::deltaR2(obj, *jet);
                  // std::cout << dR2 << std::endl;

                  // Find index of matched jet
                  if (dR2 < minDR2){
                        minDR2 = dR2;
                        matchedJetIndex = JetIndex;
                        // std::cout << "New Min : " << minDR2 << std::endl;
                  }
                  ++JetIndex;
            }

            if (minDR2 < 0.3 * 0.3) isMatched = true;
            if (!isMatched) continue;

            // See if trigger object passes filter and store distributions in histograms
            passFilter1=false;
            passFilter2=false;
            passFilter3=false;
           for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {

                  if ( hadronicleg_ == "SingleTop" ){
                        // std::cout << obj.filterLabels()[i] << std::endl;
                        if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {

                              // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h] << std::endl;
                              passFilter1 = true;
                              Filter1_Pt->Fill(obj.pt());
                              Filter1_Eta->Fill(obj.eta());
                              Filter1_Phi->Fill(obj.phi());

                              if (isMatched){
                                    SingleTop_Filter1_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());
                              }
                              // std::cout << "\t   Collection: " << obj.collection() << std::endl;//85 for jets, 86 for Bjets, 0 for No idea - HLTReco/interface/TriggerTypeDefs.h
                              // for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterLabels()[i] << " " << obj.filterIds()[h] << std::endl;
                        }               

                        if ( obj.filterLabels()[i].find(filter2_) != std::string::npos ) {
                              Filter2_Pt->Fill(obj.pt());
                              Filter2_Eta->Fill(obj.eta());
                              Filter2_Phi->Fill(obj.phi());

                              if (isMatched){
                                    Filter2_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());
                                    Filter2_matchedJetBTag->Fill(-log(1-jets->at(matchedJetIndex).bDiscriminator(btagger_)));
                              }
                        }
                        SingleTop_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());                 
                  }

                  if ( hadronicleg_ == "TTBarJet30" ){
                        if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {
                              Filter1_Pt->Fill(obj.pt());
                              Filter1_Eta->Fill(obj.eta());
                              Filter1_Phi->Fill(obj.phi());

                              if (isMatched){
                                    Filter1_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());
                              }
                        }
                  }

                  if ( hadronicleg_ == "TTBarJet304050" ){
                        if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {
                              Filter1_Pt->Fill(obj.pt());
                              Filter1_Eta->Fill(obj.eta());
                              Filter1_Phi->Fill(obj.phi());

                              if (isMatched){
                                    Filter1_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());
                              }
                        }
                  
                        if ( obj.filterLabels()[i].find(filter2_) != std::string::npos ) {
                              Filter2_Pt->Fill(obj.pt());
                              Filter2_Eta->Fill(obj.eta());
                              Filter2_Phi->Fill(obj.phi());

                              if (isMatched){
                                    Filter2_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());
                              }
                        }
                     
                        if ( obj.filterLabels()[i].find(filter3_) != std::string::npos ) {
                              Filter3_Pt->Fill(obj.pt());
                              Filter3_Eta->Fill(obj.eta());
                              Filter3_Phi->Fill(obj.phi());

                              if (isMatched){
                                    Filter3_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());
                              }
                        }
                  }
            }

            //Store efficiency of filter
            // if (isMatched){
                  // if ( hadronicleg_ == "SingleTop" ){
                        // std::cout << "Btag Value : " << jets->at(matchedJetIndex).bDiscriminator(btagger_) << std::endl;
                        // std::cout << "Pass or Fail : " << passFilter1 << " Matched Index : " << matchedJetIndex << " Matched Pt : " << jets->at(matchedJetIndex).pt() << std::endl;
                        // blarg = jets->at(matchedJetIndex).pt();
                        // Filter1_TurnOnCurve_Pt->Fill( passFilter1, blarg );
                  // }      // 
            // }
      }
}










// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyser::beginJob(){

      subDir_TrigDec = fileService->mkdir( "Trigger Decision" );
      CrossTriggerHist = subDir_TrigDec.make<TH1D>("CrossTrigger Trigger Decision", crosstrigger_.c_str(), 2, -0.5, 1.5);
      SingleLeptonHist = subDir_TrigDec.make<TH1D>("Single Lepton Trigger Decision", singleleptontrigger_.c_str(), 2, -0.5, 1.5);
      CrossTriggerCombinedHist = subDir_TrigDec.make<TH1D>("Added CrossTrigger Trigger Decision", CombinedTrigger.c_str(), 2, -0.5, 1.5);

      subDir_TrigDiffEff = fileService->mkdir( "Trigger Observables" );

      subDir_TrigDiffEff_Jet = subDir_TrigDiffEff.mkdir( "Jets" );
      CrossTrigger_Pass_JetPtHist = subDir_TrigDiffEff_Jet.make<TH1D>("CrossTrigger_Pass_JetPt", "CrossTriggerPass_Pt", 100, 0, 200);
      CrossTrigger_Total_JetPtHist = subDir_TrigDiffEff_Jet.make<TH1D>("CrossTrigger_Total_JetPt", "Total_Pt", 100, 0, 200);
      CrossTrigger_Pass_JetEtaHist = subDir_TrigDiffEff_Jet.make<TH1D>("CrossTrigger_Pass_JetEta", "CrossTriggerPass_Eta", 100, -3, 3);
      CrossTrigger_Total_JetEtaHist = subDir_TrigDiffEff_Jet.make<TH1D>("CrossTrigger_Total_JetEta", "Total_Eta", 100, -3, 3);
      CrossTrigger_Pass_JetMultiplicity = subDir_TrigDiffEff_Jet.make<TH1D>("CrossTrigger_Pass_JetMultiplicity", "Pass_JetMultiplicity", 15, 0, 15);
      CrossTrigger_Total_JetMultiplicity = subDir_TrigDiffEff_Jet.make<TH1D>("CrossTrigger_Total_JetMultiplicity", "Total_JetMultiplicity", 15, 0, 15);
      CrossTrigger_Pass_hltHT = subDir_TrigDiffEff_Jet.make<TH1D>("CrossTrigger_Pass_hltHT", "Pass_hltHT", 200, 0, 500);
      CrossTrigger_Total_hltHT = subDir_TrigDiffEff_Jet.make<TH1D>("CrossTrigger_Total_hltHT", "Total_hltHT", 200, 0, 500);
      CrossTrigger_Pass_greatestBtag = subDir_TrigDiffEff_Jet.make<TH1D>("CrossTrigger_Pass_GreatestBtag", "Pass_GreatestBtag", 100, 0, 20);
      CrossTrigger_Total_greatestBtag = subDir_TrigDiffEff_Jet.make<TH1D>("CrossTrigger_Total_GreatestBtag", "Total_GreatestBtag", 100, 0, 20);

      subDir_TrigDiffEff_MET = subDir_TrigDiffEff.mkdir( "MET" );
      CrossTrigger_Pass_METEnergyHist = subDir_TrigDiffEff_MET.make<TH1D>("CrossTrigger_Pass_METEnergy", "CrossTriggerPass_Energy", 100, 0, 300);
      CrossTrigger_Total_METEnergyHist = subDir_TrigDiffEff_MET.make<TH1D>("CrossTrigger_Total_METEnergy", "Total_Energy", 100, 0, 300);

      subDir_TrigDiffEff_Lepton = subDir_TrigDiffEff.mkdir( "Leading_Lepton" );
      CrossTrigger_Pass_LeptonPtHist = subDir_TrigDiffEff_Lepton.make<TH1D>("CrossTrigger_Pass_LeptonPt", "CrossTriggerPass_Pt", 100, 0, 200);
      CrossTrigger_Total_LeptonPtHist = subDir_TrigDiffEff_Lepton.make<TH1D>("CrossTrigger_Total_LeptonPt", "CrossTriggerTotal_Pt", 100, 0, 200);
      CrossTrigger_Pass_LeptonEtaHist = subDir_TrigDiffEff_Lepton.make<TH1D>("CrossTrigger_Pass_LeptonEta", "CrossTriggerPass_Eta", 100, -3, 3);
      CrossTrigger_Total_LeptonEtaHist = subDir_TrigDiffEff_Lepton.make<TH1D>("CrossTrigger_Total_LeptonEta", "CrossTriggerTotal_Eta", 100, -3, 3);
      CrossTrigger_Pass_LeptonEnergyHist = subDir_TrigDiffEff_Lepton.make<TH1D>("CrossTrigger_Pass_LeptonEnergy", "CrossTriggerPass_Energy", 100, 0, 200);
      CrossTrigger_Total_LeptonEnergyHist = subDir_TrigDiffEff_Lepton.make<TH1D>("CrossTrigger_Total_LeptonEnergy", "CrossTriggerTotal_Energy", 100, 0, 200);

     if ( hadronicleg_ == "SingleTop" ){

            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Filter Observables" );
            Filter1_Pt = subDir_Filter1_Observables.make<TH1D>("Pt", "Pt", 100, 0, 200);
            Filter1_Eta = subDir_Filter1_Observables.make<TH1D>("Eta", "Eta", 100, -5, 5);
            Filter1_Phi = subDir_Filter1_Observables.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            SingleTop_Filter1_matchedJetPt = subDir_Filter1_MatchedJetObservables.make<TH1D>("Filtered matched Jet Pt", "Filtered matched Jet Pt", 100, 0, 200);
            SingleTop_matchedJetPt = subDir_Filter1_MatchedJetObservables.make<TH1D>("matched Jet Pt", "matched Jet Pt", 100, 0, 200);

            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );
            // Filter1_TurnOnCurve_Pt = subDir_Filter1_TurnOnCurves.make<TEfficiency>("Efficiency vs Pt", "Efficiency vs Pt", 100, 0, 200);
            // Filter1_TurnOnCurve_Pt = subDir_Filter1_TurnOnCurves.make<TEfficiency>(*SingleTop_Filter1_matchedJetPt,*SingleTop_matchedJetPt);

            subDir_Filter2 = fileService->mkdir( filter2_.c_str() );

            subDir_Filter2_Observables = subDir_Filter2.mkdir( "Filter Observables" );
            Filter2_Pt = subDir_Filter2_Observables.make<TH1D>("Pt", "Pt", 100, 0, 200);
            Filter2_Eta = subDir_Filter2_Observables.make<TH1D>("Eta", "Eta", 100, -5, 5);
            Filter2_Phi = subDir_Filter2_Observables.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);
            subDir_Filter2_MatchedJetObservables = subDir_Filter2.mkdir( "matched RECO Jet Observables" );
            Filter2_matchedJetPt = subDir_Filter2_MatchedJetObservables.make<TH1D>("matched Jet Pt", "matched Jet Pt", 100, 0, 200);
            Filter2_matchedJetBTag = subDir_Filter2_MatchedJetObservables.make<TH1D>("matched Jet CSV", "matched Jet CSV", 100, 0, 15);
      }
      
      if ( hadronicleg_ == "TTBarJet30" ){

            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Filter Observables" );
            Filter1_Pt = subDir_Filter1_Observables.make<TH1D>("Pt", "Pt", 100, 0, 200);
            Filter1_Eta = subDir_Filter1_Observables.make<TH1D>("Eta", "Eta", 100, -5, 5);
            Filter1_Phi = subDir_Filter1_Observables.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);
            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            Filter1_matchedJetPt = subDir_Filter1_MatchedJetObservables.make<TH1D>("matched Jet Pt", "matched Jet Pt", 100, 0, 200);

      }
      
      if ( hadronicleg_ == "TTBarJet304050" ){

            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Filter Observables" );
            Filter1_Pt = subDir_Filter1_Observables.make<TH1D>("Pt", "Pt", 100, 0, 200);
            Filter1_Eta = subDir_Filter1_Observables.make<TH1D>("Eta", "Eta", 100, -5, 5);
            Filter1_Phi = subDir_Filter1_Observables.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);
            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            Filter1_matchedJetPt = subDir_Filter1_MatchedJetObservables.make<TH1D>("matched Jet Pt", "matched Jet Pt", 100, 0, 200);

            subDir_Filter2 = fileService->mkdir( filter2_.c_str() );

            subDir_Filter2_Observables = subDir_Filter2.mkdir( "Filter Observables" );
            Filter2_Pt = subDir_Filter2_Observables.make<TH1D>("Pt", "Pt", 100, 0, 200);
            Filter2_Eta = subDir_Filter2_Observables.make<TH1D>("Eta", "Eta", 100, -5, 5);
            Filter2_Phi = subDir_Filter2_Observables.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);
            subDir_Filter2_MatchedJetObservables = subDir_Filter2.mkdir( "matched RECO Jet Observables" );
            Filter2_matchedJetPt = subDir_Filter2_MatchedJetObservables.make<TH1D>("matched Jet Pt", "matched Jet Pt", 100, 0, 200);

            subDir_Filter3 = fileService->mkdir( filter3_.c_str() );

            subDir_Filter3_Observables = subDir_Filter3.mkdir( "Filter Observables" );
            Filter3_Pt = subDir_Filter3_Observables.make<TH1D>("Pt", "Pt", 100, 0, 200);
            Filter3_Eta = subDir_Filter3_Observables.make<TH1D>("Eta", "Eta", 100, -5, 5);
            Filter3_Phi = subDir_Filter3_Observables.make<TH1D>("Phi", "Phi", 100, -3.5, 3.5);
            subDir_Filter3_MatchedJetObservables = subDir_Filter3.mkdir( "matched RECO Jet Observables" );
            Filter3_matchedJetPt = subDir_Filter3_MatchedJetObservables.make<TH1D>("matched Jet Pt", "matched Jet Pt", 100, 0, 200);
      }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyser::endJob(){

      subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );

      if ( hadronicleg_ == "SingleTop" ){
 
            std::cout << "Checking Consistency" << std::endl;


            if (TEfficiency::CheckConsistency(*SingleTop_Filter1_matchedJetPt,*SingleTop_matchedJetPt)){
                  // Filter1_TurnOnCurve_Pt = new TEfficiency(*SingleTop_Filter1_matchedJetPt,*SingleTop_matchedJetPt);
                  std::cout << "Check Consistency Passed" << std::endl;
                  TCanvas * Canvas1 = new TCanvas("ETTCuts","ETTCuts", 0, 0, 800, 600);

                  Filter1_TurnOnCurve_Pt = subDir_Filter1_TurnOnCurves.make<TEfficiency>(*SingleTop_Filter1_matchedJetPt,*SingleTop_matchedJetPt);
                  Filter1_TurnOnCurve_Pt->Write();
                  Filter1_TurnOnCurve_Pt->Draw("AP");
                  gPad->Update();
                  Canvas1->Update();
                  Canvas1->SaveAs("TurnOn.png");

            }
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
