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

            // remove non-jet trigger objects (85=Jet, 86=BJet)
            isJet = false; 
            for (unsigned h = 0; h < obj.filterIds().size(); ++h){
                  if (obj.filterIds()[h] == 85 || obj.filterIds()[h] == 86 ) isJet = true;
            }
            if (!isJet) continue;

            // match trigger jet object to RECO jet and continue if no matching is found
            double minDR2 = 9999;
            int JetIndex = 0;
            isMatched = false;

            for( auto jet = jets->begin(); jet != jets->end(); ++jet ){
   
                  double const dR2 = reco::deltaR2(obj, *jet);

                  // Find index of best matched jet
                  if (dR2 < minDR2){
                        minDR2 = dR2;
                        matchedJetIndex = JetIndex;
                  }
                  ++JetIndex;
            }

            if (minDR2 < 0.3 * 0.3) isMatched = true;
            if (!isMatched) continue;

            // See if trigger jet object passes filter and store filter distributions in histograms
            // Store matched RECO jet pt / btag
           for (unsigned int i = 0, n = obj.filterLabels().size(); i < n; ++i) {

                  if ( hadronicleg_ == "SingleTop" ){
                        if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {

                              Filter1_Pt->Fill(obj.pt());
                              Filter1_Eta->Fill(obj.eta());
                              Filter1_Phi->Fill(obj.phi());
                              Filter1_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());
                        }               

                        if ( obj.filterLabels()[i].find(filter2_) != std::string::npos ) {
                              Filter2_Pt->Fill(obj.pt());
                              Filter2_Eta->Fill(obj.eta());
                              Filter2_Phi->Fill(obj.phi());
                              Filter2_matchedJetBTag->Fill(-log(1-jets->at(matchedJetIndex).bDiscriminator(btagger_)));
                        }

                        Total_matchedJetPt->Fill(jets->at(matchedJetIndex).pt()); 
                        Total_matchedJetBTag->Fill(-log(1-jets->at(matchedJetIndex).bDiscriminator(btagger_)));               
                  }


                  if ( hadronicleg_ == "TTBarJet30" ){
                        if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {
                              Filter1_Pt->Fill(obj.pt());
                              Filter1_Eta->Fill(obj.eta());
                              Filter1_Phi->Fill(obj.phi());
                              Filter1_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());
                        }
                        Total_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());                 
                  }


                  if ( hadronicleg_ == "TTBarJet304050" ){
                        if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {
                              Filter1_Pt->Fill(obj.pt());
                              Filter1_Eta->Fill(obj.eta());
                              Filter1_Phi->Fill(obj.phi());
                              Filter1_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());
                        }
                  
                        if ( obj.filterLabels()[i].find(filter2_) != std::string::npos ) {
                              Filter2_Pt->Fill(obj.pt());
                              Filter2_Eta->Fill(obj.eta());
                              Filter2_Phi->Fill(obj.phi());
                              Filter2_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());
                        }
                     
                        if ( obj.filterLabels()[i].find(filter3_) != std::string::npos ) {
                              Filter3_Pt->Fill(obj.pt());
                              Filter3_Eta->Fill(obj.eta());
                              Filter3_Phi->Fill(obj.phi());
                              Filter3_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());
                        }
                        Total_matchedJetPt->Fill(jets->at(matchedJetIndex).pt());                                         
                  }
            }
      }
}

// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyser::beginJob(){


      // INITIALISE DIRECTORIES AND HISTOGRAMS ------------------------------------------------------ //

      subDir_TrigDec = fileService->mkdir( "Trigger Decision" );
      CrossTriggerHist = subDir_TrigDec.make<TH1D>("CrossTrigger Trigger Decision", crosstrigger_.c_str(), 2, -0.5, 1.5);
      SingleLeptonHist = subDir_TrigDec.make<TH1D>("Single Lepton Trigger Decision", singleleptontrigger_.c_str(), 2, -0.5, 1.5);
      CrossTriggerCombinedHist = subDir_TrigDec.make<TH1D>("Added CrossTrigger Trigger Decision", CombinedTrigger.c_str(), 2, -0.5, 1.5);

      subDir_Observables = fileService->mkdir( "Trigger Observables" );

      subDir_Observables_Jet = subDir_Observables.mkdir( "Jets" );
      CrossTrigger_Pass_JetPtHist = subDir_Observables_Jet.make<TH1D>("CrossTrigger_Pass_JetPt", "Pass_Pt", 100, 0, 200);
      CrossTrigger_Total_JetPtHist = subDir_Observables_Jet.make<TH1D>("CrossTrigger_Total_JetPt", "Total_Pt", 100, 0, 200);
      CrossTrigger_Pass_JetEtaHist = subDir_Observables_Jet.make<TH1D>("CrossTrigger_Pass_JetEta", "Eta", 100, -3, 3);
      CrossTrigger_Total_JetEtaHist = subDir_Observables_Jet.make<TH1D>("CrossTrigger_Total_JetEta", "Total_Eta", 100, -3, 3);
      CrossTrigger_Pass_JetMultiplicity = subDir_Observables_Jet.make<TH1D>("CrossTrigger_Pass_JetMultiplicity", "Pass_JetMultiplicity", 15, 0, 15);
      CrossTrigger_Total_JetMultiplicity = subDir_Observables_Jet.make<TH1D>("CrossTrigger_Total_JetMultiplicity", "Total_JetMultiplicity", 15, 0, 15);
      CrossTrigger_Pass_hltHT = subDir_Observables_Jet.make<TH1D>("CrossTrigger_Pass_hltHT", "Pass_hltHT", 200, 0, 500);
      CrossTrigger_Total_hltHT = subDir_Observables_Jet.make<TH1D>("CrossTrigger_Total_hltHT", "Total_hltHT", 200, 0, 500);
      CrossTrigger_Pass_greatestBtag = subDir_Observables_Jet.make<TH1D>("CrossTrigger_Pass_GreatestBtag", "Pass_GreatestBtag", 100, 0, 20);
      CrossTrigger_Total_greatestBtag = subDir_Observables_Jet.make<TH1D>("CrossTrigger_Total_GreatestBtag", "Total_GreatestBtag", 100, 0, 20);

      subDir_Observables_MET = subDir_Observables.mkdir( "MET" );
      CrossTrigger_Pass_METEnergyHist = subDir_Observables_MET.make<TH1D>("CrossTrigger_Pass_METEnergy", "Pass_Energy", 100, 0, 300);
      CrossTrigger_Total_METEnergyHist = subDir_Observables_MET.make<TH1D>("CrossTrigger_Total_METEnergy", "Total_Energy", 100, 0, 300);

      subDir_Observables_Lepton = subDir_Observables.mkdir( "Leading_Lepton" );
      CrossTrigger_Pass_LeptonPtHist = subDir_Observables_Lepton.make<TH1D>("CrossTrigger_Pass_LeptonPt", "Pass_Pt", 100, 0, 200);
      CrossTrigger_Total_LeptonPtHist = subDir_Observables_Lepton.make<TH1D>("CrossTrigger_Total_LeptonPt", "Total_Pt", 100, 0, 200);
      CrossTrigger_Pass_LeptonEtaHist = subDir_Observables_Lepton.make<TH1D>("CrossTrigger_Pass_LeptonEta", "Pass_Eta", 100, -3, 3);
      CrossTrigger_Total_LeptonEtaHist = subDir_Observables_Lepton.make<TH1D>("CrossTrigger_Total_LeptonEta", "Total_Eta", 100, -3, 3);
      CrossTrigger_Pass_LeptonEnergyHist = subDir_Observables_Lepton.make<TH1D>("CrossTrigger_Pass_LeptonEnergy", "Pass_Energy", 100, 0, 200);
      CrossTrigger_Total_LeptonEnergyHist = subDir_Observables_Lepton.make<TH1D>("CrossTrigger_Total_LeptonEnergy", "Total_Energy", 100, 0, 200);

     if ( hadronicleg_ == "SingleTop" ){

            // Single Lepton Filter
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Filter Observables" );
            Filter1_Pt = subDir_Filter1_Observables.make<TH1D>("Pt", "Pt", 50, 0, 200);
            Filter1_Eta = subDir_Filter1_Observables.make<TH1D>("Eta", "Eta", 50, -5, 5);
            Filter1_Phi = subDir_Filter1_Observables.make<TH1D>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            Filter1_matchedJetPt = subDir_Filter1_MatchedJetObservables.make<TH1D>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt", 50, 0, 200);
            Total_matchedJetPt = subDir_Filter1_MatchedJetObservables.make<TH1D>("matched Jet Pt", "matched Jet Pt", 50, 0, 200);

            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );

            // CSV Filter
            subDir_Filter2 = fileService->mkdir( filter2_.c_str() );

            subDir_Filter2_Observables = subDir_Filter2.mkdir( "Filter Observables" );
            Filter2_Pt = subDir_Filter2_Observables.make<TH1D>("Pt", "Pt", 50, 0, 200);
            Filter2_Eta = subDir_Filter2_Observables.make<TH1D>("Eta", "Eta", 50, -5, 5);
            Filter2_Phi = subDir_Filter2_Observables.make<TH1D>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter2_MatchedJetObservables = subDir_Filter2.mkdir( "matched RECO Jet Observables" );
            Filter2_matchedJetBTag = subDir_Filter2_MatchedJetObservables.make<TH1D>("Filter 2 matched Jet CSV", "Filter 2 matched Jet CSV", 50, 0, 15);
            Total_matchedJetBTag = subDir_Filter2_MatchedJetObservables.make<TH1D>("matched Jet CSV", "matched Jet CSV", 50, 0, 15);

            subDir_Filter2_TurnOnCurves = subDir_Filter2.mkdir( "Turn On Curves" );
      }
      
      if ( hadronicleg_ == "TTBarJet30" ){

            // 30GeV Jet Filter
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Filter Observables" );
            Filter1_Pt = subDir_Filter1_Observables.make<TH1D>("Pt", "Pt", 50, 0, 200);
            Filter1_Eta = subDir_Filter1_Observables.make<TH1D>("Eta", "Eta", 50, -5, 5);
            Filter1_Phi = subDir_Filter1_Observables.make<TH1D>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            Filter1_matchedJetPt = subDir_Filter1_MatchedJetObservables.make<TH1D>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt", 50, 0, 200);
            Total_matchedJetPt = subDir_Filter1_MatchedJetObservables.make<TH1D>("matched Jet Pt", "matched Jet Pt", 50, 0, 200);

            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );
      }
      
      if ( hadronicleg_ == "TTBarJet304050" ){

            // 30GeV Jet Filter
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Filter Observables" );
            Filter1_Pt = subDir_Filter1_Observables.make<TH1D>("Pt", "Pt", 50, 0, 200);
            Filter1_Eta = subDir_Filter1_Observables.make<TH1D>("Eta", "Eta", 50, -5, 5);
            Filter1_Phi = subDir_Filter1_Observables.make<TH1D>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            Filter1_matchedJetPt = subDir_Filter1_MatchedJetObservables.make<TH1D>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt", 50, 0, 200);
            Total_matchedJetPt = subDir_Filter1_MatchedJetObservables.make<TH1D>("matched Jet Pt", "matched Jet Pt", 50, 0, 200);

            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );

            // 40GeV Jet Filter
            subDir_Filter2 = fileService->mkdir( filter2_.c_str() );

            subDir_Filter2_Observables = subDir_Filter2.mkdir( "Filter Observables" );
            Filter2_Pt = subDir_Filter2_Observables.make<TH1D>("Pt", "Pt", 50, 0, 200);
            Filter2_Eta = subDir_Filter2_Observables.make<TH1D>("Eta", "Eta", 50, -5, 5);
            Filter2_Phi = subDir_Filter2_Observables.make<TH1D>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter2_MatchedJetObservables = subDir_Filter2.mkdir( "matched RECO Jet Observables" );
            Filter2_matchedJetPt = subDir_Filter2_MatchedJetObservables.make<TH1D>("Filter 2 matched Jet Pt", "Filter 2 matched Jet Pt", 50, 0, 200);
            Total_matchedJetPt = subDir_Filter2_MatchedJetObservables.make<TH1D>("matched Jet Pt", "matched Jet Pt", 50, 0, 200);

            subDir_Filter2_TurnOnCurves = subDir_Filter2.mkdir( "Turn On Curves" );

            // 50GeV Jet Filter
            subDir_Filter3 = fileService->mkdir( filter3_.c_str() );

            subDir_Filter3_Observables = subDir_Filter3.mkdir( "Filter Observables" );
            Filter3_Pt = subDir_Filter3_Observables.make<TH1D>("Pt", "Pt", 50, 0, 200);
            Filter3_Eta = subDir_Filter3_Observables.make<TH1D>("Eta", "Eta", 50, -5, 5);
            Filter3_Phi = subDir_Filter3_Observables.make<TH1D>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter3_MatchedJetObservables = subDir_Filter3.mkdir( "matched RECO Jet Observables" );
            Filter3_matchedJetPt = subDir_Filter3_MatchedJetObservables.make<TH1D>("Filter 3 matched Jet Pt", "Filter 3 matched Jet Pt", 50, 0, 200);
            Total_matchedJetPt = subDir_Filter3_MatchedJetObservables.make<TH1D>("matched Jet Pt", "matched Jet Pt", 50, 0, 200);

            subDir_Filter3_TurnOnCurves = subDir_Filter3.mkdir( "Turn On Curves" );
      }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyser::endJob(){


      // CREATE TURN ON CURVES ---------------------------------------------------------------------- //

      if ( hadronicleg_ == "SingleTop" ){

            Filter1_TurnOnCurve_Pt = subDir_Filter1_TurnOnCurves.make<TGraphAsymmErrors>(Filter1_matchedJetPt,Total_matchedJetPt);
            Filter1_TurnOnCurve_Pt->SetTitle("TGraphAsymmErrors Example");
            Filter1_TurnOnCurve_Pt->SetMarkerColor(4);
            Filter1_TurnOnCurve_Pt->SetMarkerStyle(21);
            Filter1_TurnOnCurve_Pt->Draw("ALP");
            Filter1_TurnOnCurve_Pt->Write();

            Filter2_TurnOnCurve_BTag = subDir_Filter2_TurnOnCurves.make<TGraphAsymmErrors>(Filter2_matchedJetBTag,Total_matchedJetBTag);
            Filter2_TurnOnCurve_BTag->SetTitle("TGraphAsymmErrors Example");
            Filter2_TurnOnCurve_BTag->SetMarkerColor(4);
            Filter2_TurnOnCurve_BTag->SetMarkerStyle(21);
            Filter2_TurnOnCurve_BTag->Draw("ALP");
            Filter2_TurnOnCurve_BTag->Write();
      }

      if ( hadronicleg_ == "TTBarJet30" ){

            Filter1_TurnOnCurve_Pt = subDir_Filter1_TurnOnCurves.make<TGraphAsymmErrors>(Filter1_matchedJetPt,Total_matchedJetPt);
            Filter1_TurnOnCurve_Pt->SetTitle("TGraphAsymmErrors Example");
            Filter1_TurnOnCurve_Pt->SetMarkerColor(4);
            Filter1_TurnOnCurve_Pt->SetMarkerStyle(21);
            Filter1_TurnOnCurve_Pt->Draw("ALP");
            Filter1_TurnOnCurve_Pt->Write();
      }

      if ( hadronicleg_ == "TTBarJet304050" ){

            Filter1_TurnOnCurve_Pt = subDir_Filter1_TurnOnCurves.make<TGraphAsymmErrors>(Filter1_matchedJetPt,Total_matchedJetPt);
            Filter1_TurnOnCurve_Pt->SetTitle("TGraphAsymmErrors Example");
            Filter1_TurnOnCurve_Pt->SetMarkerColor(4);
            Filter1_TurnOnCurve_Pt->SetMarkerStyle(21);
            Filter1_TurnOnCurve_Pt->Draw("ALP");
            Filter1_TurnOnCurve_Pt->Write();

            Filter2_TurnOnCurve_Pt = subDir_Filter2_TurnOnCurves.make<TGraphAsymmErrors>(Filter2_matchedJetPt,Total_matchedJetPt);
            Filter2_TurnOnCurve_Pt->SetTitle("TGraphAsymmErrors Example");
            Filter2_TurnOnCurve_Pt->SetMarkerColor(4);
            Filter2_TurnOnCurve_Pt->SetMarkerStyle(21);
            Filter2_TurnOnCurve_Pt->Draw("ALP");
            Filter2_TurnOnCurve_Pt->Write();

            Filter3_TurnOnCurve_Pt = subDir_Filter3_TurnOnCurves.make<TGraphAsymmErrors>(Filter3_matchedJetPt,Total_matchedJetPt);
            Filter3_TurnOnCurve_Pt->SetTitle("TGraphAsymmErrors Example");
            Filter3_TurnOnCurve_Pt->SetMarkerColor(4);
            Filter3_TurnOnCurve_Pt->SetMarkerStyle(21);
            Filter3_TurnOnCurve_Pt->Draw("ALP");
            Filter3_TurnOnCurve_Pt->Write();
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
