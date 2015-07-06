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

            if (SingleLeptonTrigDecision == true){
                  // If passes single lepton trigger find and store trigger results for the cross trigger in histogram
                  if ( crossIndex < triggerResults->size() ) {
                        CrossTriggerCombinedTrigDecision = triggerResults->accept(crossIndex);
                        histContainer_["CrossTriggerCombinedHist"]->Fill(CrossTriggerCombinedTrigDecision);
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
      }
      else {
            std::cout << "Looking for : " << crosstrigger_ << " but failed" << std::endl;
      }


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

      jetMultiplicity = 0;
      hltHT = 0;
      jetCSV = 0;
      forwardjeteta = 0;

      // Select jets and store distributions
      for( auto jet = jets->begin(); jet != jets->end(); ++jet ){ 
            // skip jets with low pT or outside the tracker acceptance
            if( jet->pt()<30. || std::abs(jet->eta())>3.0 ) continue; 
 
            // std::cout << "Btag Value : " << jet->bDiscriminator(btagger_) << std::endl;
            ++jetMultiplicity;
            hltHT += jet->pt();

            // Most forward jet
            if (std::abs(jet->eta()) >= forwardjeteta) forwardjeteta = std::abs(jet->eta());

            // Greatest CSV value
            if (jet->bDiscriminator(btagger_) >= jetCSV) jetCSV = jet->bDiscriminator(btagger_);

            // Store Jet Pt, Eta, Higest CSV, HT and Multiplicity in histograms 
            if (CrossTriggerTrigDecision==true){
                  histContainer_["CrossTrigger_Pass_JetPtHist"]->Fill(jet->pt());
                  histContainer_["CrossTrigger_Pass_JetEtaHist"]->Fill(jet->eta());
                  histContainer_["CrossTrigger_Pass_JetCSVHist"]->Fill(jet->bDiscriminator(btagger_));
            }
            histContainer_["CrossTrigger_Total_JetPtHist"]->Fill(jet->pt());
            histContainer_["CrossTrigger_Total_JetEtaHist"]->Fill(jet->eta());
            histContainer_["CrossTrigger_Total_JetCSVHist"]->Fill(jet->bDiscriminator(btagger_));

      }

      if(CrossTriggerTrigDecision==true){
            histContainer_["CrossTrigger_Pass_JetMultiplicity"]->Fill(jetMultiplicity);
            histContainer_["CrossTrigger_Pass_hltHT"]->Fill(hltHT);
            histContainer_["CrossTrigger_Pass_greatestBtag"]->Fill(-log(1-jetCSV));
            histContainer_["CrossTrigger_Pass_forwardJetEta"]->Fill(forwardjeteta);
      }

      histContainer_["CrossTrigger_Total_JetMultiplicity"]->Fill(jetMultiplicity);
      histContainer_["CrossTrigger_Total_hltHT"]->Fill(hltHT);
      histContainer_["CrossTrigger_Total_greatestBtag"]->Fill(-log(1-jetCSV));
      histContainer_["CrossTrigger_Total_forwardJetEta"]->Fill(forwardjeteta);


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
      if ( leptonicleg_ == "Ele" ){

            //Leading lepton at (0)
            if (electrons->size() != 0){ 
                  auto lepton = electrons->at(0);

                  leptonPt = lepton.pt();
                  leptonEta = lepton.eta();
                  leptonEnergy = lepton.energy();

                  if(CrossTriggerTrigDecision==true){
                        histContainer_["CrossTrigger_Pass_LeptonPtHist"]->Fill(leptonPt);
                        histContainer_["CrossTrigger_Pass_LeptonEtaHist"]->Fill(leptonEta);
                        histContainer_["CrossTrigger_Pass_LeptonEnergyHist"]->Fill(leptonEnergy);
                  }
                  histContainer_["CrossTrigger_Total_LeptonPtHist"]->Fill(leptonPt);
                  histContainer_["CrossTrigger_Total_LeptonEtaHist"]->Fill(leptonEta);
                  histContainer_["CrossTrigger_Total_LeptonEnergyHist"]->Fill(leptonEnergy);
            }
      }
      else if ( leptonicleg_ == "Mu" ){
            if (muons->size() != 0){
                  auto lepton = muons->at(0);//leading lepton

                  leptonPt = lepton.pt();
                  leptonEta = lepton.eta();
                  leptonEnergy = lepton.energy();

                  if(CrossTriggerTrigDecision==true){
                        histContainer_["CrossTrigger_Pass_LeptonPtHist"]->Fill(leptonPt);
                        histContainer_["CrossTrigger_Pass_LeptonEtaHist"]->Fill(leptonEta);
                        histContainer_["CrossTrigger_Pass_LeptonEnergyHist"]->Fill(leptonEnergy);
                  }
                  histContainer_["CrossTrigger_Total_LeptonPtHist"]->Fill(leptonPt);
                  histContainer_["CrossTrigger_Total_LeptonEtaHist"]->Fill(leptonEta);
                  histContainer_["CrossTrigger_Total_LeptonEnergyHist"]->Fill(leptonEnergy);
            }
      }


      // FILTERS: TRIGGER OBJECTS ------------------------------------------------------------------- //
      // std::cout << "++++++++++++++++++" << std::endl;

      // Get all filter objects in an event - MOVE JETS FIRST?
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 

            // pat::TriggerObject OBJ = obj.triggerObject();//same as above

            // hltAK4PFJetsCorrected::HLT
            // remove non-jet trigger objects (85=Jet, 86=BJet)
            isJet = false; 
            for (unsigned h = 0; h < obj.filterIds().size(); ++h){
                  if (obj.filterIds()[h] == 85 || obj.filterIds()[h] == 86 ) isJet = true;
            }
            if (!isJet) continue;
            if ( (obj.collection() != "hltAK4PFJetsCorrected::HLT") && (obj.collection() != "hltPFJetForBtag::HLT" ) ) continue;



            // std::cout << "Jet Pt : " << obj.pt() << std::endl;
            // std::cout << "Collection = " << obj.collection() << std::endl;
            // for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
            //       std::cout << "FliterLabel : " << obj.filterLabels()[h] << std::endl;
            // }            
            // std::cout << std::endl;



            // match trigger jet object to RECO jet and continue if no matching is found
            double minDR2 = 9999;
            int JetIndex = 0;
            isMatched = false;

            for( auto jet = jets->begin(); jet != jets->end(); ++jet ){
   
                  double const dR2 = reco::deltaR2(obj, *jet);

                  // histContainer_["RecoJet_Pt"]->Fill(jet.pt());


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
                              distContainer_["TriggerObject_Reco_Pt"]->Fill(obj.pt(),jets->at(matchedJetIndex).pt());

                        if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {

                              histContainer_["Filter1_Pt"]->Fill(obj.pt());
                              histContainer_["Filter1_Eta"]->Fill(obj.eta());
                              histContainer_["Filter1_Phi"]->Fill(obj.phi());
                              histContainer_["Filter1_matchedJetPt"]->Fill(jets->at(matchedJetIndex).pt());
                              distContainer_["Filter1_TriggerObject_Reco_Pt"]->Fill(obj.pt(),jets->at(matchedJetIndex).pt());

                        }               

                        if ( obj.filterLabels()[i].find(filter2_) != std::string::npos ) {
                              histContainer_["Filter2_Pt"]->Fill(obj.pt());
                              histContainer_["Filter2_Eta"]->Fill(obj.eta());
                              histContainer_["Filter2_Phi"]->Fill(obj.phi());
                              histContainer_["Filter2_matchedJetBTag"]->Fill(-log(1-jets->at(matchedJetIndex).bDiscriminator(btagger_)));
                        }

                  if (obj.collection() == "hltAK4PFJetsCorrected::HLT") histContainer_["Total_matchedJetPt"]->Fill(jets->at(matchedJetIndex).pt()); //NOT CONVINCED
                  if (obj.collection() == "hltPFJetForBtag::HLT") histContainer_["Total_matchedJetBTag"]->Fill(-log(1-jets->at(matchedJetIndex).bDiscriminator(btagger_)));               
                  }


                  if ( hadronicleg_ == "TTBarJet30" ){
                        if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {
                              histContainer_["Filter1_Pt"]->Fill(obj.pt());
                              histContainer_["Filter1_Eta"]->Fill(obj.eta());
                              histContainer_["Filter1_Phi"]->Fill(obj.phi());
                              
                              histContainer_["Filter1_matchedJetPt"]->Fill(jets->at(matchedJetIndex).pt());
                        }
                        
                        histContainer_["Total_matchedJetPt"]->Fill(jets->at(matchedJetIndex).pt());                 
                  }


                  if ( hadronicleg_ == "TTBarJet304050" ){
                        if ( obj.filterLabels()[i].find(filter1_) != std::string::npos ) {
                              histContainer_["Filter1_Pt"]->Fill(obj.pt());
                              histContainer_["Filter1_Eta"]->Fill(obj.eta());
                              histContainer_["Filter1_Phi"]->Fill(obj.phi());
                              if (isMatched) histContainer_["Filter1_matchedJetPt"]->Fill(jets->at(matchedJetIndex).pt());
                        }
                  
                        if ( obj.filterLabels()[i].find(filter2_) != std::string::npos ) {
                              histContainer_["Filter2_Pt"]->Fill(obj.pt());
                              histContainer_["Filter2_Eta"]->Fill(obj.eta());
                              histContainer_["Filter2_Phi"]->Fill(obj.phi());
                              if (isMatched) histContainer_["Filter2_matchedJetPt"]->Fill(jets->at(matchedJetIndex).pt());
                        }
                     
                        if ( obj.filterLabels()[i].find(filter3_) != std::string::npos ) {
                              histContainer_["Filter3_Pt"]->Fill(obj.pt());
                              histContainer_["Filter3_Eta"]->Fill(obj.eta());
                              histContainer_["Filter3_Phi"]->Fill(obj.phi());
                              if (isMatched) histContainer_["Filter3_matchedJetPt"]->Fill(jets->at(matchedJetIndex).pt());
                        }
                        if (isMatched) histContainer_["Total_matchedJetPt"]->Fill(jets->at(matchedJetIndex).pt());                                         
                  }
            }
      }
}

// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyser::beginJob(){


      // INITIALISE DIRECTORIES AND HISTOGRAMS ------------------------------------------------------ //

      subDir_TrigDec = fileService->mkdir( "Trigger Decision" );
      histContainer_["CrossTriggerHist"] = subDir_TrigDec.make<TH1F>("CrossTrigger Trigger Decision", crosstrigger_.c_str(), 2, -0.5, 1.5);
      histContainer_["SingleLeptonHist"] = subDir_TrigDec.make<TH1F>("Single Lepton Trigger Decision", singleleptontrigger_.c_str(), 2, -0.5, 1.5);
      histContainer_["CrossTriggerCombinedHist"] = subDir_TrigDec.make<TH1F>("Added CrossTrigger Trigger Decision", (crosstrigger_ + " and " + singleleptontrigger_).c_str(), 2, -0.5, 1.5);

      subDir_TrigDec_TurnOnCurves = subDir_TrigDec.mkdir( "Trigger Turn On Curves" );

      subDir_Observables = fileService->mkdir( "Trigger Observables" );
      subDir_Observables_Jet = subDir_Observables.mkdir( "Jets" );
      histContainer_["CrossTrigger_Pass_JetPtHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetPt", "Pass_Pt", 50, 0, 120);
      histContainer_["CrossTrigger_Total_JetPtHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetPt", "Total_Pt", 50, 0, 120);
      histContainer_["CrossTrigger_Pass_JetEtaHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetEta", "Pass_Eta", 100, -3, 3);
      histContainer_["CrossTrigger_Total_JetEtaHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetEta", "Total_Eta", 100, -3, 3);
      histContainer_["CrossTrigger_Pass_JetCSVHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetCSV", "Pass_JetCSV", 10, 0, 1);
      histContainer_["CrossTrigger_Total_JetCSVHist"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetCSV", "Total_JetCSV", 10, 0, 1);       
      histContainer_["CrossTrigger_Pass_JetMultiplicity"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetMultiplicity", "Pass_JetMultiplicity", 10, 0, 10);
      histContainer_["CrossTrigger_Total_JetMultiplicity"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetMultiplicity", "Total_JetMultiplicity", 10, 0, 10);
      histContainer_["CrossTrigger_Pass_hltHT"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_hltHT", "Pass_hltHT", 100, 0, 700);
      histContainer_["CrossTrigger_Total_hltHT"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_hltHT", "Total_hltHT", 100, 0, 700);
      histContainer_["CrossTrigger_Pass_greatestBtag"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_GreatestBtag", "Pass_GreatestBtag", 100, 0, 8);
      histContainer_["CrossTrigger_Total_greatestBtag"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_GreatestBtag", "Total_GreatestBtag", 100, 0, 8);     
      histContainer_["CrossTrigger_Pass_forwardJetEta"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_forwardJetEta", "Pass_ForwardJetEta", 100, -3, 3);
      histContainer_["CrossTrigger_Total_forwardJetEta"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_forwardJetEta", "Total_ForwardJetEta", 100, -3, 3);

      subDir_Observables_Vertices = subDir_Observables.mkdir( "Vertices" );
      histContainer_["CrossTrigger_Pass_VertexMultiplicityHist"] = subDir_Observables_Vertices.make<TH1F>("CrossTrigger_Pass_VertexMultiplicityHist", "Pass_Vertex", 25, 0, 25);
      histContainer_["CrossTrigger_Total_VertexMultiplicityHist"] = subDir_Observables_Vertices.make<TH1F>("CrossTrigger_Total_VertexMultiplicityHist", "Total_Vertex", 25, 0, 25);

      subDir_Observables_MET = subDir_Observables.mkdir( "MET" );
      histContainer_["CrossTrigger_Pass_METHist"] = subDir_Observables_MET.make<TH1F>("CrossTrigger_Pass_MET", "Pass_Energy", 100, 0, 200);
      histContainer_["CrossTrigger_Total_METHist"] = subDir_Observables_MET.make<TH1F>("CrossTrigger_Total_MET", "Total_Energy", 100, 0, 200);

      subDir_Observables_Lepton = subDir_Observables.mkdir( "Leading_Lepton" );
      histContainer_["CrossTrigger_Pass_LeptonPtHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Pass_LeptonPt", "Pass_Pt", 50, 0, 120);
      histContainer_["CrossTrigger_Total_LeptonPtHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Total_LeptonPt", "Total_Pt", 50, 0, 120);
      histContainer_["CrossTrigger_Pass_LeptonEtaHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Pass_LeptonEta", "Pass_Eta", 100, -3, 3);
      histContainer_["CrossTrigger_Total_LeptonEtaHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Total_LeptonEta", "Total_Eta", 100, -3, 3);
      histContainer_["CrossTrigger_Pass_LeptonEnergyHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Pass_LeptonEnergy", "Pass_Energy", 100, 0, 200);
      histContainer_["CrossTrigger_Total_LeptonEnergyHist"] = subDir_Observables_Lepton.make<TH1F>("CrossTrigger_Total_LeptonEnergy", "Total_Energy", 100, 0, 200);

     if ( hadronicleg_ == "SingleTop" ){

            // Single Lepton Filter
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Trigger Object Observables" );
            histContainer_["Filter1_Pt"] = subDir_Filter1_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter1_Eta"] = subDir_Filter1_Observables.make<TH1F>("Eta", "Eta", 50, -5, 5);
            histContainer_["Filter1_Phi"] = subDir_Filter1_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter1_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt", 50, 0, 100);
            histContainer_["Total_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "matched Jet Pt", 50, 0, 100);
            distContainer_["Filter1_TriggerObject_Reco_Pt"] = subDir_Filter1_MatchedJetObservables.make<TH2F>("Matching", "Filter Object, Matched RECO Jet Pt; RECO Jet Pt (GeV); Filter Obj Jet Pt (GeV)", 50, 0, 200, 50, 0, 200);//cndjschdjkshcjkdnsk


            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );
            distContainer_["TriggerObject_Reco_Pt"] = subDir_Filter1_TurnOnCurves.make<TH2F>("Pt", "Pt", 50, 0, 200, 50, 0, 200);

            // CSV Filter
            subDir_Filter2 = fileService->mkdir( filter2_.c_str() );

            subDir_Filter2_Observables = subDir_Filter2.mkdir( "Trigger Object Observables" );
            histContainer_["Filter2_Pt"] = subDir_Filter2_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter2_Eta"] = subDir_Filter2_Observables.make<TH1F>("Eta", "Eta", 50, -5, 5);
            histContainer_["Filter2_Phi"] = subDir_Filter2_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter2_MatchedJetObservables = subDir_Filter2.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter2_matchedJetBTag"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet CSV", "Filter 2 matched Jet CSV", 50, 0, 15);
            histContainer_["Total_matchedJetBTag"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("matched Jet CSV", "matched Jet CSV", 50, 0, 15);

            subDir_Filter2_TurnOnCurves = subDir_Filter2.mkdir( "Turn On Curves" );
      }
      
      if ( hadronicleg_ == "TTBarJet30" ){

            // 30GeV Jet Filter
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Trigger Object Observables" );
            histContainer_["Filter1_Pt"] = subDir_Filter1_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter1_Eta"] = subDir_Filter1_Observables.make<TH1F>("Eta", "Eta", 50, -5, 5);
            histContainer_["Filter1_Phi"] = subDir_Filter1_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter1_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt", 50, 0, 200);
            histContainer_["Total_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "matched Jet Pt", 50, 0, 200);

            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );
      }
      
      if ( hadronicleg_ == "TTBarJet304050" ){

            // 30GeV Jet Filter
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Filter Observables" );
            histContainer_["Filter1_Pt"] = subDir_Filter1_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter1_Eta"] = subDir_Filter1_Observables.make<TH1F>("Eta", "Eta", 50, -5, 5);
            histContainer_["Filter1_Phi"] = subDir_Filter1_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter1_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt", 50, 0, 200);
            histContainer_["Total_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "matched Jet Pt", 50, 0, 200);

            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );

            // 40GeV Jet Filter
            subDir_Filter2 = fileService->mkdir( filter2_.c_str() );

            subDir_Filter2_Observables = subDir_Filter2.mkdir( "Filter Observables" );
            histContainer_["Filter2_Pt"] = subDir_Filter2_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter2_Eta"] = subDir_Filter2_Observables.make<TH1F>("Eta", "Eta", 50, -5, 5);
            histContainer_["Filter2_Phi"] = subDir_Filter2_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter2_MatchedJetObservables = subDir_Filter2.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter2_matchedJetPt"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Pt", "Filter 2 matched Jet Pt", 50, 0, 200);

            subDir_Filter2_TurnOnCurves = subDir_Filter2.mkdir( "Turn On Curves" );

            // 50GeV Jet Filter
            subDir_Filter3 = fileService->mkdir( filter3_.c_str() );

            subDir_Filter3_Observables = subDir_Filter3.mkdir( "Filter Observables" );
            histContainer_["Filter3_Pt"] = subDir_Filter3_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter3_Eta"] = subDir_Filter3_Observables.make<TH1F>("Eta", "Eta", 50, -5, 5);
            histContainer_["Filter3_Phi"] = subDir_Filter3_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter3_MatchedJetObservables = subDir_Filter3.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter3_matchedJetPt"] = subDir_Filter3_MatchedJetObservables.make<TH1F>("Filter 3 matched Jet Pt", "Filter 3 matched Jet Pt", 50, 0, 200);

            subDir_Filter3_TurnOnCurves = subDir_Filter3.mkdir( "Turn On Curves" );
      }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyser::endJob(){
      //reco vs trigger object
      // CREATE TRIGGER TURN ON CURVES -------------------------------------------------------------- //

      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMultiplicity"],histContainer_["CrossTrigger_Total_JetMultiplicity"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity"]->SetTitle("TGraphAsymmErrors Example");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity"]->Draw("ALP");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity"]->Write();

      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_HT"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_hltHT"],histContainer_["CrossTrigger_Total_hltHT"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_HT"]->SetTitle("TGraphAsymmErrors Example");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_HT"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_HT"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_HT"]->Draw("ALP");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_HT"]->Write();

      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_METHist"],histContainer_["CrossTrigger_Total_METHist"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->SetTitle("TGraphAsymmErrors Example");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->Draw("ALP");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->Write();

      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_VertexMultiplicityHist"],histContainer_["CrossTrigger_Total_VertexMultiplicityHist"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->SetTitle("TGraphAsymmErrors Example");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->Draw("ALP");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->Write();


      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetPt"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetPtHist"],histContainer_["CrossTrigger_Total_JetPtHist"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetPt"]->SetTitle("TGraphAsymmErrors Example");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetPt"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetPt"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetPt"]->Draw("ALP");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetPt"]->Write();

      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_LeptonPt"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_LeptonPtHist"],histContainer_["CrossTrigger_Total_LeptonPtHist"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_LeptonPt"]->SetTitle("TGraphAsymmErrors Example");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_LeptonPt"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_LeptonPt"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_LeptonPt"]->Draw("ALP");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_LeptonPt"]->Write();

      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetCSV"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetCSVHist"],histContainer_["CrossTrigger_Total_JetCSVHist"]);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetCSV"]->SetTitle("TGraphAsymmErrors Example");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetCSV"]->SetMarkerColor(4);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetCSV"]->SetMarkerStyle(21);
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetCSV"]->Draw("ALP");
      turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetCSV"]->Write();

      // CREATE FILTER TURN ON CURVES --------------------------------------------------------------- //

      if ( hadronicleg_ == "SingleTop" ){

            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"] = subDir_Filter1_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter1_matchedJetPt"] ,histContainer_["Total_matchedJetPt"]);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetTitle("TGraphAsymmErrors Example");
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerColor(4);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerStyle(21);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Draw("ALP");
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Write();

            turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"] = subDir_Filter2_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter2_matchedJetBTag"],histContainer_["Total_matchedJetBTag"]);
            turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetTitle("TGraphAsymmErrors Example");
            turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetMarkerColor(4);
            turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetMarkerStyle(21);
            turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->Draw("ALP");
            turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->Write();
      }

      if ( hadronicleg_ == "TTBarJet30" ){

            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"] = subDir_Filter1_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter1_matchedJetPt"],histContainer_["Total_matchedJetPt"]);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetTitle("TGraphAsymmErrors Example");
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerColor(4);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerStyle(21);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Draw("ALP");
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Write();
      }

      if ( hadronicleg_ == "TTBarJet304050" ){

            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"] = subDir_Filter1_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter1_matchedJetPt"],histContainer_["Total_matchedJetPt"]);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetTitle("TGraphAsymmErrors Example");
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerColor(4);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerStyle(21);
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Draw("ALP");
            turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Write();

            turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"] = subDir_Filter2_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter2_matchedJetPt"],histContainer_["Total_matchedJetPt"]);
            turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->SetTitle("TGraphAsymmErrors Example");
            turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->SetMarkerColor(4);
            turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->SetMarkerStyle(21);
            turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->Draw("ALP");
            turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->Write();

            turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"] = subDir_Filter3_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter3_matchedJetPt"],histContainer_["Total_matchedJetPt"]);
            turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->SetTitle("TGraphAsymmErrors Example");
            turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->SetMarkerColor(4);
            turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->SetMarkerStyle(21);
            turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->Draw("ALP");
            turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->Write();
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
