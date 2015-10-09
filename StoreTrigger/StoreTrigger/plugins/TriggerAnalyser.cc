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
#include <cmath>
#include <sstream>

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
    leptonfilter_(iConfig.getParameter <std::string> ("LeptonFilterInput")),
    btagger_(iConfig.getParameter <std::string> ("BTagger")),
    electronID_(iConfig.getParameter <std::string> ("ElectronID")),
    hadronicleg_(iConfig.getParameter <std::string> ("HadronicLeg")),
    leptonicleg_(iConfig.getParameter <std::string> ("LeptonicLeg")),
    leptontype_(iConfig.getParameter <std::string> ("LeptonType")){
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

      // using namespace edm;

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


      // EVENT SELECTION AND INDEX STORING ---------------------------------------------------------- //

      // Perform Event Selection
      passEventSelection = false;

      tightRecoElectrons.clear();
      tightRecoMuons.clear();
      cleanedRecoJets.clear();
      cleanedRecoBJets.clear();

      // Tight Electrons
      if (electrons->size() != 0){
            for (auto lepton = electrons->begin(); lepton != electrons->end(); ++lepton){
                  if (leptontype_ == "Mu20") ptcut=30;//change pt cuts depening on triggers
                  if (leptontype_ == "Mu24") ptcut=30;
                  if (leptontype_ == "Ele27") ptcut=30;
                  if (leptontype_ == "Ele32") ptcut=35;

                  if (lepton->pt() > ptcut && abs(lepton->eta()) < 2.1 && lepton->electronID(electronID_)){
                        tightRecoElectrons.push_back(*lepton);
                  }
            }
      }

      // Tight Isolated Muons
      if (vertices->size() != 0){
            auto vertex = vertices->at(0);
            if (vertex.isValid()) {
                  if (muons->size() != 0){
                        for (auto lepton = muons->begin(); lepton != muons->end(); ++lepton){

                              if (leptontype_ == "Mu20") ptcut=22;
                              if (leptontype_ == "Mu24") ptcut=27;
                              if (leptontype_ == "Ele27") ptcut=22;
                              if (leptontype_ == "Ele32") ptcut=22;

                              if ( lepton->pt() > ptcut && abs(lepton->eta() && lepton->isTightMuon(vertex) && isIsolated(*lepton)) < 2.1 ){
                                    tightRecoMuons.push_back(*lepton);
                              }
                        }
                  }
            }
      }

      // Cleaned Jets
      if (jets->size() != 0){
            for( auto jet = jets->begin(); jet != jets->end(); ++jet ){
                  isMatchedToLepton = false;
                  if (!isGoodJet(*jet)) continue;
      
                  if ( tightRecoElectrons.size() != 0 ){
                        if ( reco::deltaR(tightRecoElectrons[0], *jet) < 0.3 ) isMatchedToLepton = true;
                  }
            
                  if ( tightRecoMuons.size() != 0 ){
                        if ( reco::deltaR(tightRecoMuons[0], *jet) < 0.3) isMatchedToLepton = true;
                  }

                  if (isMatchedToLepton) continue;
                  cleanedRecoJets.push_back(*jet);
                  
                  if (jet->bDiscriminator(btagger_) < 0.890 ) continue;// 2 med bjet value (https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X50ns#Supported_Algorithms_and_Operati)
                  cleanedRecoBJets.push_back(*jet);
            }
      }

      // Perform an offline event selection
      if (leptonicleg_ == "Ele"){
            if (hadronicleg_ == "SingleTop"){
                  if (cleanedRecoJets.size() >= 2 && cleanedRecoBJets.size() >= 1 && tightRecoMuons.size() == 0 && tightRecoElectrons.size() == 1) passEventSelection = true;
            }
             if (hadronicleg_ == "TTBarJet30" || hadronicleg_ == "TTBarJet304050"){
                  if (cleanedRecoJets.size() >= 4 && cleanedRecoBJets.size() >= 2 && tightRecoMuons.size() == 0 && tightRecoElectrons.size() == 1) passEventSelection = true;
            }
      }
      if (leptonicleg_ == "Mu"){
            if (hadronicleg_ == "SingleTop"){
                  if (cleanedRecoJets.size() >= 2 && cleanedRecoBJets.size() >= 1 && tightRecoMuons.size() == 1 && tightRecoElectrons.size() == 0) passEventSelection = true;
            }
             if (hadronicleg_ == "TTBarJet30" || hadronicleg_ == "TTBarJet304050"){
                  if (cleanedRecoJets.size() >= 4 && cleanedRecoBJets.size() >= 2 && tightRecoMuons.size() == 1 && tightRecoElectrons.size() == 0) passEventSelection = true;
            }
      }


      // TRIGGER NAMES ------------------------------------------------------------------------------------ //

      SL_TrigDecision = false;
      X_TrigDecision = false;

      // Find and store indices of the input triggers
      for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) {

            // std::cout << "Trigger " << i << " in Menu : " << TrigNames.triggerName(i) << std::endl;
            
            if ( TrigNames.triggerName(i).find(singleleptontrigger_) != std::string::npos ) {
                  if(triggerResults->accept(i) ) SL_TrigDecision = true;
            }
            if ( TrigNames.triggerName(i).find(crosstrigger_) != std::string::npos ) {
                  if(triggerResults->accept(i) ) X_TrigDecision = true;
            }
      }

      TypeOfEvent=0;
      histContainer_["TypeOfEvent"]->Fill(TypeOfEvent); // Number of Events


      // DIFFERENTIAL EFFICIENCIES ------------------------------------------------------------------------ //
      if (passEventSelection){

            TypeOfEvent=1;
            histContainer_["TypeOfEvent"]->Fill(TypeOfEvent);

            if (SL_TrigDecision == true){
                  TypeOfEvent=2;
                  histContainer_["TypeOfEvent"]->Fill(TypeOfEvent);
            }

            if (X_TrigDecision == true){
                  TypeOfEvent=3;
                  histContainer_["TypeOfEvent"]->Fill(TypeOfEvent);
            }            

            // VERTICES ----------------------------------------------------------------------------------- //

            vertexMult = 0;
            for( auto vertex = vertices->begin(); vertex != vertices->end(); ++vertex ){ 
                  if (vertex->isValid()) vertexMult += 1;
            }

            if (X_TrigDecision==true){
                  histContainer_["CrossTrigger_Pass_VertexMultHist"]->Fill(vertexMult);
            }
            histContainer_["CrossTrigger_Total_VertexMultHist"]->Fill(vertexMult);


            // JETS --------------------------------------------------------------------------------------- //

            // Select jets and store distributions
            jetMult_20 = jetMult_30 = jetMult_40 = jetMult_50 = 0;
            forwardjeteta = jetCSV = 0;
            HT = ST = 0;
            
            for (uint Jet = 0; Jet < cleanedRecoJets.size(); ++Jet){

                  // Number of jets in an event and HT depending on jet pt cut
                  if ( cleanedRecoJets[Jet].pt()>20) ++jetMult_20;
                  if ( cleanedRecoJets[Jet].pt()>30) ++jetMult_30;
                  if ( cleanedRecoJets[Jet].pt()>40) ++jetMult_40;
                  if ( cleanedRecoJets[Jet].pt()>50) ++jetMult_50;
                  
                  // HT is the sum of cleaned central jets over 20GeV (No MET or Leptons) for events that pass selection
                  HT += cleanedRecoJets[Jet].pt();

                  // Most forward jet has highest eta
                  if (std::abs(cleanedRecoJets[Jet].eta()) >= forwardjeteta) forwardjeteta = std::abs(cleanedRecoJets[Jet].eta());

                  // Jet with greatest btag value
                  if (cleanedRecoJets[Jet].bDiscriminator(btagger_) >= jetCSV) jetCSV = cleanedRecoJets[Jet].bDiscriminator(btagger_);

                  // Store Jet Pt, Eta, Higest CSV, Global_HT and Mult in histograms 
                  if (X_TrigDecision==true){
                        histContainer_["CrossTrigger_Pass_JetPtHist"]->Fill(cleanedRecoJets[Jet].pt());
                        histContainer_["CrossTrigger_Pass_JetEtaHist"]->Fill(cleanedRecoJets[Jet].eta());
                        histContainer_["CrossTrigger_Pass_JetPhiHist"]->Fill(cleanedRecoJets[Jet].phi());
                        histContainer_["CrossTrigger_Pass_JetCSVHist"]->Fill(cleanedRecoJets[Jet].bDiscriminator(btagger_));
                  }
                  histContainer_["CrossTrigger_Total_JetPtHist"]->Fill(cleanedRecoJets[Jet].pt());
                  histContainer_["CrossTrigger_Total_JetEtaHist"]->Fill(cleanedRecoJets[Jet].eta());
                  histContainer_["CrossTrigger_Total_JetPhiHist"]->Fill(cleanedRecoJets[Jet].phi());
                  histContainer_["CrossTrigger_Total_JetCSVHist"]->Fill(cleanedRecoJets[Jet].bDiscriminator(btagger_));
            }

            if(X_TrigDecision==true){
                  histContainer_["CrossTrigger_Pass_JetMult_20"]->Fill(jetMult_20);
                  histContainer_["CrossTrigger_Pass_JetMult_30"]->Fill(jetMult_30);
                  histContainer_["CrossTrigger_Pass_JetMult_40"]->Fill(jetMult_40);
                  histContainer_["CrossTrigger_Pass_JetMult_50"]->Fill(jetMult_50);
                  histContainer_["CrossTrigger_Pass_greatestBtag"]->Fill(jetCSV);
                  histContainer_["CrossTrigger_Pass_forwardJetEta"]->Fill(forwardjeteta);
                  histContainer_["CrossTrigger_Pass_HT"]->Fill(HT);
            }
            histContainer_["CrossTrigger_Total_JetMult_20"]->Fill(jetMult_20);
            histContainer_["CrossTrigger_Total_JetMult_30"]->Fill(jetMult_30);
            histContainer_["CrossTrigger_Total_JetMult_40"]->Fill(jetMult_40);
            histContainer_["CrossTrigger_Total_JetMult_50"]->Fill(jetMult_50);
            histContainer_["CrossTrigger_Total_greatestBtag"]->Fill(jetCSV);
            histContainer_["CrossTrigger_Total_forwardJetEta"]->Fill(forwardjeteta);
            histContainer_["CrossTrigger_Total_HT"]->Fill(HT);  

            ST += HT;


            // MET ---------------------------------------------------------------------------------------- //
      
            // Store MET of event
            for( auto met = mets->begin(); met != mets->end(); ++met ){ 
                  histContainer_["CrossTrigger_Total_METHist"]->Fill(met->energy());
                  if(X_TrigDecision==true){
                        histContainer_["CrossTrigger_Pass_METHist"]->Fill(met->energy());
                  }
                  ST += met->energy();//MET = MPT
            }


            // LEADING LEPTONS ---------------------------------------------------------------------------- //
     
            // Store the leading lepton pt, eta and energy distributions in histograms
            if (leptonicleg_ == "Ele"){            
                  for (uint El = 0; El < tightRecoElectrons.size(); ++El){

                        histContainer_["CrossTrigger_Total_LeptonPtHist"]->Fill(tightRecoElectrons[El].pt());
                        histContainer_["CrossTrigger_Total_LeptonEtaHist"]->Fill(tightRecoElectrons[El].eta());
                        histContainer_["CrossTrigger_Total_LeptonPhiHist"]->Fill(tightRecoElectrons[El].phi());
                        histContainer_["CrossTrigger_Total_LeptonEnergyHist"]->Fill(tightRecoElectrons[El].energy());
                        if(X_TrigDecision==true){
                              histContainer_["CrossTrigger_Pass_LeptonPtHist"]->Fill(tightRecoElectrons[El].pt());
                              histContainer_["CrossTrigger_Pass_LeptonEtaHist"]->Fill(tightRecoElectrons[El].eta());
                              histContainer_["CrossTrigger_Pass_LeptonPhiHist"]->Fill(tightRecoElectrons[El].phi());
                              histContainer_["CrossTrigger_Pass_LeptonEnergyHist"]->Fill(tightRecoElectrons[El].energy());
                        }

                        ST += tightRecoElectrons[El].pt();
                  }
            }

            if ( leptonicleg_ == "Mu" ){
                  for (uint Mu = 0; Mu < tightRecoMuons.size(); ++Mu){

                        histContainer_["CrossTrigger_Total_LeptonPtHist"]->Fill(tightRecoMuons[Mu].pt());
                        histContainer_["CrossTrigger_Total_LeptonEtaHist"]->Fill(tightRecoMuons[Mu].eta());
                        histContainer_["CrossTrigger_Total_LeptonPhiHist"]->Fill(tightRecoMuons[Mu].phi());
                        histContainer_["CrossTrigger_Total_LeptonEnergyHist"]->Fill(tightRecoMuons[Mu].energy());
                        if(X_TrigDecision==true){
                              histContainer_["CrossTrigger_Pass_LeptonPtHist"]->Fill(tightRecoMuons[Mu].pt());
                              histContainer_["CrossTrigger_Pass_LeptonEtaHist"]->Fill(tightRecoMuons[Mu].eta());
                              histContainer_["CrossTrigger_Pass_LeptonPhiHist"]->Fill(tightRecoMuons[Mu].phi());
                              histContainer_["CrossTrigger_Pass_LeptonEnergyHist"]->Fill(tightRecoMuons[Mu].energy());
                        }

                        ST += tightRecoMuons[Mu].pt();
                  }
            } 
      }


      // ST ----------------------------------------------------------------------------------------- //

      // std::cout << "ST : " << ST << std::endl;
      histContainer_["CrossTrigger_Total_ST"]->Fill(ST);
      if(X_TrigDecision==true){
            histContainer_["CrossTrigger_Pass_ST"]->Fill(ST);
      }


      // ############################################################################################ //
      // ############################################################################################ //
      // ############################################################################################ //
      // ############################################################################################ //
      // ############################################################################################ //



      // FILTERS: TRIGGER OBJECTS ------------------------------------------------------------------- //

      // if (SL_TrigDecision){
      if (X_TrigDecision){

            hltJets.clear();
            hltBJets.clear();
            hltLeadingLeptons.clear();

            int N_Filter1_Tags = 0;
            int N_Filter2_Tags = 0;
            int N_Filter3_Tags = 0;

            bool passFilter1 = false;
            bool passFilter2 = false;
            bool passFilter3 = false;

            // Create hlt object collections
            for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                  if (obj.collection() == "hltPFJetForBtag::HLT") hltBJets.push_back(obj);
                  if (obj.collection() == "hltAK4PFJetsCorrected::HLT") hltJets.push_back(obj);
                  if (obj.collection() == "hltL3MuonCandidates::HLT" && obj.hasFilterLabel(leptonfilter_)) hltLeadingLeptons.push_back(obj);
                  if (obj.collection() == "hltEgammaCandidates::HLT" && obj.hasFilterLabel(leptonfilter_)) hltLeadingLeptons.push_back(obj);

                  if (obj.hasFilterLabel(filter1_)) ++N_Filter1_Tags;
                  if (obj.hasFilterLabel(filter2_)) ++N_Filter2_Tags;
                  if (obj.hasFilterLabel(filter3_)) ++N_Filter3_Tags;
            }
            
            if ((hadronicleg_ == "TTBarJet304050" || hadronicleg_ == "TTBarJet30") && (N_Filter1_Tags >= 3)) passFilter1 = true;
            if ((hadronicleg_ == "TTBarJet304050") && (N_Filter2_Tags >= 2)) passFilter2 = true;
            if ((hadronicleg_ == "TTBarJet304050") && (N_Filter3_Tags >= 1)) passFilter3 = true;

            // if (hadronicleg_ == "TTBarJet304050" ) std::cout << "Event has : " << N_Filter1_Tags << " tricentral>30 tags, " << N_Filter2_Tags << " dicentral>40 tags and " << N_Filter3_Tags << " central>50 tags." << std::endl;
            // if (hadronicleg_ == "TTBarJet304050" ) std::cout << "Filter 1 : " << passFilter1 << ", Filter 2 : " << passFilter2 << ", Filter 3 : " << passFilter3 << std::endl;


            ///////////////////////////////////
            // Print number of Reco/HLT jets //
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // std::cout << "-----------------" << std::endl;                                                                                   //
            // std::cout << "Number of hlt Jets : " << hltJets.size() << std::endl;                                                             //
            // std::cout << "Number of hlt BJets : " << hltBJets.size() << std::endl;                                                           //
            // std::cout << "Number of hlt leading leptons : " << hltLeadingLeptons.size() << std::endl;                                        //
            //                                                                                                                                  //
            // for( auto jet = jets->begin(); jet != jets->end(); ++jet ){                                                                      //
            //       std::cout << "RECO pt : " << jet->pt() << ", eta : " << jet->eta() << std::endl;                                           //
            // }                                                                                                                                //
            // for (uint x = 0; x < hltJets.size(); ++x){                                                                                       //
            //       std::cout << "hlt pt : " << hltJets[x].pt() << ", eta : " << hltJets[x].eta() << std::endl;                                //
            // }                                                                                                                                //
            // for (uint x = 0; x < hltLeadingLeptons.size(); ++x){                                                                             //
            //       std::cout << "lep pt : " << hltLeadingLeptons[x].pt() << ", eta : " << hltLeadingLeptons[x].eta() << std::endl;            //
            // }                                                                                                                                //
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // Clean the RECO jets wrt hlt objects
            cleanedJets.clear();
            for( auto jet = jets->begin(); jet != jets->end(); ++jet ){
                  if (!isGoodJet(*jet)) continue;

                  bool isJetMatched = false;
                  bool isLeptonMatched = false;

                  for (uint hltjet = 0; hltjet < hltJets.size(); ++hltjet){
                        if (reco::deltaR(hltJets[hltjet], *jet) < 0.3) isJetMatched = true;
                  }
                  for (uint hltBjet = 0; hltBjet < hltBJets.size(); ++hltBjet){
                        if (reco::deltaR(hltBJets[hltBjet], *jet) < 0.3) isJetMatched = true;
                  }

                  if (reco::deltaR(hltLeadingLeptons[0], *jet) < 0.3) isLeptonMatched = true;//Clean vs leading lepton only - still want to consider secondary leptons as jets

                  if (!isJetMatched || isLeptonMatched) continue;

                  cleanedJets.push_back(*jet);
                  // if ((hadronicleg_ == "TTBarJet304050" || hadronicleg_ == "TTBarJet30") && (N_Filter1_Tags == 2 || N_Filter1_Tags == 1 ) ) continue;

                  // Fill clean central jets that are matched to a hlt jet object and not matched to a hlt lepton object into histograms
                  histContainer_["Total_RECO_JetPt"]->Fill(jet->pt());
                  histContainer_["Total_RECO_JetBTag"]->Fill(-log10(1-jet->bDiscriminator(btagger_)));
                  histContainer_["Total_RECO_JetCSV"]->Fill(jet->bDiscriminator(btagger_));
            }

            //////////////////////////////////////////////////
            // Compare number of cleaned jets to total jets //
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // std::cout << "Jet size : " << jets->size() << std::endl;                                                                         //
            // std::cout << "Cleaned jet size : " << cleanedJets.size() << std::endl;                                                           //
            // for (uint x = 0; x < cleanedJets.size(); ++x){                                                                                   //
            //       std::cout << "Cleaned jet Pt : " << cleanedJets[x].pt() << std::endl;                                                      //
            // }                                                                                                                                //
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if (cleanedJets.size() != 0){

                  if (hadronicleg_ == "SingleTop" ){
                        for (uint recojet = 0; recojet < cleanedJets.size(); ++recojet){

                              for (uint hltjet = 0; hltjet < hltJets.size(); ++hltjet){
                                    if (reco::deltaR(cleanedJets[recojet], hltJets[hltjet]) > 0.3) continue;
                                    if ( hltJets[hltjet].hasFilterLabel(filter1_)){

                                          histContainer_["Filter1_Pt"]->Fill(hltJets[hltjet].pt());
                                          histContainer_["Filter1_Eta"]->Fill(hltJets[hltjet].eta());
                                          histContainer_["Filter1_Phi"]->Fill(hltJets[hltjet].phi());

                                          histContainer_["Filter1_matchedJetPt"]->Fill(cleanedJets[recojet].pt());
                                          histContainer_["Filter1_matchedJetEta"]->Fill(cleanedJets[recojet].eta());
                                          histContainer_["Filter1_matchedJetPhi"]->Fill(cleanedJets[recojet].phi());
                                          break;
                                    }
                              }

                              for (uint hltBjet = 0; hltBjet < hltBJets.size(); ++hltBjet){
                                    if ( reco::deltaR(cleanedJets[recojet], hltBJets[hltBjet]) > 0.3) continue;
                                    if ( hltBJets[hltBjet].hasFilterLabel(filter2_)){

                                          histContainer_["Filter2_Pt"]->Fill(hltBJets[hltBjet].pt());
                                          histContainer_["Filter2_Eta"]->Fill(hltBJets[hltBjet].eta());
                                          histContainer_["Filter2_Phi"]->Fill(hltBJets[hltBjet].phi());

                                          histContainer_["Filter2_matchedJetPt"]->Fill(cleanedJets[recojet].pt());
                                          histContainer_["Filter2_matchedJetEta"]->Fill(cleanedJets[recojet].eta());
                                          histContainer_["Filter2_matchedJetPhi"]->Fill(cleanedJets[recojet].phi());
                                    
                                          histContainer_["Filter2_matchedJetBTag"]->Fill(-log10(1-cleanedJets[recojet].bDiscriminator(btagger_)));
                                          histContainer_["Filter2_matchedJetCSV"]->Fill(cleanedJets[recojet].bDiscriminator(btagger_));
                                          break;
                                    }
                              }
                        }
                  }

                  if (hadronicleg_ == "TTBarJet30" ){
                        for (uint recojet = 0; recojet < cleanedJets.size(); ++recojet){
                              for (uint hltjet = 0; hltjet < hltJets.size(); ++hltjet){
                                    if ( reco::deltaR(cleanedJets[recojet], hltJets[hltjet]) > 0.3 ) continue;
                                    if ( !passFilter1 ) continue;
                                    if ( hltJets[hltjet].hasFilterLabel(filter1_) ){

                                          histContainer_["Filter1_Pt"]->Fill(hltJets[hltjet].pt());
                                          histContainer_["Filter1_Eta"]->Fill(hltJets[hltjet].eta());
                                          histContainer_["Filter1_Phi"]->Fill(hltJets[hltjet].phi());

                                          histContainer_["Filter1_matchedJetPt"]->Fill(cleanedJets[recojet].pt());
                                          histContainer_["Filter1_matchedJetEta"]->Fill(cleanedJets[recojet].eta());
                                          histContainer_["Filter1_matchedJetPhi"]->Fill(cleanedJets[recojet].phi());
                                          break;
                                    }
                              }
                        }
                  }

                  if (hadronicleg_ == "TTBarJet304050" ){

                        bool pF1, pF2, pF3;
                        pF1 = pF2 = pF3 = false;
                        float minF1dR2, minF2dR2, minF3dR2;
                        minF1dR2 = minF2dR2 = minF3dR2 = 9999;

                        for (uint recojet = 0; recojet < cleanedJets.size(); ++recojet){


                              histContainer_["Total_RECO_Filter1Background_JetPt"]->Fill(cleanedJets[recojet].pt());
                              if (passFilter1) histContainer_["Total_RECO_Filter2Background_JetPt"]->Fill(cleanedJets[recojet].pt());
                              if (passFilter2) histContainer_["Total_RECO_Filter3Background_JetPt"]->Fill(cleanedJets[recojet].pt());

                              if (cleanedJets[recojet].pt() > 100 && abs(cleanedJets[recojet].eta() < 2.0)){
                                    for (uint hltjet = 0; hltjet < hltJets.size(); ++hltjet){
                                          float dR2 = reco::deltaR(cleanedJets[recojet], hltJets[hltjet]);
                                          if ( hltJets[hltjet].hasFilterLabel(filter1_)){
                                                if (dR2 < minF1dR2){
                                                      minF1dR2 = dR2;
                                                      pF1 = true;
                                                }
                                          }
                                          if ( hltJets[hltjet].hasFilterLabel(filter2_)){
                                                if (dR2 < minF2dR2){
                                                      minF2dR2 = dR2;
                                                      pF2 = true;
                                                }
                                          }        
                                          if ( hltJets[hltjet].hasFilterLabel(filter3_)){
                                                if (dR2 < minF3dR2){
                                                      minF3dR2 = dR2;
                                                      pF3 = true;
                                                }
                                          }
                                    }
                              }

                              for (uint hltjet = 0; hltjet < hltJets.size(); ++hltjet){

                                    if (reco::deltaR(cleanedJets[recojet], hltJets[hltjet]) > 0.3) continue;

                                    if (!passFilter1) continue;
                                    distContainer_["Filter1_hlt_vs_Reco_Pt"]->Fill(hltJets[hltjet].pt(),cleanedJets[recojet].pt());
                                    if ( hltJets[hltjet].hasFilterLabel(filter1_)){

                                          histContainer_["Filter1_Pt"]->Fill(hltJets[hltjet].pt());
                                          histContainer_["Filter1_Eta"]->Fill(hltJets[hltjet].eta());
                                          histContainer_["Filter1_Phi"]->Fill(hltJets[hltjet].phi());

                                          histContainer_["Filter1_matchedJetPt"]->Fill(cleanedJets[recojet].pt());
                                          histContainer_["Filter1_matchedJetEta"]->Fill(cleanedJets[recojet].eta());
                                          histContainer_["Filter1_matchedJetPhi"]->Fill(cleanedJets[recojet].phi());
                                          
                                          if (!passFilter2) continue;
                                          if ( hltJets[hltjet].hasFilterLabel(filter2_)){

                                                histContainer_["Filter2_Pt"]->Fill(hltJets[hltjet].pt());
                                                histContainer_["Filter2_Eta"]->Fill(hltJets[hltjet].eta());
                                                histContainer_["Filter2_Phi"]->Fill(hltJets[hltjet].phi());

                                                histContainer_["Filter2_matchedJetPt"]->Fill(cleanedJets[recojet].pt());
                                                histContainer_["Filter2_matchedJetEta"]->Fill(cleanedJets[recojet].eta());
                                                histContainer_["Filter2_matchedJetPhi"]->Fill(cleanedJets[recojet].phi());

                                                
                                                if (!passFilter3) continue;
                                                if ( hltJets[hltjet].hasFilterLabel(filter3_)){

                                                      histContainer_["Filter3_Pt"]->Fill(hltJets[hltjet].pt());
                                                      histContainer_["Filter3_Eta"]->Fill(hltJets[hltjet].eta());
                                                      histContainer_["Filter3_Phi"]->Fill(hltJets[hltjet].phi());

                                                      histContainer_["Filter3_matchedJetPt"]->Fill(cleanedJets[recojet].pt());
                                                      histContainer_["Filter3_matchedJetEta"]->Fill(cleanedJets[recojet].eta());
                                                      histContainer_["Filter3_matchedJetPhi"]->Fill(cleanedJets[recojet].phi());
                                          

                                                }
                                          }
                                    break;
                                    }
                              }
                        }
                  
                        if (pF1){
                              if (minF1dR2 > 0.3) minF1dR2 = 0.305; // fill overflow bin
                              histContainer_["Filter1_dR"]->Fill(minF1dR2);
                        } 
                        if (pF2){
                              if (minF2dR2 > 0.3) minF2dR2 = 0.305;
                              histContainer_["Filter2_dR"]->Fill(minF2dR2);
                        }                         
                        if (pF3){
                              if (minF3dR2 > 0.3) minF3dR2 = 0.305;
                              histContainer_["Filter3_dR"]->Fill(minF3dR2);
                        }                   
                  }
            }

            //////////////////////////////////
            // Print hlt object information //
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //                                                                                                                                  //
            // if  (hadronicleg_ == "TTBarJet304050" ){                                                                                         //
            //       std::cout << std::endl;                                                                                                    //
            //       std::cout << "________________________________________" << std::endl;                                                      //
            //       std::cout << "_______________NEW EVENT________________" << std::endl;                                                      //
            //                                                                                                                                  //
            //       std::cout << "Event has : " << hltLeadingLeptons.size() << " leptons passing the filter" << std::endl;                       //                                                                                                                     //
            //       for (uint hltlep = 0; hltlep < hltLeadingLeptons.size(); ++hltlep){
            //             std::cout << "Lepton Pt : " << hltLeadingLeptons[hltlep].pt() << ", Lepton Eta : " << hltLeadingLeptons[hltlep].eta() 
            //             << ", Lepton Phi : " << hltLeadingLeptons[hltlep].phi() << std::endl;          //
            //       }          
            //                                                                                                      //
            //       std::cout << "Event has : " << cleanedJets.size() << " cleaned RECO Jets against all leptons: " << std::endl;            //
            //       for (uint recojet = 0; recojet<cleanedJets.size(); ++recojet){
            //             std::cout << "RECOjet Pt : " << cleanedJets[recojet].pt() << ", RECOjet Eta : " << cleanedJets[recojet].eta() 
            //             << ", RECOjet Phi : " << cleanedJets[recojet].phi() << std::endl;    
            //       }                                                                                              //
            //                                                                                                                                  //
            //       int count1 = 0, count2 = 0, count3 = 0, count4 = 0, count5 = 0, count6 = 0;                                                //
            //                                                                                                                                  //
            //       for (pat::TriggerObjectStandAlone obj : *triggerObjects) {                                                                 //
            //                                                                                                                                  //
            //             isJetCollection = false;                                                                                             //
            //             if (obj.collection() == "hltAK4PFJetsCorrected::HLT") isJetCollection = true;                                        //
            //                                                                                                                                  //
            //             // Only list for first run through of jets.                                                                          //
            //             // if (isJetCollection && std::abs(obj.eta()) < 2.6 && obj.pt() > 29.9){  
            //             if (isJetCollection){                                                //
            //                                                                                                                                  //
            //                   std::cout << "__________NEW HLT FILTER OBJ____________" << std::endl;                                          //
            //                                                                                                                                  //
            //                   std::cout << "Obj Pt is : " << obj.pt() << ", Obj Eta is : " << obj.eta() << ", Obj Phi is : " << obj.phi() << std::endl;                        //
            //                   std::cout << "Obj Collection is : " << obj.collection() << std::endl;                                          //
            //                   std::cout << "Obj Filter Ids are : ";                                                                          //
            //                   for (unsigned int i = 0; i < obj.filterIds().size(); ++i) std::cout << " " << obj.filterIds()[i];              //
            //                   std::cout << std::endl;                                                                                        //
            //                                                                                                                                  //
            //                   if (obj.pt() > 30) ++count1;                                                                                   //
            //                   if (obj.pt() > 40) ++count2;                                                                                   //
            //                   if (obj.pt() > 50) ++count3;                                                                                   //
            //                                                                                                                                  //
            //                   if ( obj.hasFilterLabel(filter1_) ) std::cout << "Found filter : " << filter1_                                 //
            // << "filter" << std::endl, ++count4;                                                                                              //
            //                   if ( obj.hasFilterLabel(filter2_) ) std::cout << "Found filter : " << filter2_                                 //
            // << "filter" << std::endl, ++count5;                                                                                              //
            //                   if ( obj.hasFilterLabel(filter3_) ) std::cout << "Found filter : " << filter3_                                 //
            // << "filter" << std::endl, ++count6;                                                                                              //
            //                                                                                                                                  //
            //             }                                                                                                                    //
            //       }                                                                                                                          //
            //                                                                                                                                  //
            //       std::cout << "________________________________________" << std::endl;                                                      //
            //       std::cout << "Number of hlt obj jets above 30 GeV : " << count1 << std::endl;                                              //
            //       std::cout << "Number of 30GeV filter instances found : " << count4 << std::endl;                                           //
            //       std::cout << "Number of hlt obj jets above 40 GeV : " << count2 << std::endl;                                              //
            //       std::cout << "Number of 40GeV filter instances found : " << count5 << std::endl;                                           //
            //       std::cout << "Number of hlt obj jets above 50 GeV : " << count3 << std::endl;                                              //
            //       std::cout << "Number of 50GeV filter instances found : " << count6 << std::endl;                                           //
            // }   
            // if (hltLeadingLeptons.size() > 1) system("read");                                                                                                                             //
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      }


}



bool TriggerAnalyser::isGoodJet(const pat::Jet& jet) const {

      bool passesJetID(false);
      bool passesPtAndEta(jet.pt() > 20 && fabs(jet.eta()) < 2.6);

      bool passNHF = jet.neutralHadronEnergyFraction() < 0.99;
      bool passNEMF = jet.neutralEmEnergyFraction() < 0.99;
      // bool passMUF = jet.muonEnergyFraction() < 0.8;
      bool passNumConst = (jet.chargedMultiplicity()+jet.neutralMultiplicity()) > 1;

      passesJetID = passNHF && passNEMF && passNumConst; // && passMUF
      
      bool passCEMF = false;
      bool passCHF = false;
      bool passCHM = false;
      bool passNumNeutConst = false;

      if (fabs(jet.eta()) < 2.4) {
            passCEMF = jet.chargedEmEnergyFraction() < 0.99;
            passCHF = jet.chargedHadronEnergyFraction() > 0;
            passCHM = jet.chargedMultiplicity() > 0;
            passesJetID = passNHF && passNEMF && passNumConst && passCHF && passCHM && passCEMF; // && passMUF
      
      }
      if (fabs(jet.eta()) > 3.0) {
            passNEMF = jet.neutralEmEnergyFraction() < 0.90;
            passNumNeutConst = jet.neutralMultiplicity() > 10;
            passesJetID = passNHF && passNEMF && passNumConst && passCHF && passCHM && passCEMF && passNumNeutConst;      
      }

      ///////////////////////////////
      // Whats wrong with the jets //
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // if (!passesPtAndEta) std::cout << "Does not pass jet pt/eta requirements" << std::endl;                                                                        //
      // if (!passesJetID){                                                                                                                                             //   
      //       std::cout << "Does not pass jet ID because : ";                                                                                                          //
      //       if (!passNHF) std::cout << "Has not passed Neutral Hadron Energy Fraction < 0.99 : " << jet.neutralHadronEnergyFraction() << std::endl;                  //
      //       if (!passNEMF) std::cout << "Has not passed Neutral EM Energy Fraction < 0.99 : " << jet.neutralEmEnergyFraction() << std::endl;                         //
      //       if (!passMUF) std::cout << "Has not passed muon Energy Fraction < 0.8 : " << jet.muonEnergyFraction() << std::endl;                                      //
      //       if (!passNumConst) std::cout << "Has not passed Number of Constituents > 1: " << jet.chargedMult()+jet.neutralMult() << std::endl;       //
      //       if (!passCHF) std::cout << "Has not passed Charged Hadron Energy Fraction > 0. : " << jet.chargedHadronEnergyFraction() << std::endl;                    //
      //       if (!passCHM) std::cout << "Has not passed Charged Mult > 0 : " << jet.chargedMult() << std::endl;                                       //
      //       if (!passCEMF) std::cout << "Has not passed Charged EM Energy Fraction < 0.99 : " << jet.chargedEmEnergyFraction() << std::endl;                         //
      // }                                                                                                                                                              //
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      return passesPtAndEta && passesJetID;    
}

// is the muon isolated?
bool TriggerAnalyser::isIsolated(const pat::Muon& muon) const {
      bool passesIsolation(false);
      float isolation = ( (muon.chargedHadronIso() + std::max( muon.neutralHadronIso() + muon.photonIso() - 0.5 * muon.puChargedHadronIso(), 0.) ) / muon.pt() );
      if (isolation < 0.12) passesIsolation = true;
      return passesIsolation;    
}

// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyser::beginJob(){
      // a= b = 0;
      // INITIALISE DIRECTORIES AND HISTOGRAMS ------------------------------------------------------ //

      subDir_TrigDec = fileService->mkdir( "Trigger Decision" );

      histContainer_["TypeOfEvent"] = subDir_TrigDec.make<TH1F>("TypeOfEvent", "Type of Event; Total/Selection/SingleLepton/Cross; Number of Events", 4, -0.5, 3.5);
      histContainer_["TypeOfEvent"]->GetXaxis()->SetBinLabel(1, "Total Events");
      histContainer_["TypeOfEvent"]->GetXaxis()->SetBinLabel(3, "Pass Single Lepton Trigger");
      histContainer_["TypeOfEvent"]->GetXaxis()->SetBinLabel(4, "Pass Cross Trigger");
      histContainer_["TypeOfEvent"]->GetXaxis()->SetBinLabel(2, "Pass Event Selection");
      histContainer_["CondEff"] = subDir_TrigDec.make<TH1F>("ConditionalEff", "Conditional Efficiency; Efficiency; ", 100, 0, 1);

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
      histContainer_["CrossTrigger_Pass_JetMult_20"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetMult_20", "Pass_JetMult_20; RECO Jet Mult; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Total_JetMult_20"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetMult_20", "Total_JetMult_20; RECO Jet Mult; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Pass_JetMult_30"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetMult_30", "Pass_JetMult_30; RECO Jet Mult; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Total_JetMult_30"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetMult_30", "Total_JetMult_30; RECO Jet Mult; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Pass_JetMult_40"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetMult_40", "Pass_JetMult_40; RECO Jet Mult; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Total_JetMult_40"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetMult_40", "Total_JetMult_40; RECO Jet Mult; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Pass_JetMult_50"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_JetMult_50", "Pass_JetMult_50; RECO Jet Mult; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Total_JetMult_50"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_JetMult_50", "Total_JetMult_50; RECO Jet Mult; Number of Events", 15, 0, 15);
      histContainer_["CrossTrigger_Pass_HT"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_HT", "Pass_HT; RECO Jet HT (GeV); Number of Events", 50, 0, 700);
      histContainer_["CrossTrigger_Total_HT"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_HT", "Total_HT; RECO Jet HT (GeV); Number of Events", 50, 0, 700);
      histContainer_["CrossTrigger_Pass_ST"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_ST", "Pass_ST; RECO Jet ST (GeV); Number of Events", 50, 0, 1000);
      histContainer_["CrossTrigger_Total_ST"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_ST", "Total_ST; RECO Jet ST (GeV); Number of Events", 50, 0, 1000);
      histContainer_["CrossTrigger_Pass_greatestBtag"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_GreatestBtag", "Pass_GreatestBtag; Greatest RECO Jet pfCSVv2; Number of Events", 100, 0, 1);
      histContainer_["CrossTrigger_Total_greatestBtag"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_GreatestBtag", "Total_GreatestBtag; Greatest RECO Jet pfCSVv2; Number of Events;", 100, 0, 1);     
      histContainer_["CrossTrigger_Pass_forwardJetEta"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_forwardJetEta", "Pass_ForwardJetEta; Forward RECO Jet Eta; Number of Events", 100, 0, 3);
      histContainer_["CrossTrigger_Total_forwardJetEta"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_forwardJetEta", "Total_ForwardJetEta; Forward RECO Jet Eta; Number of Events", 100, 0, 3);
      histContainer_["Total_RECO_JetPt"] = subDir_Observables_Jet.make<TH1F>("Total_RECO_Jet_Pt", "RECO Jet Pt; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);
      histContainer_["Total_RECO_JetBTag"] = subDir_Observables_Jet.make<TH1F>("Total_RECO_Jet_BTag", "RECO Jet BTag; RECO Jet BTag pfCSVv2 -log10(1-CSV); Number of Events", 50, 0, 12);
      histContainer_["Total_RECO_JetCSV"] = subDir_Observables_Jet.make<TH1F>("Total_RECO_Jet_CSV", "RECO Jet BTag; RECO Jet BTag pfCSVv2; Number of Events", 50, 0, 1);

      subDir_Observables_Vertices = subDir_Observables.mkdir( "Vertices" );
      histContainer_["CrossTrigger_Pass_VertexMultHist"] = subDir_Observables_Vertices.make<TH1F>("CrossTrigger_Pass_VertexMultHist", "Pass_Vertex; Vertex Mult; Number of Events", 30, 0, 30);
      histContainer_["CrossTrigger_Total_VertexMultHist"] = subDir_Observables_Vertices.make<TH1F>("CrossTrigger_Total_VertexMultHist", "Total_Vertex; Vertex Mult; Number of Events", 30, 0, 30);

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
            histContainer_["Filter1_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);
            histContainer_["Filter1_matchedJetEta"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Eta", "Filter 1 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter1_matchedJetPhi"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Phi", "Filter 1 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);
            
            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );

            // CSV Filter
            subDir_Filter2 = fileService->mkdir( filter2_.c_str() );

            subDir_Filter2_Observables = subDir_Filter2.mkdir( "Trigger Object Observables" );
            histContainer_["Filter2_Pt"] = subDir_Filter2_Observables.make<TH1F>("Pt", "Pt", 50, 0, 100);
            histContainer_["Filter2_Eta"] = subDir_Filter2_Observables.make<TH1F>("Eta", "Eta", 50, -3, 3);
            histContainer_["Filter2_Phi"] = subDir_Filter2_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter2_MatchedJetObservables = subDir_Filter2.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter2_matchedJetBTag"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet CSV -log10(1-CSV)", "Filter 2 matched Jet CSV; RECO Jet CSV (pfCSVv2) -log10(1-CSV); Number of Events", 50, 0, 12);
            histContainer_["Filter2_matchedJetCSV"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet CSV", "Filter 2 matched Jet CSV; RECO Jet CSV (pfCSVv2); Number of Events", 50, 0, 1);
            histContainer_["Filter2_matchedJetPt"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Pt", "Filter 2 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);
            histContainer_["Filter2_matchedJetEta"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Eta", "Filter 2 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter2_matchedJetPhi"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Phi", "Filter 2 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);
            subDir_Filter2_TurnOnCurves = subDir_Filter2.mkdir( "Turn On Curves" );
      }
      
      if ( hadronicleg_ == "TTBarJet30" ){

            // 30GeV Jet Filter
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Trigger Object Observables" );
            histContainer_["Filter1_Pt"] = subDir_Filter1_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter1_Eta"] = subDir_Filter1_Observables.make<TH1F>("Eta", "Eta", 50, -3, 3);
            histContainer_["Filter1_Phi"] = subDir_Filter1_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter1_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);
            histContainer_["Filter1_matchedJetEta"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Eta", "Filter 1 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter1_matchedJetPhi"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Phi", "Filter 1 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);

            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );
      }
      
      if ( hadronicleg_ == "TTBarJet304050" ){

            // 30GeV Jet Filter
            subDir_Filter1 = fileService->mkdir( filter1_.c_str() );

            // const Int_t xbins = 4;
            // Double_t xbinedges[xbins + 1] = {0, 30, 40, 50, 200};

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Filter Observables" );
            histContainer_["Filter1_Pt"] = subDir_Filter1_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter1_Eta"] = subDir_Filter1_Observables.make<TH1F>("Eta", "Eta", 50, -3, 3);
            histContainer_["Filter1_Phi"] = subDir_Filter1_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter1_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);
            histContainer_["Filter1_matchedJetEta"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Eta", "Filter 1 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter1_matchedJetPhi"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Phi", "Filter 1 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);
            histContainer_["Total_RECO_Filter1Background_JetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Background Filter 1", "Background Filter 1; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);

            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );
            histContainer_["Filter1_dR"] = subDir_Filter1_TurnOnCurves.make<TH1F>("Filter1_dR","Filter1_dR",310, 0 ,0.31);
            distContainer_["Filter1_hlt_vs_Reco_Pt"] = subDir_Filter1_TurnOnCurves.make<TH2F>("Matching", "Filter Object, Matched RECO Jet Pt; Filter Obj Jet Pt (GeV); RECO Jet Pt (GeV)", 50, 0, 200, 50, 0, 200);

            // 40GeV Jet Filter
            subDir_Filter2 = fileService->mkdir( filter2_.c_str() );

            subDir_Filter2_Observables = subDir_Filter2.mkdir( "Filter Observables" );
            histContainer_["Filter2_Pt"] = subDir_Filter2_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter2_Eta"] = subDir_Filter2_Observables.make<TH1F>("Eta", "Eta", 50, -3, 3);
            histContainer_["Filter2_Phi"] = subDir_Filter2_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter2_MatchedJetObservables = subDir_Filter2.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter2_matchedJetPt"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Pt", "Filter 2 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);
            histContainer_["Filter2_matchedJetEta"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Eta", "Filter 2 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter2_matchedJetPhi"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Filter 2 matched Jet Phi", "Filter 2 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);
            histContainer_["Total_RECO_Filter2Background_JetPt"] = subDir_Filter2_MatchedJetObservables.make<TH1F>("Background Filter 2", "Background Filter 2; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);

            subDir_Filter2_TurnOnCurves = subDir_Filter2.mkdir( "Turn On Curves" );
            histContainer_["Filter2_dR"] = subDir_Filter2_TurnOnCurves.make<TH1F>("Filter2_dR","Filter2_dR",310, 0 ,0.31);

            // 50GeV Jet Filter
            subDir_Filter3 = fileService->mkdir( filter3_.c_str() );

            subDir_Filter3_Observables = subDir_Filter3.mkdir( "Filter Observables" );
            histContainer_["Filter3_Pt"] = subDir_Filter3_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter3_Eta"] = subDir_Filter3_Observables.make<TH1F>("Eta", "Eta", 50, -3, 3);
            histContainer_["Filter3_Phi"] = subDir_Filter3_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter3_MatchedJetObservables = subDir_Filter3.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter3_matchedJetPt"] = subDir_Filter3_MatchedJetObservables.make<TH1F>("Filter 3 matched Jet Pt", "Filter 3 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);
            histContainer_["Filter3_matchedJetEta"] = subDir_Filter3_MatchedJetObservables.make<TH1F>("Filter 3 matched Jet Eta", "Filter 3 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter3_matchedJetPhi"] = subDir_Filter3_MatchedJetObservables.make<TH1F>("Filter 3 matched Jet Phi", "Filter 3 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);
            histContainer_["Total_RECO_Filter3Background_JetPt"] = subDir_Filter3_MatchedJetObservables.make<TH1F>("Background Filter 3", "Background Filter 3; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);

            subDir_Filter3_TurnOnCurves = subDir_Filter3.mkdir( "Turn On Curves" );
            histContainer_["Filter3_dR"] = subDir_Filter3_TurnOnCurves.make<TH1F>("Filter3_dR","Filter3_dR",310, 0 ,0.31);
      }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyser::endJob(){
      //reco vs trigger object
     
      histContainer_["CondEff"]->Fill(histContainer_["TypeOfEvent"]->GetBinContent(3)/histContainer_["TypeOfEvent"]->GetBinContent(2));


      // CREATE TRIGGER TURN ON CURVES -------------------------------------------------------------- //

      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_JetMult_20"],*histContainer_["CrossTrigger_Total_JetMult_20"]) ){

            TCanvas *c1 = new TCanvas((crosstrigger_ + "_JetMult_20").c_str(),"c1",600,400);
            c1->SetGrid();
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_20"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMult_20"],histContainer_["CrossTrigger_Total_JetMult_20"]);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_20"]->SetTitle("Efficiency Plot Jet Mult > 20 GeV; Jet Mult; Efficiency");
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_20"]->SetMarkerColor(4);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_20"]->SetMarkerStyle(21);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_20"]->SetMaximum(1);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_20"]->Draw("AP");
            gPad->Update();
            c1->Update();
            c1->Write();
      }

      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_JetMult_30"],*histContainer_["CrossTrigger_Total_JetMult_30"]) ){

            TCanvas *c2 = new TCanvas((crosstrigger_ + "_JetMult_30").c_str(),"c2",600,400);
            c2->SetGrid();
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_30"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMult_30"],histContainer_["CrossTrigger_Total_JetMult_30"]);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_30"]->SetTitle("Efficiency Plot Jet Mult > 30 GeV; Jet Mult; Efficiency");
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_30"]->SetMarkerColor(4);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_30"]->SetMarkerStyle(21);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_30"]->SetMaximum(1);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_30"]->Draw("AP");
            gPad->Update();
            c2->Update();
            c2->Write();
      }

      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_JetMult_40"],*histContainer_["CrossTrigger_Total_JetMult_40"]) ){

            TCanvas *c3 = new TCanvas((crosstrigger_ + "_JetMult_40").c_str(),"c2",600,400);
            c3->SetGrid();
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_40"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMult_40"],histContainer_["CrossTrigger_Total_JetMult_40"]);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_40"]->SetTitle("Efficiency Plot Jet Mult > 40 GeV; Jet Mult; Efficiency");
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_40"]->SetMarkerColor(4);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_40"]->SetMarkerStyle(21);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_40"]->SetMaximum(1);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_40"]->Draw("AP");
            gPad->Update();
            c3->Update();
            c3->Write();
      }

      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_JetMult_50"],*histContainer_["CrossTrigger_Total_JetMult_50"]) ){

            TCanvas *c4 = new TCanvas((crosstrigger_ + "_JetMult_50").c_str(),"c4",600,400);
            c4->SetGrid();
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_50"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMult_50"],histContainer_["CrossTrigger_Total_JetMult_50"]);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_50"]->SetTitle("Efficiency Plot Jet Mult > 50 GeV; Jet Mult; Efficiency");
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_50"]->SetMarkerColor(4);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_50"]->SetMarkerStyle(21);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_50"]->SetMaximum(1);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMult_50"]->Draw("AP");
            gPad->Update();
            c4->Update();
            c4->Write();
      }

      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_HT"],*histContainer_["CrossTrigger_Total_HT"]) ){

            TCanvas *c5 = new TCanvas((crosstrigger_ + "_Global_HT_20").c_str(),"c5",600,400);
            c5->SetGrid();
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Global_HT_20"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_HT"],histContainer_["CrossTrigger_Total_HT"]);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Global_HT_20"]->SetTitle("Efficiency Plot HT (Sum Jets > 20 GeV); Global_HT_20 (GeV); Efficiency");
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Global_HT_20"]->SetMarkerColor(4);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Global_HT_20"]->SetMarkerStyle(21);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Global_HT_20"]->SetMaximum(1);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Global_HT_20"]->Draw("AP");
            gPad->Update();
            c5->Update();
            c5->Write();
      }


      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_METHist"],*histContainer_["CrossTrigger_Total_METHist"]) ){

            TCanvas *c6 = new TCanvas((crosstrigger_ + "_MET").c_str(),"c6",600,400);
            c6->SetGrid();
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_METHist"],histContainer_["CrossTrigger_Total_METHist"]);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->SetTitle("Efficiency Plot MET; MET (GeV); Efficiency");
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->SetMarkerColor(4);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->SetMarkerStyle(21);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->SetMaximum(1);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_MET"]->Draw("AP");
            gPad->Update();
            c6->Update();
            c6->Write();
      }

      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_VertexMultHist"],*histContainer_["CrossTrigger_Total_VertexMultHist"]) ){

            TCanvas *c7 = new TCanvas((crosstrigger_ + "_VertMult").c_str(),"c7",600,400);
            c7->SetGrid();
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_VertexMultHist"],histContainer_["CrossTrigger_Total_VertexMultHist"]);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->SetTitle("Efficiency Plot Vertex Mult; Vertex Mult; Efficiency");
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->SetMarkerColor(4);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->SetMarkerStyle(21);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->SetMaximum(1);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->Draw("AP");
            gPad->Update();
            c7->Update();
            c7->Write();
      }

      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_LeptonPtHist"],*histContainer_["CrossTrigger_Total_LeptonPtHist"]) ){

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
      }


      // CREATE FILTER TURN ON CURVES --------------------------------------------------------------- //

      gStyle->SetOptFit(0111);

      TF1* f1 = new TF1("f1","(0.5*[0]*(1+TMath::Erf((x-[1])/(sqrt(x)*[2]))))",0,100);//pt
      TF1* f2 = new TF1("f2","(0.5*[0]*(1+TMath::Erf((x-[1])/(sqrt(x)*[2]))))",0,10);//btag

      f1->SetParName(0,"Plateau Efficiency");
      f1->SetParName(1,"Pt at Half Plateau Efficiency");
      f1->SetParName(2,"Slope at Half Plateau Efficiency");
      f2->SetParName(0,"Plateau Efficiency");
      f2->SetParName(1,"CSV at Half Plateau Efficiency");
      f2->SetParName(2,"Slope at Half Plateau Efficiency");


      if ( hadronicleg_ == "SingleTop" ){

            if ( TEfficiency::CheckConsistency(*histContainer_["Filter1_matchedJetPt"],*histContainer_["Total_RECO_JetPt"]) ){
                  TCanvas *c9 = new TCanvas((crosstrigger_ + "_FilterTurnOn_Pt30").c_str(),"c9",600,400);
                  c9->SetGrid();
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"] = subDir_Filter1_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter1_matchedJetPt"],histContainer_["Total_RECO_JetPt"]);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetTitle( (filter1_ + " Turn On Jet Pt; RECO Jet Pt (GeV); Efficiency").c_str() );
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerColor(4);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerStyle(21);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMaximum(1.1);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Draw("AP");

                  f1->SetParameters(1.0,30,2.);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Fit(f1);
       
                  gPad->Update();
        
                  // TPaveStats *stats9 = (TPaveStats*)(turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->GetListOfFunctions()->FindObject("stats"));
                  // stats9->SetX1NDC(0.5);
                  // stats9->SetX2NDC(0.9);
                  // stats9->SetY1NDC(0.1);
                  // stats9->SetY2NDC(0.3);

                  c9->Update();
                  c9->Write();
            }

            if ( TEfficiency::CheckConsistency(*histContainer_["Filter2_matchedJetBTag"],*histContainer_["Total_RECO_JetBTag"]) ){

                  TCanvas *c10 = new TCanvas((crosstrigger_ + "_FilterTurnOn_CSV -log10(1-CSV)").c_str(),"c10",600,400);
                  c10->SetGrid();
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"] = subDir_Filter2_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter2_matchedJetBTag"],histContainer_["Total_RECO_JetBTag"]);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetTitle( (filter2_ + " Turn On Jet CSV; RECO Jet CSV; Efficiency").c_str() );
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetMarkerColor(4);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetMarkerStyle(21);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetMaximum(1.1);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->Draw("AP");

                  f2->SetParameters(1.0,2,1.5);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->Fit(f2);

                  gPad->Update();

                  // TPaveStats *stats10 = (TPaveStats*)(turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->GetListOfFunctions()->FindObject("stats"));
                  // stats10->SetX1NDC(0.5);
                  // stats10->SetX2NDC(0.9);
                  // stats10->SetY1NDC(0.1);
                  // stats10->SetY2NDC(0.3);

                  c10->Update();
                  c10->Write();
            }

            if ( TEfficiency::CheckConsistency(*histContainer_["Filter2_matchedJetCSV"],*histContainer_["Total_RECO_JetCSV"]) ){

                  TCanvas *c11 = new TCanvas((crosstrigger_ + "_FilterTurnOn_CSV").c_str(),"c11",600,400);
                  c11->SetGrid();
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"] = subDir_Filter2_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter2_matchedJetCSV"],histContainer_["Total_RECO_JetCSV"]);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetTitle( (filter2_ + " Turn On Jet CSV; RECO Jet CSV; Efficiency").c_str() );
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetMarkerColor(4);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetMarkerStyle(21);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->SetMaximum(1.1);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->Draw("AP");

                  f2->SetParameters(1.0,2,1.5);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->Fit(f2);

                  gPad->Update();

                  // TPaveStats *stats11 = (TPaveStats*)(turnOnCurveContainer_["Filter2_TurnOnCurve_BTag"]->GetListOfFunctions()->FindObject("stats"));
                  // stats11->SetX1NDC(0.5);
                  // stats11->SetX2NDC(0.9);
                  // stats11->SetY1NDC(0.1);
                  // stats11->SetY2NDC(0.3);

                  c11->Update();
                  c11->Write();
            }
      }

      if ( hadronicleg_ == "TTBarJet30" ){

            if ( TEfficiency::CheckConsistency(*histContainer_["Filter1_matchedJetPt"],*histContainer_["Total_RECO_JetPt"]) ){

                  TCanvas *c9 = new TCanvas((crosstrigger_ + "_FilterTurnOn_Pt30").c_str(),"c9",600,400);
                  c9->SetGrid();
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"] = subDir_Filter1_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter1_matchedJetPt"],histContainer_["Total_RECO_JetPt"]);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetTitle( (filter1_ + " Turn On Jet Pt; RECO Jet Pt (GeV); Efficiency").c_str() );
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerColor(4);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerStyle(21);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMaximum(1.1);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Draw("AP");

                  f1->SetParameters(1.0,30,2.);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Fit(f1);

                  gPad->Update();

                  // TPaveStats *stats9 = (TPaveStats*)(turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->GetListOfFunctions()->FindObject("stats"));
                  // stats9->SetX1NDC(0.5);
                  // stats9->SetX2NDC(0.9);
                  // stats9->SetY1NDC(0.1);
                  // stats9->SetY2NDC(0.3);

                  c9->Update();
                  c9->Write();
            }
      }

      if ( hadronicleg_ == "TTBarJet304050" ){

            if ( TEfficiency::CheckConsistency(*histContainer_["Filter1_matchedJetPt"],*histContainer_["Total_RECO_Filter1Background_JetPt"]) ){

                  TCanvas *c9 = new TCanvas((crosstrigger_ + "_FilterTurnOn_Pt30").c_str(),"c9",600,400);
                  c9->SetGrid();
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"] = subDir_Filter1_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter1_matchedJetPt"],histContainer_["Total_RECO_Filter1Background_JetPt"]);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetTitle( (filter1_ + " Turn On Jet Pt; RECO Jet Pt (GeV); Efficiency").c_str() );
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerColor(4);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMarkerStyle(21);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->SetMaximum(1.1);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Draw("AP");

                  f1->SetParameters(1.0,30,2.);
                  turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->Fit(f1);

                  gPad->Update();

                  // TPaveStats *stats9 = (TPaveStats*)(turnOnCurveContainer_["Filter1_TurnOnCurve_Pt"]->GetListOfFunctions()->FindObject("stats"));
                  // stats9->SetX1NDC(0.5);
                  // stats9->SetX2NDC(0.9);
                  // stats9->SetY1NDC(0.1);
                  // stats9->SetY2NDC(0.3);

                  c9->Update();
                  c9->Write();
            }


            if ( TEfficiency::CheckConsistency(*histContainer_["Filter2_matchedJetPt"],*histContainer_["Total_RECO_Filter2Background_JetPt"]) ){

                  TCanvas *c10 = new TCanvas((crosstrigger_ + "_FilterTurnOn_Pt40").c_str(),"c10",600,400);
                  c10->SetGrid();
                  turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"] = subDir_Filter2_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter2_matchedJetPt"],histContainer_["Total_RECO_Filter2Background_JetPt"]);//CHANGED
                  turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->SetTitle( (filter2_ + " Turn On Jet Pt; RECO Jet Pt (GeV); Efficiency").c_str() );
                  turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->SetMarkerColor(4);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->SetMarkerStyle(21);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->SetMaximum(1.1);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->Draw("AP");

                  f1->SetParameters(1.0,40,1.5);
                  turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->Fit(f1);

                  gPad->Update();

                  // TPaveStats *stats10 = (TPaveStats*)(turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"]->GetListOfFunctions()->FindObject("stats"));
                  // stats10->SetX1NDC(0.5);
                  // stats10->SetX2NDC(0.9);
                  // stats10->SetY1NDC(0.1);
                  // stats10->SetY2NDC(0.3);

                  c10->Update();
                  c10->Write();
            }


            if ( TEfficiency::CheckConsistency(*histContainer_["Filter3_matchedJetPt"],*histContainer_["Total_RECO_Filter3Background_JetPt"]) ){

                  TCanvas *c11 = new TCanvas((crosstrigger_ + "_FilterTurnOn_Pt50").c_str(),"c11",600,400);
                  c11->SetGrid();
                  turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"] = subDir_Filter3_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter3_matchedJetPt"],histContainer_["Total_RECO_Filter3Background_JetPt"]);//CHANGED
                  turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->SetTitle( (filter3_ + " Turn On Jet Pt; RECO Jet Pt (GeV); Efficiency").c_str() );
                  turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->SetMarkerColor(4);
                  turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->SetMarkerStyle(21);
                  turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->SetMaximum(1.1);
                  turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->Draw("AP");

                  f1->SetParameters(1.0,50,1.5);
                  turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->Fit(f1);
                  gPad->Update();

                  // TPaveStats *stats11 = (TPaveStats*)(turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"]->GetListOfFunctions()->FindObject("stats"));
                  // stats11->SetX1NDC(0.5);
                  // stats11->SetX2NDC(0.9);
                  // stats11->SetY1NDC(0.1);
                  // stats11->SetY2NDC(0.3);
                  
                  c11->Update();
                  c11->Write();
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


