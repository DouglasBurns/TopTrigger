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


      // EVENT SELECTION AND INDEX STORING ---------------------------------------------------------- //


      // Perform Event Selection
      passMockEventSelection = false;
      leadingElectronIndex.clear();
      leadingMuonIndex.clear();
      cleanedJetIndex.clear();
      cleanedBJetIndex.clear();

      // Tight Electrons
      Index = 0;
      if (electrons->size() != 0){
            for  (auto lepton = electrons->begin(); lepton != electrons->end(); ++lepton){
                  if (leptontype_ == "Mu20") ptcut=30;
                  if (leptontype_ == "Mu24") ptcut=30;
                  if (leptontype_ == "Ele27") ptcut=30;
                  if (leptontype_ == "Ele32") ptcut=35;

                  if (lepton->pt() < ptcut && abs(lepton->eta()) ) continue;
                  if (lepton->electronID(electronID_)){
                        leadingElectronIndex.push_back(Index);
                  }
                  ++Index;
            }
      }

      // Tight Isolated Muons
      Index = 0;
      // std::cout << "Number of Vertices : " << vertices->size() << std::endl;
      if (vertices->size() != 0){
            auto vertex = vertices->at(0);
            if (vertex.isValid()) {
                   // std::cout << "Number of Muons : " << muons->size() << std::endl;

                  if (muons->size() != 0){
                        for (auto lepton = muons->begin(); lepton != muons->end(); ++lepton){
                              isolated_muon = true;
                              if (leptontype_ == "Mu20") ptcut=22;
                              if (leptontype_ == "Mu24") ptcut=27;
                              if (leptontype_ == "Ele27") ptcut=22;
                              if (leptontype_ == "Ele32") ptcut=22;

                              if (lepton->pt() > ptcut && abs(lepton->eta()) < 2.1){

                                    if (lepton->isTightMuon(vertex)){
                                          isolation = ( (lepton->chargedHadronIso() + std::max( lepton->neutralHadronIso() + lepton->photonIso() - 0.5 * lepton->puChargedHadronIso(), 0.) ) / lepton->pt() );
                                          // std::cout << isolation << std::endl;
                                          
                                          // // std::cout << "Muon is Tight" << std::endl;
                                          // for( auto jet = jets->begin(); jet != jets->end(); ++jet ){

                                          //       if (!isGoodJet(*jet)) continue;

                                          //       if (reco::deltaR2(*lepton, *jet) < 0.09) isolated_muon = false;
                                                
                                          // }
                                          // if (isolated_muon) leadingMuonIndex.push_back(Index);
                                          if (isolation < 0.12) leadingMuonIndex.push_back(Index);
                                    }
                              }
                              ++Index;
                        }
                  }
            }
      }

      // std::cout << "Size : " << leadingElectronIndex.size() << "  Indices : ";
      // for (uint x = 0; x<leadingElectronIndex.size(); ++x) std::cout << " " << leadingElectronIndex[x];
      // std::cout << std::endl;

      // std::cout << "Muon Index Size : " << leadingMuonIndex.size() << "  Indices : ";
      // for (uint x = 0; x<leadingMuonIndex.size(); ++x) std::cout << " " << leadingMuonIndex[x];
      // std::cout << std::endl;

      Index = 0;
      if (jets->size() != 0){
            for( auto jet = jets->begin(); jet != jets->end(); ++jet ){
                  isMatchedToLepton = false;
                  if (isGoodJet(*jet)){
      
                        if ( leadingElectronIndex.size() != 0 ){
                              for (uint x = 0; x<leadingElectronIndex.size(); ++x){
                                    if ( leadingElectronIndex.size() <= electrons->size() ) continue;
                                    auto lepton = electrons->at(leadingElectronIndex[x]);                              
                                    double const dR2 = reco::deltaR2(lepton, *jet);
                                    if (dR2 < 0.3 * 0.3) isMatchedToLepton = true;
                                    // if (dR2 < 0.3 * 0.3) std::cout << "Woooo! Matched with deltaR : " << dR2 << ", with ele pt of : " << lepton.pt() << " to reco jet pt of : " << jet->pt() << std::endl;

                              }
                        }
                  
                        if ( leadingMuonIndex.size() != 0 ){
                              for (uint x = 0; x<leadingMuonIndex.size(); ++x){
                                    if ( leadingMuonIndex.size() <= muons->size() ) continue;
                                    auto lepton = muons->at(leadingMuonIndex[x]);                              
                                    double const dR2 = reco::deltaR2(lepton, *jet);
                                    if (dR2 < 0.3 * 0.3) isMatchedToLepton = true;
                              }
                        }

                        if (!isMatchedToLepton) {
                              cleanedJetIndex.push_back(Index);
                              if (jet->bDiscriminator(btagger_) > 0.890 ){// 2 med bjet value (https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X50ns#Supported_Algorithms_and_Operati)
                                    cleanedBJetIndex.push_back(Index);
                              }
                        }
                  }
                  ++Index;
            }
      }

      // std::cout << "Jet Index Size : " << cleanedJetIndex.size() << "  Indices : ";
      // for (uint x = 0; x<cleanedJetIndex.size(); ++x) std::cout << " " << cleanedJetIndex[x];
      // std::cout << " Pt : ";
      // for (uint x = 0; x<cleanedJetIndex.size(); ++x){
      // auto jet = jets->at(cleanedJetIndex[x]);
      // std::cout << " " << jet.pt();
      // } 
      // std::cout << std::endl;

      // std::cout << "Size : " << cleanedBJetIndex.size() << "  Indices : ";
      // for (uint x = 0; x<cleanedBJetIndex.size(); ++x) std::cout << " " << cleanedBJetIndex[x];
      // std::cout << std::endl;


      // Perform an offline event selection
      if (leptonicleg_ == "Ele"){
            if (hadronicleg_ == "SingleTop"){
                  if (cleanedJetIndex.size() >= 2 && cleanedBJetIndex.size() >= 1 && leadingMuonIndex.size() == 0 && leadingElectronIndex.size() == 1) passMockEventSelection = true;
            }
             if (hadronicleg_ == "TTBarJet30" || hadronicleg_ == "TTBarJet304050"){
                  if (cleanedJetIndex.size() >= 4 && cleanedBJetIndex.size() >= 2 && leadingMuonIndex.size() == 0 && leadingElectronIndex.size() == 1) passMockEventSelection = true;
            }
      }
      if (leptonicleg_ == "Mu"){
            if (hadronicleg_ == "SingleTop"){
                  if (cleanedJetIndex.size() >= 2 && cleanedBJetIndex.size() >= 1 && leadingMuonIndex.size() == 1 && leadingElectronIndex.size() == 0) passMockEventSelection = true;
            }
             if (hadronicleg_ == "TTBarJet30" || hadronicleg_ == "TTBarJet304050"){
                  if (cleanedJetIndex.size() >= 4 && cleanedBJetIndex.size() >= 2 && leadingMuonIndex.size() == 1 && leadingElectronIndex.size() == 0) passMockEventSelection = true;
            }
      }


      // std::cout << "Pass Event Selection : " << passMockEventSelection << std::endl;

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

      SingleLeptonTrigDecision = false;
      CrossTriggerTrigDecision = false;

      TypeOfEvent=0;
      histContainer_["TypeOfEvent"]->Fill(TypeOfEvent);

      if (passMockEventSelection){
            TypeOfEvent=1;
            histContainer_["TypeOfEvent"]->Fill(TypeOfEvent);

            if ( singleleptonIndex < triggerResults->size() ) {
                  if(triggerResults->accept(singleleptonIndex)){
                        SingleLeptonTrigDecision = true;
                        TypeOfEvent=2;
                        histContainer_["TypeOfEvent"]->Fill(TypeOfEvent);
                  }
            }
            else std::cout << "Exception : Looking for " << singleleptontrigger_ << " but failed" << std::endl;
            
            if ( crossIndex < triggerResults->size() ) {
                  if(triggerResults->accept(crossIndex)){
                        CrossTriggerTrigDecision = true;
                        TypeOfEvent=3;
                        histContainer_["TypeOfEvent"]->Fill(TypeOfEvent);
                  }
            }
            else std::cout << "Exception : Looking for " << crosstrigger_ << " but failed" << std::endl;


            // DIFFERENTIAL EFFICIENCIES ------------------------------------------------------------------ //
      

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

            for (uint Index = 0; Index<cleanedJetIndex.size(); ++Index){
                  auto jet = jets->at(cleanedJetIndex[Index]);   

                  // Number of jets in an event and HT depending on jet pt cut
                  if ( jet.pt()>20){
                        ++jetMultiplicity_20;
                  } 
                  if ( jet.pt()>30){
                        ++jetMultiplicity_30;
                  }
                  if ( jet.pt()>40){
                        ++jetMultiplicity_40;
                  } 
                  if ( jet.pt()>50){
                        ++jetMultiplicity_50;
                  } 

                  // HT is the sum of cleaned central jets over 20GeV (No MET or Leptons) for events that pass selection
                  HT += jet.pt();

                  // Most forward jet has highest eta
                  if (std::abs(jet.eta()) >= forwardjeteta) forwardjeteta = std::abs(jet.eta());

                  // Jet with greatest btag value
                  if (jet.bDiscriminator(btagger_) >= jetCSV) jetCSV = jet.bDiscriminator(btagger_);

                  // Store Jet Pt, Eta, Higest CSV, Global_HT and Multiplicity in histograms 
                  if (CrossTriggerTrigDecision==true){
                        histContainer_["CrossTrigger_Pass_JetPtHist"]->Fill(jet.pt());
                        histContainer_["CrossTrigger_Pass_JetEtaHist"]->Fill(jet.eta());
                        histContainer_["CrossTrigger_Pass_JetPhiHist"]->Fill(jet.phi());
                        histContainer_["CrossTrigger_Pass_JetCSVHist"]->Fill(jet.bDiscriminator(btagger_));
                  }
                  histContainer_["CrossTrigger_Total_JetPtHist"]->Fill(jet.pt());
                  histContainer_["CrossTrigger_Total_JetEtaHist"]->Fill(jet.eta());
                  histContainer_["CrossTrigger_Total_JetPhiHist"]->Fill(jet.phi());
                  histContainer_["CrossTrigger_Total_JetCSVHist"]->Fill(jet.bDiscriminator(btagger_));
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

                  metEnergy = met->energy();

                  if(CrossTriggerTrigDecision==true){
                        histContainer_["CrossTrigger_Pass_METHist"]->Fill(metEnergy);
                  }
            histContainer_["CrossTrigger_Total_METHist"]->Fill(metEnergy);
            }


            // LEADING LEPTONS ---------------------------------------------------------------------------- //
     
            // Store the leading lepton pt, eta and energy distributions in histograms
            if (leptonicleg_ == "Ele"){            
                  for (uint Index = 0; Index<leadingElectronIndex.size(); ++Index){
                  auto lepton = electrons->at(leadingElectronIndex[Index]);  

                        if(CrossTriggerTrigDecision==true){
                              histContainer_["CrossTrigger_Pass_LeptonPtHist"]->Fill(lepton.pt());
                              histContainer_["CrossTrigger_Pass_LeptonEtaHist"]->Fill(lepton.eta());
                              histContainer_["CrossTrigger_Pass_LeptonPhiHist"]->Fill(lepton.phi());
                              histContainer_["CrossTrigger_Pass_LeptonEnergyHist"]->Fill(lepton.energy());
                        }
                        histContainer_["CrossTrigger_Total_LeptonPtHist"]->Fill(lepton.pt());
                        histContainer_["CrossTrigger_Total_LeptonEtaHist"]->Fill(lepton.eta());
                        histContainer_["CrossTrigger_Total_LeptonPhiHist"]->Fill(lepton.phi());
                        histContainer_["CrossTrigger_Total_LeptonEnergyHist"]->Fill(lepton.energy());
                  }
                  if (leadingElectronIndex.size() > 1) std::cout << "Warning : More than 1 lepton has passed event selection" << std::endl;
            }

            if ( leptonicleg_ == "Mu" ){
                  for (uint Index = 0; Index<leadingMuonIndex.size(); ++Index){
                  auto lepton = muons->at(leadingMuonIndex[Index]); 

                        if(CrossTriggerTrigDecision==true){
                              histContainer_["CrossTrigger_Pass_LeptonPtHist"]->Fill(lepton.pt());
                              histContainer_["CrossTrigger_Pass_LeptonEtaHist"]->Fill(lepton.eta());
                              histContainer_["CrossTrigger_Pass_LeptonPhiHist"]->Fill(lepton.phi());
                              histContainer_["CrossTrigger_Pass_LeptonEnergyHist"]->Fill(lepton.energy());
                        }
                        histContainer_["CrossTrigger_Total_LeptonPtHist"]->Fill(lepton.pt());
                        histContainer_["CrossTrigger_Total_LeptonEtaHist"]->Fill(lepton.eta());
                        histContainer_["CrossTrigger_Total_LeptonPhiHist"]->Fill(lepton.phi());
                        histContainer_["CrossTrigger_Total_LeptonEnergyHist"]->Fill(lepton.energy());
                  }
                  if (leadingMuonIndex.size() > 1) std::cout << "Warning : More than 1 lepton has passed event selection" << std::endl;      
            } 
      }

      // FILTERS: TRIGGER OBJECTS ------------------------------------------------------------------- //
      // std::cout << "++++++++++++++++++" << std::endl;
      //If it doesnt pass the mock event selection then it is still considered in filter.
      if ( singleleptonIndex < triggerResults->size() ) {
            if(triggerResults->accept(singleleptonIndex)) {
                  SingleLeptonTrigDecision = true;
            }
      }

      if (SingleLeptonTrigDecision){

            cleanedRECOJetIndex.clear();
            isLeadingLepton = false;
            std::vector<float> leadingLeptonPt;
            std::vector<float> leadingLeptonEta;
            leadingLeptonPt.clear();
            leadingLeptonEta.clear();

            for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
                  if (obj.collection() != "hltL3MuonCandidates::HLT" && obj.collection() != "hltEgammaCandidates::HLT") continue;
                  if ( obj.hasFilterLabel(leptonfilter_)){
                        leadingLeptonPt.push_back(obj.pt());
                        leadingLeptonEta.push_back(obj.eta());
                  } 
            }



            Index = 0;
            // bool isMatch;
            for( auto jet = jets->begin(); jet != jets->end(); ++jet ){
                  if (isGoodJet(*jet)){
                        // Looking for HLT Leptons that pass the lepton filter              
                        bool isMatched=false;
                        bool isMatchedToHLTJet=false;

                        for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
                              if ( obj.collection() == "hltL3MuonCandidates::HLT" || obj.collection() == "hltEgammaCandidates::HLT"){
                                    if ( obj.hasFilterLabel(leptonfilter_)) {
                                          if (reco::deltaR2(obj, *jet) < (0.3 * 0.3)) isMatched=true;
                                    }
                              }

                              else if ( obj.collection() == "hltPFJetForBtag::HLT" || obj.collection() == "hltAK4PFJetsCorrected::HLT"){
                                    if (reco::deltaR2(obj, *jet) < (0.3 * 0.3)) isMatchedToHLTJet=true;
                              } //Only include RECO jets that have a matching to an hlt jet filter object.
                        }
                        if (!isMatched && isMatchedToHLTJet) cleanedRECOJetIndex.push_back(Index);
                        // else{
                        //       std::cout << "Is matched to Lepton : " << isMatched << std::endl 
                        //       << "Is matched to Jet object : " << isMatchedToHLTJet << std::endl;
                        // }
                  }
                  ++Index;
            }

            // // Looking for HLT Leptons that pass the lepton filter              
            // for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 

            //       // Only clean vs leading lepton
            //       // if (isLeadingLepton) break;
            //       if (obj.collection() != "hltL3MuonCandidates::HLT" && obj.collection() != "hltEgammaCandidates::HLT") continue;

            //       // Finds a lepton passing filter
            //       if ( obj.hasFilterLabel(leptonfilter_)) {
            //             // std::cout << "Found the leading lepton in the Collection : " << obj.collection() << std::endl;
            //             // std::cout << "Leading lepton filter object Pt : " << obj.pt() << ", Eta : " << obj.eta() << std::endl;

            //             // Once found stop looking for more leptons
            //             // isLeadingLepton=true;
            //             Index = 0;
            //             for( auto jet = jets->begin(); jet != jets->end(); ++jet ){
            //                   if (isGoodJet(*jet)){

            //                         // Match leading lepton to current jet and if not matched store jet
            //                         if (reco::deltaR2(obj, *jet) > (0.3 * 0.3)){
            //                               // if a jet is not matched (i.e. clean) then add to list of cleaned jets
            //                               cleanedRECOJetIndex.push_back(Index);
            //                         }
            //                         // else{
            //                         // std::cout << "For : " << leptontype_ << ", Jet of Index : " << Index <<  " is matched with deltaR : " << reco::deltaR2(obj, *jet) << " to obj pt of : " << obj.pt() << ", with reco pt of : " << jet->pt() << std::endl;
            //                         // }
            //                   }
            //                   // else std::cout << "Jet of Index : " << Index << " is not good, Pt : " << jet->pt() << ", Eta : " << jet->eta() << std::endl;
            //                   ++Index;                                    
            //             }
            //       }
            // }



            // std::cout << "Size : " << cleanedRECOJetIndex.size() << "  Indices : ";
            // for (uint x = 0; x<cleanedRECOJetIndex.size(); ++x) std::cout << " " << cleanedRECOJetIndex[x];
            // std::cout << std::endl;




            if (leadingLeptonPt.size() == 0){
                  std::cout << "Print Out : No leading lepton found in this event." << std::endl;
            }

            // if (cleanedRECOJetIndex.size() < 1){
            //       std::cout << "Print Out : There are no cleaned jets left" << std::endl;
            // }

            if (cleanedRECOJetIndex.size() != 0){

                  for (uint Index = 0; Index<cleanedRECOJetIndex.size(); ++Index){
                        auto jet = jets->at(cleanedRECOJetIndex[Index]); 

                        // if (jet.pt() > 30 && jet.pt() < 40) histContainer_["NoJets"]->Fill(0);
                        // if (jet.pt() > 40 && jet.pt() < 50) histContainer_["NoJets"]->Fill(1);
                        // if (jet.pt() > 50) histContainer_["NoJets"]->Fill(2);

                        // Fill clean central jets into histograms
                        histContainer_["Total_RECO_JetPt"]->Fill(jet.pt());
                        histContainer_["Total_RECO_JetBTag"]->Fill(-log10(1-jet.bDiscriminator(btagger_)));
                        histContainer_["Total_RECO_JetCSV"]->Fill(jet.bDiscriminator(btagger_));

                        passesFilter1 = false, passesFilter2 = false, passesFilter3 = false;
                        float minF1dR2 = 9999, minF2dR2 = 9999, minF3dR2 = 9999;
                        bool pF1 = false, pF2 = false, pF3 = false;
                        // Find matching filter object,
                        for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 

                              isJetCollection = false;
                              isBJetCollection = false;

                              if (obj.collection() == "hltPFJetForBtag::HLT") isBJetCollection = true;
                              if (obj.collection() == "hltAK4PFJetsCorrected::HLT") isJetCollection = true;


                              if (!isJetCollection && !isBJetCollection && (obj.hasFilterLabel(filter1_) 
                                    || obj.hasFilterLabel(filter2_) || obj.hasFilterLabel(filter3_) ) ) 
                                    std::cout << "Print Out : Filter is in incorrect collection " << obj.collection() << std::endl;

                              if (!isJetCollection && !isBJetCollection) continue;
                              // By looking through the filter labels (which filters are associated with filter object)
 
                              // Make into function?????
                              if ( hadronicleg_ == "SingleTop" ){

                                    // Keeping the filter objects that match the filters we are interested in
                                    if ( obj.hasFilterLabel(filter1_) ) {

                                          // Finally the matching takes place and histograms are filled
                                          double const dR2 = reco::deltaR2(obj, jet);
                                          histContainer_["DeltaR2 Matching"]->Fill(dR2);

                                          if (dR2 < 0.3 * 0.3){
                                                // std::cout << "Woooo! Matched with deltaR : " << dR2 << ", with obj pt of : " << 
                                                // obj.pt() << " to reco pt of : " << jet.pt() << std::endl;
                                                
                                                // Store RECO jet only once if it is matched 
                                                if (!passesFilter1){
                                                      histContainer_["Filter1_Pt"]->Fill(obj.pt());
                                                      histContainer_["Filter1_Eta"]->Fill(obj.eta());
                                                      histContainer_["Filter1_Phi"]->Fill(obj.phi());

                                                      histContainer_["Filter1_matchedJetPt"]->Fill(jet.pt());
                                                      histContainer_["Filter1_matchedJetEta"]->Fill(jet.eta());
                                                      histContainer_["Filter1_matchedJetPhi"]->Fill(jet.phi());

                                                      distContainer_["Filter1_TriggerObject_Reco_Pt"]->Fill(obj.pt(),jet.pt());
                                                                                                            
                                                }
                                                passesFilter1 = true;
                                          }
                                    }
                              
                                    if ( obj.hasFilterLabel(filter2_) ) {
                                          double const dR2 = reco::deltaR2(obj, jet);
                                          if (dR2 < 0.3 * 0.3){
                                                if (!passesFilter2){
                                                      histContainer_["Filter2_Pt"]->Fill(obj.pt());
                                                      histContainer_["Filter2_Eta"]->Fill(obj.eta());
                                                      histContainer_["Filter2_Phi"]->Fill(obj.phi());

                                                      histContainer_["Filter2_matchedJetPt"]->Fill(jet.pt());
                                                      histContainer_["Filter2_matchedJetEta"]->Fill(jet.eta());
                                                      histContainer_["Filter2_matchedJetPhi"]->Fill(jet.phi());
                                                      histContainer_["Filter2_matchedJetBTag"]->Fill(-log10(1-jet.bDiscriminator(btagger_)));
                                                      histContainer_["Filter2_matchedJetCSV"]->Fill(jet.bDiscriminator(btagger_));
                                                }
                                                passesFilter2 = true;
                                          }
                                    }

                              }
                  

                              if ( hadronicleg_ == "TTBarJet30" ){

                                    if ( obj.hasFilterLabel(filter1_) ) {

                                          double const dR2 = reco::deltaR2(obj, jet);
                                          if (dR2 < 0.3 * 0.3){
                                                if (!passesFilter1){
                                                      histContainer_["Filter1_Pt"]->Fill(obj.pt());
                                                      histContainer_["Filter1_Eta"]->Fill(obj.eta());
                                                      histContainer_["Filter1_Phi"]->Fill(obj.phi());

                                                      histContainer_["Filter1_matchedJetPt"]->Fill(jet.pt());
                                                      histContainer_["Filter1_matchedJetEta"]->Fill(jet.eta());
                                                      histContainer_["Filter1_matchedJetPhi"]->Fill(jet.phi());
                                                }
                                                passesFilter1 = true;
                                          }
                                    }
                              }
                              

                              if (reco::deltaR2(obj, jet) < 0.09 ){
                                    if (obj.pt() > 30 ) histContainer_["NoJets"]->Fill(0);
                                    if (obj.pt() > 40 ) histContainer_["NoJets"]->Fill(2);
                                    if (obj.pt() > 50 ) histContainer_["NoJets"]->Fill(4);
                              }
                              if ( hadronicleg_ == "TTBarJet304050" ){

                                    if ( obj.hasFilterLabel(filter1_)) {

                                          double const dR2 = reco::deltaR2(obj, jet);

                                          if (jet.pt() > 100 && abs(jet.eta()<2.0) ){
                                                if (dR2 < minF1dR2){
                                                      minF1dR2 = dR2;
                                                      pF1 = true;
                                                }
                                          }

                                          if (dR2 < 0.3 * 0.3){

                                                if (!passesFilter1){
                                                      histContainer_["NoJets"]->Fill(1);
                                                      histContainer_["Filter1_Pt"]->Fill(obj.pt());
                                                      histContainer_["Filter1_Eta"]->Fill(obj.eta());
                                                      histContainer_["Filter1_Phi"]->Fill(obj.phi());

                                                      histContainer_["Filter1_matchedJetPt"]->Fill(jet.pt());
                                                      histContainer_["Filter1_matchedJetEta"]->Fill(jet.eta());
                                                      histContainer_["Filter1_matchedJetPhi"]->Fill(jet.phi());
                                                }
                                                passesFilter1 = true;
                                         }
                                    }

                                    if (  obj.hasFilterLabel(filter2_) ) {

                                          double const dR2 = reco::deltaR2(obj, jet);

                                          if (jet.pt() > 100 && abs(jet.eta()<2.0) ){

                                                if (dR2 < minF2dR2){
                                                      minF2dR2 = dR2;
                                                      pF2 = true;
                                                }
                                          }

                                          if (dR2 < 0.3 * 0.3){
                                                if (!passesFilter2){
                                                      histContainer_["NoJets"]->Fill(3);
                                                      histContainer_["Filter2_Pt"]->Fill(obj.pt());
                                                      histContainer_["Filter2_Eta"]->Fill(obj.eta());
                                                      histContainer_["Filter2_Phi"]->Fill(obj.phi());

                                                      histContainer_["Filter2_matchedJetPt"]->Fill(jet.pt());
                                                      histContainer_["Filter2_matchedJetEta"]->Fill(jet.eta());
                                                      histContainer_["Filter2_matchedJetPhi"]->Fill(jet.phi());
                                                }
                                                passesFilter2 = true;
                                         }
                                    }

                                    if ( obj.hasFilterLabel(filter3_) ) {

                                          double const dR2 = reco::deltaR2(obj, jet);

                                          if (jet.pt() > 100 && abs(jet.eta()<2.0) ){

                                                if (dR2 < minF3dR2){
                                                      minF3dR2 = dR2;
                                                      pF3 = true;
                                                }
                                          }

                                          if (dR2 < 0.3 * 0.3){
                                                if (!passesFilter3){
                                                      histContainer_["NoJets"]->Fill(5);
                                                      histContainer_["Filter3_Pt"]->Fill(obj.pt());
                                                      histContainer_["Filter3_Eta"]->Fill(obj.eta());
                                                      histContainer_["Filter3_Phi"]->Fill(obj.phi());

                                                      histContainer_["Filter3_matchedJetPt"]->Fill(jet.pt());
                                                      histContainer_["Filter3_matchedJetEta"]->Fill(jet.eta());
                                                      histContainer_["Filter3_matchedJetPhi"]->Fill(jet.phi());
                                                }
                                                passesFilter3 = true;
                                         }
                                    }
                              }       
                        }

                  
                        
                        if (pF1) histContainer_["Filter1_dR"]->Fill(minF1dR2);
                        if (pF2) histContainer_["Filter2_dR"]->Fill(minF2dR2);
                        if (pF3) histContainer_["Filter3_dR"]->Fill(minF3dR2);
                  }     
            }




            // if  (hadronicleg_ == "TTBarJet304050" ){

            //       std::cout << std::endl << "________________________________________" << std::endl;
            //       std::cout << "_______________NEW EVENT________________" << std::endl;

            //       std::cout << "Event has : " << leadingLeptonPt.size() << " leptons passing the filter" << std::endl;
            //       for( uint i=0;i<leadingLeptonPt.size();i++){
            //             std::cout << "Lepton Pt : " << leadingLeptonPt[i] << ", Lepton Eta : " << leadingLeptonEta[i] << std::endl;
            //       }

            //       std::cout << "Size of Cleaned RECO Jets against all leptons: " << cleanedRECOJetIndex.size() << "  Indices : ";
            //       for (uint x = 0; x<cleanedRECOJetIndex.size(); ++x) std::cout << " " << cleanedRECOJetIndex[x];
            //       std::cout << std::endl;

            //       int count1 = 0, count2 = 0, count3 = 0, count4 = 0, count5 = 0, count6 = 0;

            //       for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 

            //             isJetCollection = false;
            //             if (obj.collection() == "hltAK4PFJetsCorrected::HLT") isJetCollection = true;

            //             // Only list for first run through of jets.
            //             if (isJetCollection && std::abs(obj.eta()) < 2.6 && obj.pt() > 20){

            //                   std::cout << "__________NEW HLT FILTER OBJ____________" << std::endl;

            //                   std::cout << "Obj Pt is : " << obj.pt() << ", Obj Eta is : " << obj.eta() << std::endl;
            //                   std::cout << "Obj Collection is : " << obj.collection() << std::endl;
            //                   std::cout << "Obj Filter Ids are : ";
            //                   for (unsigned int i = 0; i < obj.filterIds().size(); ++i) std::cout << " " << obj.filterIds()[i];
            //                   std::cout << std::endl;

            //                   if (obj.pt() > 30) ++count1;
            //                   if (obj.pt() > 40) ++count2;
            //                   if (obj.pt() > 50) ++count3;
                                    
            //                   if ( obj.hasFilterLabel(filter1_) ) std::cout << "Found filter : " << filter1_ << "filter" << std::endl, ++count4;
            //                   if ( obj.hasFilterLabel(filter2_) ) std::cout << "Found filter : " << filter2_ << "filter" << std::endl, ++count5;
            //                   if ( obj.hasFilterLabel(filter3_) ) std::cout << "Found filter : " << filter3_ << "filter" << std::endl, ++count6;
                                        
            //             }
            //       }

            //       std::cout << "________________________________________" << std::endl;
            //       std::cout << "Number of hlt obj jets above 30 GeV : " << count1 << std::endl;
            //       std::cout << "Number of 30GeV filter instances found : " << count4 << std::endl;
            //       std::cout << "Number of hlt obj jets above 40 GeV : " << count2 << std::endl;
            //       std::cout << "Number of 40GeV filter instances found : " << count5 << std::endl;
            //       std::cout << "Number of hlt obj jets above 50 GeV : " << count3 << std::endl;
            //       std::cout << "Number of 50GeV filter instances found : " << count6 << std::endl;
            // }
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


      // if (!passesPtAndEta) std::cout << "Does not pass jet pt/eta requirements" << std::endl;
      // if (!passesJetID){
      //       std::cout << "Does not pass jet ID because : ";
      //       if (!passNHF) std::cout << "Has not passed Neutral Hadron Energy Fraction < 0.99 : " << jet.neutralHadronEnergyFraction() << std::endl;
      //       if (!passNEMF) std::cout << "Has not passed Neutral EM Energy Fraction < 0.99 : " << jet.neutralEmEnergyFraction() << std::endl;
      //       if (!passMUF) std::cout << "Has not passed muon Energy Fraction < 0.8 : " << jet.muonEnergyFraction() << std::endl;
      //       if (!passNumConst) std::cout << "Has not passed Number of Constituents > 1: " << jet.chargedMultiplicity()+jet.neutralMultiplicity() << std::endl;
      //       if (!passCHF) std::cout << "Has not passed Charged Hadron Energy Fraction > 0. : " << jet.chargedHadronEnergyFraction() << std::endl;
      //       if (!passCHM) std::cout << "Has not passed Charged Multiplicity > 0 : " << jet.chargedMultiplicity() << std::endl;
      //       if (!passCEMF) std::cout << "Has not passed Charged EM Energy Fraction < 0.99 : " << jet.chargedEmEnergyFraction() << std::endl;
      // }
      return passesPtAndEta && passesJetID;
     
}

// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyser::beginJob(){





      // INITIALISE DIRECTORIES AND HISTOGRAMS ------------------------------------------------------ //

      subDir_TrigDec = fileService->mkdir( "Trigger Decision" );

      histContainer_["TypeOfEvent"] = subDir_TrigDec.make<TH1F>("TypeOfEvent", "Type of Event; Total/Selection/SingleLepton/Cross; Number of Events", 4, -0.5, 3.5);
      // histContainer_["TypeOfEvent"] = subDir_TrigDec.make<TH1F>("TypeOfEvent", "Type of Event; Total/Selection/SingleLepton/Cross; Number of Events", 6, -2.5, 3.5);
      histContainer_["CondEff"] = subDir_TrigDec.make<TH1F>("ConditionalEff", "Conditional Efficiency; Efficiency; ", 100, 0, 1);
      histContainer_["NoJets"] = subDir_TrigDec.make<TH1F>("NoJets", "NoJets; Jet Pt Category; ", 6, -0.5, 5.5);
      histContainer_["NoJets"]->GetXaxis()->SetBinLabel(1, "Jets >30");
      histContainer_["NoJets"]->GetXaxis()->SetBinLabel(2, "Jets >30+PF");
      histContainer_["NoJets"]->GetXaxis()->SetBinLabel(3, "Jets >40");
      histContainer_["NoJets"]->GetXaxis()->SetBinLabel(4, "Jets >40+PF");
      histContainer_["NoJets"]->GetXaxis()->SetBinLabel(5, "Jets >50");
      histContainer_["NoJets"]->GetXaxis()->SetBinLabel(6, "Jets >50+PF");

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
      histContainer_["CrossTrigger_Total_greatestBtag"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_GreatestBtag", "Total_GreatestBtag; Greatest RECO Jet pfCSVv2; Number of Events;", 100, 0, 1);     
      histContainer_["CrossTrigger_Pass_forwardJetEta"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Pass_forwardJetEta", "Pass_ForwardJetEta; Forward RECO Jet Eta; Number of Events", 100, 0, 3);
      histContainer_["CrossTrigger_Total_forwardJetEta"] = subDir_Observables_Jet.make<TH1F>("CrossTrigger_Total_forwardJetEta", "Total_ForwardJetEta; Forward RECO Jet Eta; Number of Events", 100, 0, 3);
      histContainer_["Total_RECO_JetPt"] = subDir_Observables_Jet.make<TH1F>("Total_RECO_Jet_Pt", "RECO Jet Pt; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);
      histContainer_["Total_RECO_JetBTag"] = subDir_Observables_Jet.make<TH1F>("Total_RECO_Jet_BTag", "RECO Jet BTag; RECO Jet BTag pfCSVv2 -log10(1-CSV); Number of Events", 50, 0, 12);
      histContainer_["Total_RECO_JetCSV"] = subDir_Observables_Jet.make<TH1F>("Total_RECO_Jet_CSV", "RECO Jet BTag; RECO Jet BTag pfCSVv2; Number of Events", 50, 0, 1);

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
            histContainer_["Filter1_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);
            histContainer_["Filter1_matchedJetEta"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Eta", "Filter 1 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter1_matchedJetPhi"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Phi", "Filter 1 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);
            
            distContainer_["Filter1_TriggerObject_Reco_Pt"] = subDir_Filter1_MatchedJetObservables.make<TH2F>("Matching", "Filter Object, Matched RECO Jet Pt; Filter Obj Jet Pt (GeV); RECO Jet Pt (GeV)", 50, 0, 200, 50, 0, 200);//cndjschdjkshcjkdnsk
            histContainer_["DeltaR2 Matching"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("DeltaR","DeltaR2 Distribution; DeltaR2; Number of Events", 150, 0, 4);

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

            subDir_Filter1_Observables = subDir_Filter1.mkdir( "Filter Observables" );
            histContainer_["Filter1_Pt"] = subDir_Filter1_Observables.make<TH1F>("Pt", "Pt", 50, 0, 200);
            histContainer_["Filter1_Eta"] = subDir_Filter1_Observables.make<TH1F>("Eta", "Eta", 50, -3, 3);
            histContainer_["Filter1_Phi"] = subDir_Filter1_Observables.make<TH1F>("Phi", "Phi", 50, -3.5, 3.5);

            subDir_Filter1_MatchedJetObservables = subDir_Filter1.mkdir( "matched RECO Jet Observables" );
            histContainer_["Filter1_matchedJetPt"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Pt", "Filter 1 matched Jet Pt; RECO Jet Pt (GeV); Number of Events", 100, 0, 200);
            histContainer_["Filter1_matchedJetEta"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Eta", "Filter 1 matched Jet Eta; RECO Jet Eta; Number of Events", 50, -3, 3);
            histContainer_["Filter1_matchedJetPhi"] = subDir_Filter1_MatchedJetObservables.make<TH1F>("Filter 1 matched Jet Phi", "Filter 1 matched Jet Phi; RECO Jet Phi; Number of Events", 50, -3.5, 3.5);

            subDir_Filter1_TurnOnCurves = subDir_Filter1.mkdir( "Turn On Curves" );
            histContainer_["Filter1_dR"] = subDir_Filter1_TurnOnCurves.make<TH1F>("Filter1_dR","Filter1_dR",100, 0 ,0.010);
            
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

            subDir_Filter2_TurnOnCurves = subDir_Filter2.mkdir( "Turn On Curves" );
            histContainer_["Filter2_dR"] = subDir_Filter2_TurnOnCurves.make<TH1F>("Filter2_dR","Filter2_dR",100, 0 ,0.010);

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

            subDir_Filter3_TurnOnCurves = subDir_Filter3.mkdir( "Turn On Curves" );
            histContainer_["Filter3_dR"] = subDir_Filter3_TurnOnCurves.make<TH1F>("Filter3_dR","Filter2_dR",100, 0 ,0.010);

      }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyser::endJob(){
      //reco vs trigger object

      histContainer_["CondEff"]->Fill(histContainer_["TypeOfEvent"]->GetBinContent(4)/histContainer_["TypeOfEvent"]->GetBinContent(3));
      // histContainer_["CondEff"]->Fill(histContainer_["TypeOfEvent"]->GetBinContent(5)/histContainer_["TypeOfEvent"]->GetBinContent(4));

      // CREATE TRIGGER TURN ON CURVES -------------------------------------------------------------- //


      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_JetMultiplicity_20"],*histContainer_["CrossTrigger_Total_JetMultiplicity_20"]) ){

            TCanvas *c1 = new TCanvas((crosstrigger_ + "_JetMult_20").c_str(),"c1",600,400);
            c1->SetGrid();
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_20"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMultiplicity_20"],histContainer_["CrossTrigger_Total_JetMultiplicity_20"]);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_20"]->SetTitle("Efficiency Plot Jet Multiplicity > 20 GeV; Jet Multiplicity; Efficiency");
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_20"]->SetMarkerColor(4);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_20"]->SetMarkerStyle(21);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_20"]->SetMaximum(1);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_20"]->Draw("AP");
            gPad->Update();
            c1->Update();
            c1->Write();
      }

      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_JetMultiplicity_30"],*histContainer_["CrossTrigger_Total_JetMultiplicity_30"]) ){

            TCanvas *c2 = new TCanvas((crosstrigger_ + "_JetMult_30").c_str(),"c2",600,400);
            c2->SetGrid();
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_30"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMultiplicity_30"],histContainer_["CrossTrigger_Total_JetMultiplicity_30"]);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_30"]->SetTitle("Efficiency Plot Jet Multiplicity > 30 GeV; Jet Multiplicity; Efficiency");
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_30"]->SetMarkerColor(4);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_30"]->SetMarkerStyle(21);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_30"]->SetMaximum(1);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_30"]->Draw("AP");
            gPad->Update();
            c2->Update();
            c2->Write();
      }

      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_JetMultiplicity_40"],*histContainer_["CrossTrigger_Total_JetMultiplicity_40"]) ){

            TCanvas *c3 = new TCanvas((crosstrigger_ + "_JetMult_40").c_str(),"c2",600,400);
            c3->SetGrid();
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_40"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMultiplicity_40"],histContainer_["CrossTrigger_Total_JetMultiplicity_40"]);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_40"]->SetTitle("Efficiency Plot Jet Multiplicity > 40 GeV; Jet Multiplicity; Efficiency");
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_40"]->SetMarkerColor(4);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_40"]->SetMarkerStyle(21);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_40"]->SetMaximum(1);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_40"]->Draw("AP");
            gPad->Update();
            c3->Update();
            c3->Write();
      }

      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_JetMultiplicity_50"],*histContainer_["CrossTrigger_Total_JetMultiplicity_50"]) ){

            TCanvas *c4 = new TCanvas((crosstrigger_ + "_JetMult_50").c_str(),"c4",600,400);
            c4->SetGrid();
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_50"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_JetMultiplicity_50"],histContainer_["CrossTrigger_Total_JetMultiplicity_50"]);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_50"]->SetTitle("Efficiency Plot Jet Multiplicity > 50 GeV; Jet Multiplicity; Efficiency");
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_50"]->SetMarkerColor(4);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_50"]->SetMarkerStyle(21);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_50"]->SetMaximum(1);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_JetMultiplicity_50"]->Draw("AP");
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

      if ( TEfficiency::CheckConsistency(*histContainer_["CrossTrigger_Pass_VertexMultiplicityHist"],*histContainer_["CrossTrigger_Total_VertexMultiplicityHist"]) ){

            TCanvas *c7 = new TCanvas((crosstrigger_ + "_VertMult").c_str(),"c7",600,400);
            c7->SetGrid();
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"] = subDir_TrigDec_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["CrossTrigger_Pass_VertexMultiplicityHist"],histContainer_["CrossTrigger_Total_VertexMultiplicityHist"]);
            turnOnCurveContainer_["CrossTrigger_TurnOnCurve_Vertices"]->SetTitle("Efficiency Plot Vertex Multiplicity; Vertex Multiplicity; Efficiency");
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

            if ( TEfficiency::CheckConsistency(*histContainer_["Filter2_matchedJetBTag"],*histContainer_["Total_RECO_JetBTag"]) ){

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


            if ( TEfficiency::CheckConsistency(*histContainer_["Filter2_matchedJetPt"],*histContainer_["Filter1_matchedJetPt"]) ){

                  TCanvas *c10 = new TCanvas((crosstrigger_ + "_FilterTurnOn_Pt40").c_str(),"c10",600,400);
                  c10->SetGrid();
                  turnOnCurveContainer_["Filter2_TurnOnCurve_Pt"] = subDir_Filter2_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter2_matchedJetPt"],histContainer_["Filter1_matchedJetPt"]);//CHANGED
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


            if ( TEfficiency::CheckConsistency(*histContainer_["Filter3_matchedJetPt"],*histContainer_["Filter2_matchedJetPt"]) ){

                  TCanvas *c11 = new TCanvas((crosstrigger_ + "_FilterTurnOn_Pt50").c_str(),"c11",600,400);
                  c11->SetGrid();
                  turnOnCurveContainer_["Filter3_TurnOnCurve_Pt"] = subDir_Filter3_TurnOnCurves.make<TGraphAsymmErrors>(histContainer_["Filter3_matchedJetPt"],histContainer_["Filter2_matchedJetPt"]);//CHANGED
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


