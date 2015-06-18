// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/Math/interface/deltaR.h>

#include <TString.h>
#include <string>
#include <TH1D.h>
#include <limits>

//
// class declaration
//

class TriggerAnalyser : public edm::EDAnalyzer {
   public:
      explicit TriggerAnalyser(const edm::ParameterSet&);
      ~TriggerAnalyser();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

// edm::Association< TriggerObjectStandAloneCollection > TriggerObjectStandAloneMatch:

      edm::EDGetTokenT<edm::TriggerResults> triggerResults_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<std::vector<reco::GenJet>> genjets_;
      edm::EDGetTokenT<std::vector<pat::Jet>> jets_;
      edm::EDGetTokenT<std::vector<pat::MET>> mets_;
      edm::EDGetTokenT<std::vector<pat::Electron>> electrons_;
      edm::EDGetTokenT<std::vector<pat::Muon>> muons_;
      const std::string singleleptontrigger_;
      const std::string crosstrigger_;
      const std::string filter1_;
      const std::string filter2_;
      const std::string filter3_;
      const std::string btagger_;
      const std::string hadronicleg_;
      const std::string leptonicleg_;
      std::string CombinedTrigger = "";

      // edm::Association< pat::TriggerObjectStandAloneCollection > TriggerObjectMatch

      edm::Service<TFileService> fileService;
      TH1D *SingleLeptonHist, *CrossTriggerHist, *CrossTriggerCombinedHist;
      TH1D *Filter1_Pt, *Filter1_Eta, *Filter1_Phi;
      TH1D *Filter2_Pt, *Filter2_Eta, *Filter2_Phi;
      TH1D *Filter3_Pt, *Filter3_Eta, *Filter3_Phi;
      TH1D *CrossTrigger_Pass_JetPtHist, *CrossTrigger_Pass_JetEtaHist, *CrossTrigger_Total_JetPtHist, *CrossTrigger_Total_JetEtaHist, *CrossTrigger_Pass_JetMultiplicity, *CrossTrigger_Total_JetMultiplicity, *CrossTrigger_Pass_hltHT, *CrossTrigger_Total_hltHT;
      TH1D *CrossTrigger_Pass_METPtHist, *CrossTrigger_Total_METPtHist, *CrossTrigger_Pass_METEnergyHist, *CrossTrigger_Total_METEnergyHist;
      TH1D *CrossTrigger_Pass_LeptonPtHist, *CrossTrigger_Pass_LeptonEnergyHist, *CrossTrigger_Pass_LeptonEtaHist, *CrossTrigger_Total_LeptonPtHist, *CrossTrigger_Total_LeptonEnergyHist, *CrossTrigger_Total_LeptonEtaHist;
      TH1D *Filter1_matchedJetPt;
      TFileDirectory subDir_TrigDec, subDir_TrigDiffEff, subDir_TrigDiffEff_Jet, subDir_TrigDiffEff_MET, subDir_TrigDiffEff_Lepton, subDir_Filter1, subDir_Filter2, subDir_Filter3;
      TFileDirectory subDir_Filter1_Observables, subDir_Filter1_MatchedJetObservables, subDir_Filter2_Observables, subDir_Filter2_MatchedJetObservables, subDir_Filter3_Observables, subDir_Filter3_MatchedJetObservables;
      // double minDR2 = numeric_limits<double>::infinity();

      bool isJet, isMatched, SingleLeptonTrigDecision, CrossTriggerTrigDecision, CrossTriggerCombinedTrigDecision = false;
      unsigned int crossIndex, singleleptonIndex = 9999;
      float jetPt, jetEta, hltHT = 0;
      float metPt, metEnergy = 0;
      float leptonPt, leptonEta, leptonEnergy = 0;
      int jetMultiplicity, matchedJetIndex = 0;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};




//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
// up to heere can go in .h file?