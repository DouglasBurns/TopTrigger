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

#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include <TString.h>
#include <string>
#include <TH1D.h>

// Declaration of class
class SingleTop : public edm::EDAnalyzer {
   public:
      explicit SingleTop(const edm::ParameterSet&);
      ~SingleTop();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void StoreJets(const edm::Event&, bool);
      virtual void endJob() override;

      edm::EDGetTokenT<edm::TriggerResults> hltInputTag_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<std::vector<pat::Jet>> jets_;
      const std::string singleleptontrigger_;
      const std::string singletoptrigger_;
      const std::string singletopfilter_;
      const std::string bjetfilter_;
      std::string CombinedTrigger = "";

      edm::Service<TFileService> fileService;
      TH1D *SingleLeptonHist, *SingleTopHist, *SingleTopCombinedHist;
      TH1D *SingleTopFilterHist_Pt, *SingleTopFilterHist_Eta, *SingleTopFilterHist_Phi, *SingleTopFilterHist_Mass;
      TH1D *BJetFilterHist_Pt, *BJetFilterHist_Eta, *BJetFilterHist_Phi, *BJetFilterHist_Mass;
      TH1D *SingleTop_Pass_JetPtHist,*SingleTop_Pass_JetEtaHist,*SingleTop_Total_JetPtHist,*SingleTop_Total_JetEtaHist;

      TFileDirectory subDir_TrigDec, subDir_SingleTopFilter, subDir_BJetFilter, subDir_TrigDiffEff;
      bool SingleLeptonTrigDecision, SingleTopTrigDecision, SingleTopCombinedTrigDecision = false;
      unsigned int singletopIndex, singleleptonIndex = 9999;
      float pT, eta;

};
