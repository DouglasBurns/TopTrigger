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


      // const edm::InputTag hltInputTag_;
      edm::EDGetTokenT<edm::TriggerResults> hltInputTag_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<std::vector<pat::Jet>> jets_;
      const std::string singleleptontrigger_;
      const std::string crosstrigger_;
      const std::string filter1_;
      const std::string filter2_;
      const std::string filter3_;
      const std::string hadronicleg_;

      std::string CombinedTrigger = "";

      edm::Service<TFileService> fileService;
      TH1D *SingleLeptonHist, *CrossTriggerHist, *CrossTriggerCombinedHist;
      TH1D *Filter1_Pt, *Filter1_Eta, *Filter1_Phi;
      TH1D *Filter2_Pt, *Filter2_Eta, *Filter2_Phi;
      TH1D *Filter3_Pt, *Filter3_Eta, *Filter3_Phi;
      TH1D *CrossTrigger_Pass_JetPtHist, *CrossTrigger_Pass_JetEtaHist, *CrossTrigger_Total_JetPtHist, *CrossTrigger_Total_JetEtaHist;

      TFileDirectory subDir_TrigDec, subDir_TrigDiffEff, subDir_Filter1, subDir_Filter2, subDir_Filter3;

      bool SingleLeptonTrigDecision, CrossTriggerTrigDecision, CrossTriggerCombinedTrigDecision = false;
      unsigned int crossIndex, singleleptonIndex = 9999;
      float jetPt, jetEta;

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