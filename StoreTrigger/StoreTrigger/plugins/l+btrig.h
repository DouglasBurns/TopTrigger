// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include <DataFormats/JetReco/interface/Jet.h>
// #include <DataFormats/HLTReco/interface/TriggerTypeDefs.h>
// #include "DataFormats/HLTReco/interface/TriggerEvent.h"
// #include "DataFormats/HLTReco/interface/TriggerObject.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
// #include <DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h>
// #include <FWCore/Utilities/interface/EDGetToken.h>
// #include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include <TString.h>
#include <string>
#include <TH1D.h>

//
// class declaration
//

class SingleTop : public edm::EDAnalyzer {
   public:
      explicit SingleTop(const edm::ParameterSet&);
      ~SingleTop();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


      // const edm::InputTag hltInputTag_;
      edm::EDGetTokenT<edm::TriggerResults> hltInputTag_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      const std::string singleleptontrigger_;
      const std::string singletoptrigger_;
      const std::string singletopfilter_;
      const std::string bjetfilter_;

      edm::Service<TFileService> fileService;
      TH1D *SingleLeptonHist, *SingleTopHist;
      TH1D *SingleTopFilterHist_Pt, *SingleTopFilterHist_Eta, *SingleTopFilterHist_Phi, *SingleTopFilterHist_Mass;
      TH1D *BJetFilterHist_Pt, *BJetFilterHist_Eta, *BJetFilterHist_Phi, *BJetFilterHist_Mass;

      TFileDirectory subDir_TrigDec, subDir_SingleTopFilter, subDir_BJetFilter;
      bool SingleLeptonTrigDecision, SingleTopTrigDecision = false;
      unsigned int singletopIndex, singleleptonIndex = 9999;

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