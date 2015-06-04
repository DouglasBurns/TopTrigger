// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

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

      edm::Service<TFileService> fileService;
      TH1D *SingleLeptonHist, *SingleTopHist;
      
      bool SingleLeptonTrigDecision, SingleTopTrigDecision = false;
      std::string singleleptontrigger, singletoptrigger = "Trigger";
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