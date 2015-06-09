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

class TTBarJet304050 : public edm::EDAnalyzer {
   public:
      explicit TTBarJet304050(const edm::ParameterSet&);
      ~TTBarJet304050();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


      // const edm::InputTag hltInputTag_;
      edm::EDGetTokenT<edm::TriggerResults> hltInputTag_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      const std::string singleleptontrigger_;
      const std::string ttbarjet304050trigger_;
      const std::string asymmetricjet30filter_;
      const std::string asymmetricjet40filter_;
      const std::string asymmetricjet50filter_;
      std::string CombinedTrigger = "";

      edm::Service<TFileService> fileService;
      TH1D *SingleLeptonHist, *TTBarJet304050Hist, *TTBarJet304050CombinedHist;
      TH1D *TTBarJet30Hist_Pt, *TTBarJet30Hist_Eta, *TTBarJet30Hist_Phi;
      TH1D *TTBarJet40Hist_Pt, *TTBarJet40Hist_Eta, *TTBarJet40Hist_Phi;
      TH1D *TTBarJet50Hist_Pt, *TTBarJet50Hist_Eta, *TTBarJet50Hist_Phi;

      TFileDirectory subDir_TrigDec, subDir_AsymmetricJet30Filter, subDir_AsymmetricJet40Filter, subDir_AsymmetricJet50Filter;

      bool SingleLeptonTrigDecision, TTBarJet304050TrigDecision, TTBarJet304050CombinedTrigDecision = false;
      unsigned int ttbarjet304050Index, singleleptonIndex = 9999;
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