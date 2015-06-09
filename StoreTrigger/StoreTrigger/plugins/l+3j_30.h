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

class TTBarJet30 : public edm::EDAnalyzer {
   public:
      explicit TTBarJet30(const edm::ParameterSet&);
      ~TTBarJet30();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


      // const edm::InputTag hltInputTag_;
      edm::EDGetTokenT<edm::TriggerResults> hltInputTag_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      const std::string singleleptontrigger_;
      const std::string ttbarjet30trigger_;
      const std::string symmetricjetfilter_;
      std::string CombinedTrigger = "";

      edm::Service<TFileService> fileService;
      TH1D *SingleLeptonHist, *TTBarJet30Hist, *TTBarJet30CombinedHist;
      TH1D *TTBarJet30Hist_Pt, *TTBarJet30Hist_Eta, *TTBarJet30Hist_Phi;

      TFileDirectory subDir_TrigDec, subDir_SymmetricJetFilter;
      bool SingleLeptonTrigDecision, TTBarJet30TrigDecision, TTBarJet30CombinedTrigDecision = false;
      unsigned int ttbarjet30Index, singleleptonIndex = 9999;

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