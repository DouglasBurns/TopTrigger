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

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <map>
#include <TString.h>
#include <TH1.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TAxis.h>
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
      virtual bool isGoodJet(const pat::Jet& jet) const;

      edm::EDGetTokenT<edm::TriggerResults> triggerResults_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<std::vector<pat::Jet>> jets_;
      edm::EDGetTokenT<std::vector<pat::MET>> mets_;
      edm::EDGetTokenT<std::vector<pat::Electron>> electrons_;
      edm::EDGetTokenT<std::vector<pat::Muon>> muons_;
      edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_;
      const std::string singleleptontrigger_;
      const std::string crosstrigger_;
      const std::string filter1_;
      const std::string filter2_;
      const std::string filter3_;
      const std::string leptonfilter_;
      const std::string btagger_;
      const std::string electronID_;
      const std::string hadronicleg_;
      const std::string leptonicleg_;
      const std::string leptontype_;

      edm::Service<TFileService> fileService;

      std::map<std::string,TH1F*> histContainer_; 
      std::map<std::string,TH2F*> distContainer_; 
      std::map<std::string, TGraphAsymmErrors*> turnOnCurveContainer_;
   
      TFileDirectory subDir_TrigDec, subDir_TrigDec_TurnOnCurves;
      TFileDirectory subDir_Observables, subDir_Observables_Jet, subDir_Observables_Vertices, subDir_Observables_MET, subDir_Observables_Lepton;

      TFileDirectory subDir_Filter1, subDir_Filter2, subDir_Filter3;
      TFileDirectory subDir_Filter1_Observables, subDir_Filter2_Observables, subDir_Filter3_Observables;
      TFileDirectory subDir_Filter1_MatchedJetObservables, subDir_Filter2_MatchedJetObservables, subDir_Filter3_MatchedJetObservables;
      TFileDirectory subDir_Filter1_TurnOnCurves, subDir_Filter2_TurnOnCurves, subDir_Filter3_TurnOnCurves;

      std::vector<int> leadingElectronIndex, leadingMuonIndex, cleanedJetIndex, cleanedBJetIndex, cleanedRECOJetIndex, cleanedRECOJetIndex2;
      bool isJetCollection, isBJetCollection = false;
      bool passesFilter1, passesFilter2, passesFilter3, isMatchedToLepton, isLeadingLepton = false;
      bool SingleLeptonTrigDecision, CrossTriggerTrigDecision = false;
      unsigned int Index, leptonIndex, crossIndex, singleleptonIndex = 9999;
      float jetCSV, forwardjeteta, HT = 0;
      float metPt, metEnergy = 0;
      float leptonPt, leptonEta, leptonEnergy = 0;
      float ptcut, isolation = 0;
      int jetMultiplicity_20, jetMultiplicity_30, jetMultiplicity_40, jetMultiplicity_50 = 0;
      int vertexMultiplicity;
      int no_Jets, no_BJets, no_Muons, no_Electrons = 0;
      int TypeOfEvent = 0;
      bool passMockEventSelection = false;
      bool isolated_muon = true;
      float NHF, NEMF, MUF, CHF, CEMF, eta, pt = 0;
      int NumConst, CHM = 0;
      

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