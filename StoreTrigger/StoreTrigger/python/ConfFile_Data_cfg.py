import FWCore.ParameterSet.Config as cms

process = cms.Process("DEMO")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Latest MC Files - Currently 74X Release
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # # DATA Muon
        # 'root://xrootd.unl.edu//store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/168/00000/60FF8405-EA26-E511-A892-02163E01387D.root',
        # # DATA Electron
        # 'root://xrootd.unl.edu//store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/244/00000/12EE24E2-8F27-E511-80D1-02163E013793.root' 


        # DATA Muon 25ns
        # 'root://xrootd.unl.edu//store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/253/890/00000/B0BE6B62-0941-E511-8FBB-02163E014499.root',
        # 'file:StoreTrigger/data/SingleMuon25ns.root',    

        # DATA Electron 25ns 
        'root://xrootd.unl.edu//store/data/Run2015C/SingleElectron/MINIAOD/PromptReco-v1/000/254/227/00000/FC7C42C7-C845-E511-9E5F-02163E012460.root'
        # 'file:StoreTrigger/data/SingleElectron25ns.root',

    )
)  


# Input Triggers and Filters
process.Template = cms.EDAnalyzer('TriggerAnalyser',
    HLTriggerResults = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','RECO'),
    jets = cms.InputTag('slimmedJets','','RECO'),
    mets = cms.InputTag('slimmedMETs','','RECO'),
    electrons = cms.InputTag('slimmedElectrons','','RECO'),
    muons = cms.InputTag('slimmedMuons','','RECO'),
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices','','RECO'),
    SingleLeptonTriggerInput = cms.string(''),
    CrossTriggerInput = cms.string(''),
    FilterInput1 = cms.string(''),
    FilterInput2 = cms.string(''),
    FilterInput3 = cms.string(''),
    LeptonFilterInput = cms.string(''),
    BTagger = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    ElectronID = cms.string('cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight'),
    HadronicLeg = cms.string(''),
    LeptonicLeg = cms.string(''),
    LeptonType = cms.string(''),
)

process.Ele27_SingleTop = process.Template.clone(SingleLeptonTriggerInput = 'HLT_Ele27_eta2p1_WPLoose_Gsf_v1', CrossTriggerInput = 'HLT_Ele27_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v1', FilterInput1 = 'hltJetFilterSingleTopEle27', FilterInput2 = 'hltCSVFilterSingleTop', FilterInput3 = '', LeptonFilterInput = 'hltEle27WPLooseGsfTrackIsoFilter', HadronicLeg  = 'SingleTop', LeptonicLeg = 'Ele', LeptonType = 'Ele27')
process.Ele32_SingleTop = process.Template.clone(SingleLeptonTriggerInput = 'HLT_Ele32_eta2p1_WPLoose_Gsf_v1', CrossTriggerInput = 'HLT_Ele32_eta2p1_WPLoose_Gsf_CentralPFJet30_BTagCSV07_v1', FilterInput1 = 'hltJetFilterSingleTopEle32', FilterInput2 = 'hltCSVFilterSingleTop', FilterInput3 = '', LeptonFilterInput = 'hltEle32WPLooseGsfTrackIsoFilter', HadronicLeg  = 'SingleTop', LeptonicLeg = 'Ele', LeptonType = 'Ele32')
process.Mu20_SingleTop = process.Template.clone(SingleLeptonTriggerInput = 'HLT_IsoMu20_eta2p1_v2', CrossTriggerInput = 'HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v2', FilterInput1 = 'hltJetFilterSingleTopIsoMu20Eta2p1', FilterInput2 = 'hltCSVFilterSingleTop', FilterInput3 = '', LeptonFilterInput = 'hltL3crIsoL1sMu16Eta2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p09', HadronicLeg  = 'SingleTop', LeptonicLeg = 'Mu', LeptonType = 'Mu20')
process.Mu24_SingleTop = process.Template.clone(SingleLeptonTriggerInput = 'HLT_IsoMu24_eta2p1_v2', CrossTriggerInput = 'HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v2', FilterInput1 = 'hltJetFilterSingleTopIsoMu24Eta2p1', FilterInput2 = 'hltCSVFilterSingleTop', FilterInput3 = '', LeptonFilterInput = 'hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09', HadronicLeg  = 'SingleTop', LeptonicLeg = 'Mu', LeptonType = 'Mu24')
process.Ele27_TTBarJet30 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_Ele27_eta2p1_WPLoose_Gsf_v1', CrossTriggerInput = 'HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v1', FilterInput1 = 'hltEle27TriCentralPFJet30EleCleaned', FilterInput2 = '', FilterInput3 = '', LeptonFilterInput = 'hltEle27WPLooseGsfTrackIsoFilter', HadronicLeg  = 'TTBarJet30', LeptonicLeg = 'Ele', LeptonType = 'Ele27')
process.Ele32_TTBarJet30 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_Ele32_eta2p1_WPLoose_Gsf_v1', CrossTriggerInput = 'HLT_Ele32_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v1', FilterInput1 = 'hltEle32TriCentralPFJet30EleCleaned', FilterInput2 = '', FilterInput3 = '', LeptonFilterInput = 'hltEle32WPLooseGsfTrackIsoFilter', HadronicLeg  = 'TTBarJet30', LeptonicLeg = 'Ele', LeptonType = 'Ele32')
process.Mu20_TTBarJet30 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_IsoMu20_eta2p1_v2', CrossTriggerInput = 'HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2', FilterInput1 = 'hltIsoMu20Eta2p1Trk02TriCentralPFJet30MuCleaned', FilterInput2 = '', FilterInput3 = '', LeptonFilterInput = 'hltL3crIsoL1sMu16Eta2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p09', HadronicLeg  = 'TTBarJet30', LeptonicLeg = 'Mu', LeptonType = 'Mu20')
process.Mu24_TTBarJet30 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_IsoMu24_eta2p1_v2', CrossTriggerInput = 'HLT_IsoMu24_eta2p1_TriCentralPFJet30_v2', FilterInput1 = 'hltIsoMu24Eta2p1Trk02TriCentralPFJet30MuCleaned', FilterInput2 = '', FilterInput3 = '', LeptonFilterInput = 'hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09', HadronicLeg  = 'TTBarJet30', LeptonicLeg = 'Mu', LeptonType = 'Mu24')
process.Ele27_TTBarJet304050 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_Ele27_eta2p1_WPLoose_Gsf_v1', CrossTriggerInput = 'HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v1', FilterInput1 = 'hltEle27TriCentralPFJet30EleCleaned', FilterInput2 = 'hltEle27DiCentralPFJet40EleCleaned', FilterInput3 = 'hltEle27CentralPFJet50EleCleaned', LeptonFilterInput = 'hltEle27WPLooseGsfTrackIsoFilter', HadronicLeg  = 'TTBarJet304050', LeptonicLeg = 'Ele', LeptonType = 'Ele27')
process.Ele32_TTBarJet304050 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_Ele32_eta2p1_WPLoose_Gsf_v1', CrossTriggerInput = 'HLT_Ele32_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v1', FilterInput1 = 'hltEle32TriCentralPFJet30EleCleaned', FilterInput2 = 'hltEle32DiCentralPFJet40EleCleaned', FilterInput3 = 'hltEle32CentralPFJet50EleCleaned', LeptonFilterInput = 'hltEle32WPLooseGsfTrackIsoFilter', HadronicLeg  = 'TTBarJet304050', LeptonicLeg = 'Ele', LeptonType = 'Ele32')
process.Mu20_TTBarJet304050 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_IsoMu20_eta2p1_v2', CrossTriggerInput = 'HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v2', FilterInput1 = 'hltIsoMu20Eta2p1Trk02TriCentralPFJet30MuCleaned', FilterInput2 = 'hltIsoMu20Eta2p1Trk02DiCentralPFJet40MuCleaned', FilterInput3 = 'hltIsoMu20Eta2p1Trk02CentralPFJet50MuCleaned', LeptonFilterInput = 'hltL3crIsoL1sMu16Eta2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p09', HadronicLeg  = 'TTBarJet304050', LeptonicLeg = 'Mu', LeptonType = 'Mu20')
process.Mu24_TTBarJet304050 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_IsoMu24_eta2p1_v2', CrossTriggerInput = 'HLT_IsoMu24_eta2p1_TriCentralPFJet50_40_30_v2', FilterInput1 = 'hltIsoMu24Eta2p1Trk02TriCentralPFJet30MuCleaned', FilterInput2 = 'hltIsoMu24Eta2p1Trk02DiCentralPFJet40MuCleaned', FilterInput3 = 'hltIsoMu24Eta2p1Trk02CentralPFJet50MuCleaned', LeptonFilterInput = 'hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09', HadronicLeg  = 'TTBarJet304050', LeptonicLeg = 'Mu', LeptonType = 'Mu24')

process.p = cms.Path(
    process.Ele27_SingleTop*
    process.Ele32_SingleTop*
    process.Mu20_SingleTop*
    process.Mu24_SingleTop*
    process.Ele27_TTBarJet30*
    process.Ele32_TTBarJet30*
    process.Mu20_TTBarJet30*
    process.Mu24_TTBarJet30*
    process.Ele27_TTBarJet304050*
    process.Ele32_TTBarJet304050*
    process.Mu20_TTBarJet304050*
    process.Mu24_TTBarJet304050)

# OutFile
process.TFileService = cms.Service('TFileService',
    # fileName = cms.string('Trigger_Data_Muon.root'))
    # fileName = cms.string('Trigger_Data_Electron.root'))
    fileName = cms.string('Trigger_Data_25ns.root'))
