import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Latest MC Files - Currently 74X Release
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
		# MiniAODSIM 25
        'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/022B08C4-C702-E511-9995-D4856459AC30.root'
        # MiniAODSIM 50
        # 'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/00466730-F801-E511-9594-549F35AF450A.root'
        # Aligning
        # 'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v4/10000/00D2A247-2910-E511-9F3D-0CC47A4DEDD2.root'
        # AODSIM
        # 'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/00B2F04A-E301-E511-8CF3-0025905A6090.root'
        # GEN-SIM-RAW
        # 'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/GEN-SIM-RAW/AsymptFlat10to50bx25Raw_MCRUN2_74_V9-v1/10000/00BA7065-3F12-E511-AAB5-00259073E36E.root'
    )
)  


# Input Triggers and Filters
process.Template = cms.EDAnalyzer('TriggerAnalyser',
    HLTriggerResults = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    mets = cms.InputTag('slimmedMETs','','PAT'),
    electrons = cms.InputTag('slimmedElectrons','','PAT'),
    muons = cms.InputTag('slimmedMuons','','PAT'),
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices','','PAT'),
    SingleLeptonTriggerInput = cms.string(''),
    CrossTriggerInput = cms.string(''),
    FilterInput1 = cms.string(''),
    FilterInput2 = cms.string(''),
    FilterInput3 = cms.string(''),
    BTagger = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    HadronicLeg = cms.string(''),
    LeptonicLeg = cms.string(''),
)

process.Ele27_SingleTop = process.Template.clone(SingleLeptonTriggerInput = 'HLT_Ele27_eta2p1_WP75_Gsf_v1', CrossTriggerInput = 'HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1', FilterInput1 = 'hltJetFilterSingleTopEle27', FilterInput2 = 'hltCSVFilterSingleTop', FilterInput3 = '', HadronicLeg  = 'SingleTop', LeptonicLeg = 'Ele')
process.Ele32_SingleTop = process.Template.clone(SingleLeptonTriggerInput = 'HLT_Ele32_eta2p1_WP75_Gsf_v1', CrossTriggerInput = 'HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1', FilterInput1 = 'hltJetFilterSingleTopEle32', FilterInput2 = 'hltCSVFilterSingleTop', FilterInput3 = '', HadronicLeg  = 'SingleTop', LeptonicLeg = 'Ele')
process.Mu20_SingleTop = process.Template.clone(SingleLeptonTriggerInput = 'HLT_IsoMu20_eta2p1_v1', CrossTriggerInput = 'HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v1', FilterInput1 = 'hltJetFilterSingleTopIsoMu20Eta2p1', FilterInput2 = 'hltCSVFilterSingleTop', FilterInput3 = '', HadronicLeg  = 'SingleTop', LeptonicLeg = 'Mu')
process.Mu24_SingleTop = process.Template.clone(SingleLeptonTriggerInput = 'HLT_IsoMu24_eta2p1_v1', CrossTriggerInput = 'HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1', FilterInput1 = 'hltJetFilterSingleTopIsoMu24Eta2p1', FilterInput2 = 'hltCSVFilterSingleTop', FilterInput3 = '', HadronicLeg  = 'SingleTop', LeptonicLeg = 'Mu')
process.Ele27_TTBarJet30 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_Ele27_eta2p1_WP75_Gsf_v1', CrossTriggerInput = 'HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v1', FilterInput1 = 'hltEle27TriCentralPFJet30EleCleaned', FilterInput2 = '', FilterInput3 = '', HadronicLeg  = 'TTBarJet30', LeptonicLeg = 'Ele')
process.Ele32_TTBarJet30 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_Ele32_eta2p1_WP75_Gsf_v1', CrossTriggerInput = 'HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet30_v1', FilterInput1 = 'hltEle32TriCentralPFJet30EleCleaned', FilterInput2 = '', FilterInput3 = '', HadronicLeg  = 'TTBarJet30', LeptonicLeg = 'Ele')
process.Mu20_TTBarJet30 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_IsoMu20_eta2p1_v1', CrossTriggerInput = 'HLT_IsoMu20_eta2p1_TriCentralPFJet30_v1', FilterInput1 = 'hltIsoMu20Eta2p1Trk02TriCentralPFJet30MuCleaned', FilterInput2 = '', FilterInput3 = '', HadronicLeg  = 'TTBarJet30', LeptonicLeg = 'Mu')
process.Mu24_TTBarJet30 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_IsoMu24_eta2p1_v1', CrossTriggerInput = 'HLT_IsoMu24_eta2p1_TriCentralPFJet30_v1', FilterInput1 = 'hltIsoMu24Eta2p1Trk02TriCentralPFJet30MuCleaned', FilterInput2 = '', FilterInput3 = '', HadronicLeg  = 'TTBarJet30', LeptonicLeg = 'Mu')
process.Ele27_TTBarJet304050 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_Ele27_eta2p1_WP75_Gsf_v1', CrossTriggerInput = 'HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1', FilterInput1 = 'hltEle27TriCentralPFJet30EleCleaned', FilterInput2 = 'hltEle27DiCentralPFJet40EleCleaned', FilterInput3 = 'hltEle27CentralPFJet50EleCleaned', HadronicLeg  = 'TTBarJet304050', LeptonicLeg = 'Ele')
process.Ele32_TTBarJet304050 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_Ele32_eta2p1_WP75_Gsf_v1', CrossTriggerInput = 'HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1', FilterInput1 = 'hltEle32TriCentralPFJet30EleCleaned', FilterInput2 = 'hltEle32DiCentralPFJet40EleCleaned', FilterInput3 = 'hltEle32CentralPFJet50EleCleaned', HadronicLeg  = 'TTBarJet304050', LeptonicLeg = 'Ele')
process.Mu20_TTBarJet304050 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_IsoMu20_eta2p1_v1', CrossTriggerInput = 'HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v1', FilterInput1 = 'hltIsoMu20Eta2p1Trk02TriCentralPFJet30MuCleaned', FilterInput2 = 'hltIsoMu20Eta2p1Trk02DiCentralPFJet40MuCleaned', FilterInput3 = 'hltIsoMu20Eta2p1Trk02CentralPFJet50MuCleaned', HadronicLeg  = 'TTBarJet304050', LeptonicLeg = 'Mu')
process.Mu24_TTBarJet304050 = process.Template.clone(SingleLeptonTriggerInput = 'HLT_IsoMu24_eta2p1_v1', CrossTriggerInput = 'HLT_IsoMu24_eta2p1_TriCentralPFJet50_40_30_v1', FilterInput1 = 'hltIsoMu24Eta2p1Trk02TriCentralPFJet30MuCleaned', FilterInput2 = 'hltIsoMu24Eta2p1Trk02DiCentralPFJet40MuCleaned', FilterInput3 = 'hltIsoMu24Eta2p1Trk02CentralPFJet50MuCleaned', HadronicLeg  = 'TTBarJet304050', LeptonicLeg = 'Mu')

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
    fileName = cms.string('TestTrigger.root'))