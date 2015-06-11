import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

# Latest MC Files - Currently 74X Release
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
		'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/022B08C4-C702-E511-9995-D4856459AC30.root'
    )
)

# Input Triggers and Filters
process.Ele27_SingleTop = cms.EDAnalyzer('SingleTop',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    SingleLeptonTriggerInput = cms.string("HLT_Ele27_eta2p1_WP75_Gsf_v1"),
    SingleTopTriggerInput = cms.string("HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1"),
    SingleTopFilterInput = cms.string("hltJetFilterSingleTopEle27"),
    BJetFilterInput = cms.string("hltCSVFilterSingleTop"),
)

process.Ele32_SingleTop = cms.EDAnalyzer('SingleTop',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    SingleLeptonTriggerInput = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_v1"),
    SingleTopTriggerInput = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1"),
    SingleTopFilterInput = cms.string("hltJetFilterSingleTopEle32"),
	BJetFilterInput = cms.string("hltCSVFilterSingleTop"),
)

process.Mu20_SingleTop = cms.EDAnalyzer('SingleTop',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    SingleLeptonTriggerInput = cms.string("HLT_IsoMu20_eta2p1_v1"),
    SingleTopTriggerInput = cms.string("HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v1"),
    SingleTopFilterInput = cms.string("hltJetFilterSingleTopIsoMu20Eta2p1"),
    BJetFilterInput = cms.string("hltCSVFilterSingleTop"),
)

process.Mu24_SingleTop = cms.EDAnalyzer('SingleTop',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    SingleLeptonTriggerInput = cms.string("HLT_IsoMu24_eta2p1_v1"),
    SingleTopTriggerInput = cms.string("HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1"),
    SingleTopFilterInput = cms.string("hltJetFilterSingleTopIsoMu24Eta2p1"),
    BJetFilterInput = cms.string("hltCSVFilterSingleTop"),
)


process.Ele27_TTBarJet30 = cms.EDAnalyzer('TTBarJet30',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    SingleLeptonTriggerInput = cms.string("HLT_Ele27_eta2p1_WP75_Gsf_v1"),
    TTBarJet30TriggerInput = cms.string("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v1"),
    SymmetricJetFilterInput = cms.string("hltEle27TriCentralPFJet30EleCleaned"),
)

process.Ele32_TTBarJet30 = cms.EDAnalyzer('TTBarJet30',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    SingleLeptonTriggerInput = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_v1"),
    TTBarJet30TriggerInput = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet30_v1"),
    SymmetricJetFilterInput = cms.string("hltEle32TriCentralPFJet30EleCleaned"),
)

process.Mu20_TTBarJet30 = cms.EDAnalyzer('TTBarJet30',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    SingleLeptonTriggerInput = cms.string("HLT_IsoMu20_eta2p1_v1"),
    TTBarJet30TriggerInput = cms.string("HLT_IsoMu24_eta2p1_TriCentralPFJet30_v1"),
    SymmetricJetFilterInput = cms.string("hltIsoMu20Eta2p1Trk02TriCentralPFJet30MuCleaned"),
)

process.Mu24_TTBarJet30 = cms.EDAnalyzer('TTBarJet30',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    SingleLeptonTriggerInput = cms.string("HLT_IsoMu24_eta2p1_v1"),
    TTBarJet30TriggerInput = cms.string("HLT_IsoMu24_eta2p1_TriCentralPFJet30_v1"),
    SymmetricJetFilterInput = cms.string("hltIsoMu24Eta2p1Trk02TriCentralPFJet30MuCleaned"),
)


process.Ele27_TTBarJet304050 = cms.EDAnalyzer('TTBarJet304050',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    SingleLeptonTriggerInput = cms.string("HLT_Ele27_eta2p1_WP75_Gsf_v1"),
    TTBarJet304050TriggerInput = cms.string("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1"),
    AsymmetricJet30FilterInput = cms.string("hltEle27TriCentralPFJet30EleCleaned"),
    AsymmetricJet40FilterInput = cms.string("hltEle27DiCentralPFJet40EleCleaned"),
    AsymmetricJet50FilterInput = cms.string("hltEle27CentralPFJet50EleCleaned"),
)

process.Ele32_TTBarJet304050 = cms.EDAnalyzer('TTBarJet304050',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    SingleLeptonTriggerInput = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_v1"),
    TTBarJet304050TriggerInput = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1"),
    AsymmetricJet30FilterInput = cms.string("hltEle32TriCentralPFJet30EleCleaned"),
    AsymmetricJet40FilterInput = cms.string("hltEle32DiCentralPFJet40EleCleaned"),
    AsymmetricJet50FilterInput = cms.string("hltEle32CentralPFJet50EleCleaned"),
)

process.Mu20_TTBarJet304050 = cms.EDAnalyzer('TTBarJet304050',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    SingleLeptonTriggerInput = cms.string("HLT_IsoMu20_eta2p1_v1"),
    TTBarJet304050TriggerInput = cms.string("HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v1"),
    AsymmetricJet30FilterInput = cms.string("hltIsoMu20Eta2p1Trk02TriCentralPFJet30MuCleaned"),
    AsymmetricJet40FilterInput = cms.string("hltIsoMu20Eta2p1Trk02DiCentralPFJet40MuCleaned"),
    AsymmetricJet50FilterInput = cms.string("hltIsoMu20Eta2p1Trk02CentralPFJet50MuCleaned"),
)

process.Mu24_TTBarJet304050 = cms.EDAnalyzer('TTBarJet304050',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTriggerObjects = cms.InputTag('selectedPatTrigger','','PAT'),
    jets = cms.InputTag('slimmedJets','','PAT'),
    SingleLeptonTriggerInput = cms.string("HLT_IsoMu24_eta2p1_v1"),
    TTBarJet304050TriggerInput = cms.string("HLT_IsoMu24_eta2p1_TriCentralPFJet50_40_30_v1"),
    AsymmetricJet30FilterInput = cms.string("hltIsoMu24Eta2p1Trk02TriCentralPFJet30MuCleaned"),
    AsymmetricJet40FilterInput = cms.string("hltIsoMu24Eta2p1Trk02DiCentralPFJet40MuCleaned"),
    AsymmetricJet50FilterInput = cms.string("hltIsoMu24Eta2p1Trk02CentralPFJet50MuCleaned"),
)

# Which trigger analyses are needed
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
    fileName = cms.string('Trigger.root'))