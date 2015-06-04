import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'root://xrootd.unl.edu//store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU40bx25_tsg_PHYS14_25_V1-v1/00000/06E41ADB-7870-E411-8850-0025905A605E.root'
		'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/022B08C4-C702-E511-9995-D4856459AC30.root'
    )
)

process.Ele27_SingleTop = cms.EDAnalyzer('SingleTop',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    objects = cms.InputTag("selectedPatTrigger"),
    SingleLeptonTriggerInput = cms.string("HLT_Ele27_eta2p1_WP75_Gsf_v1"),
    SingleTopTriggerInput = cms.string("HLT_Ele27_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1"),
)
process.Ele32_SingleTop = cms.EDAnalyzer('SingleTop',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    objects = cms.InputTag("selectedPatTrigger"),
    SingleLeptonTriggerInput = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_v1"),
    SingleTopTriggerInput = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_CentralPFJet30_BTagCSV07_v1"),
)
process.Mu20_SingleTop = cms.EDAnalyzer('SingleTop',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    objects = cms.InputTag("selectedPatTrigger"),
    SingleLeptonTriggerInput = cms.string("HLT_IsoMu20_eta2p1_v1"),
    SingleTopTriggerInput = cms.string("HLT_IsoMu20_eta2p1_CentralPFJet30_BTagCSV07_v1"),
)
process.Mu24_SingleTop = cms.EDAnalyzer('SingleTop',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    objects = cms.InputTag("selectedPatTrigger"),
    SingleLeptonTriggerInput = cms.string("HLT_IsoMu24_eta2p1_v1"),
    SingleTopTriggerInput = cms.string("HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1"),
)


process.Ele27_TTBarJet30 = cms.EDAnalyzer('TTBarJet30',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    objects = cms.InputTag("selectedPatTrigger"),
    SingleLeptonTriggerInput = cms.string("HLT_Ele27_eta2p1_WP75_Gsf_v1"),
    TTBarJet30TriggerInput = cms.string("HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1"),
)
process.Ele32_TTBarJet30 = cms.EDAnalyzer('TTBarJet30',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    objects = cms.InputTag("selectedPatTrigger"),
    SingleLeptonTriggerInput = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_v1"),
    TTBarJet30TriggerInput = cms.string("HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1"),
)
process.Mu20_TTBarJet30 = cms.EDAnalyzer('TTBarJet30',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    objects = cms.InputTag("selectedPatTrigger"),
    SingleLeptonTriggerInput = cms.string("HLT_IsoMu20_eta2p1_v1"),
    TTBarJet30TriggerInput = cms.string("HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1"),
)
process.Mu24_TTBarJet30 = cms.EDAnalyzer('TTBarJet30',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    objects = cms.InputTag("selectedPatTrigger"),
    SingleLeptonTriggerInput = cms.string("HLT_IsoMu24_eta2p1_v1"),
    TTBarJet30TriggerInput = cms.string("HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1"),
)


process.Ele27_TTBarJet304050 = cms.EDAnalyzer('TTBarJet304050',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    objects = cms.InputTag("selectedPatTrigger"),
    SingleLeptonTriggerInput = cms.string("HLT_Ele27_eta2p1_WP75_Gsf_v1"),
    TTBarJet304050TriggerInput = cms.string("HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1"),
)
process.Ele32_TTBarJet304050 = cms.EDAnalyzer('TTBarJet304050',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    objects = cms.InputTag("selectedPatTrigger"),
    SingleLeptonTriggerInput = cms.string("HLT_Ele32_eta2p1_WP75_Gsf_v1"),
    TTBarJet304050TriggerInput = cms.string("HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1"),
)
process.Mu20_TTBarJet304050 = cms.EDAnalyzer('TTBarJet304050',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    objects = cms.InputTag("selectedPatTrigger"),
    SingleLeptonTriggerInput = cms.string("HLT_IsoMu20_eta2p1_v1"),
    TTBarJet304050TriggerInput = cms.string("HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1"),
)
process.Mu24_TTBarJet304050 = cms.EDAnalyzer('TTBarJet304050',
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    objects = cms.InputTag("selectedPatTrigger"),
    SingleLeptonTriggerInput = cms.string("HLT_IsoMu24_eta2p1_v1"),
    TTBarJet304050TriggerInput = cms.string("HLT_IsoMu24_eta2p1_CentralPFJet30_BTagCSV07_v1"),
)


process.p = cms.Path(process.Ele27_SingleTop*
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

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('TestTrigger.root'))