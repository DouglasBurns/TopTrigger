from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.requestName = 'Crab_Trigger_Sample'

config.section_('JobType')
config.JobType.psetName = 'ConfFile_cfg.py'
config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'

config.Data.unitsPerJob = 20
config.Data.splitting = 'FileBased'

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_UK_SGrid_Bristol'