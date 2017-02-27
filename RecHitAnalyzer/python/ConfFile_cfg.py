import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("DQM.Integration.config.FrontierCondition_GT_Offline_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi");

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
		# replace 'myfile.root' with the source file you want to use
		fileNames = cms.untracked.vstring(
			#'file:myfile.root'
			#'file:pickevents.root'
			'file:pickeventsRECO.root'
			)
		)

process.GlobalTag.globaltag = cms.string('80X_dataRun2_HLT_v2')
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.demo = cms.EDAnalyzer('RecHitAnalyzer'
		#, tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
		, EBRecHitCollection = cms.InputTag('ecalRecHit:EcalRecHitsEB')
		, reducedEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
		, selectedEBDigiCollection = cms.InputTag('selectDigi:selectedEcalEBDigiCollection')
		, reducedHBHERecHitCollection = cms.InputTag('reducedHcalRecHits:hbhereco')
		)

process.TFileService = cms.Service("TFileService",
		fileName = cms.string('histo.root')
		)

process.p = cms.Path(process.demo)
