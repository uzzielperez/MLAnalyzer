import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register('skipEvents', 
    default=0, 
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.int,
    info = "skipEvents")
options.parseArguments()

process = cms.Process("FEVTAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("DQM.Integration.config.FrontierCondition_GT_Offline_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
#process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet( 
    #input = cms.untracked.int32(1) 
    input = cms.untracked.int32(options.maxEvents) 
    )

print " >> Loaded",len(options.inputFiles),"input files from list."
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
      #'file:myfile.root'
      #'file:pickevents.root'
      #'file:pickeventsRECO.root'
      #'file:step3.root'
      #'file:SinglePhotonPt50_FEVTDEBUG.root'
      #'file:SingleElectronPt50_FEVTDEBUG.root'
      #'file:/eos/uscms/store/user/mba2012/FEVTDEBUG/H125GGgluonfusion_13TeV_TuneCUETP8M1_FEVTDEBUG/*/*/step_full_*.root'
      options.inputFiles
      )
    , skipEvents = cms.untracked.uint32(options.skipEvents)
    )

process.GlobalTag.globaltag = cms.string('80X_dataRun2_HLT_v12')
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.fevt = cms.EDAnalyzer('SCRegressor'
    #, EBRecHitCollection = cms.InputTag('ecalRecHit:EcalRecHitsEB')
    , gsfElectronCollection = cms.InputTag('gedGsfElectrons')
    #, photonCollection = cms.InputTag('gedPhotons')
    , photonCollection = cms.InputTag('slimmedPhotons')
    #, reducedEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
    #, reducedEERecHitCollection = cms.InputTag('reducedEcalRecHitsEE')
    #, reducedESRecHitCollection = cms.InputTag('reducedEcalRecHitsES')
    , reducedEBRecHitCollection = cms.InputTag('reducedEgamma:reducedEBRecHits')
    , reducedEERecHitCollection = cms.InputTag('reducedEgamma:reducedEERecHits')
    , reducedESRecHitCollection = cms.InputTag('reducedEgamma:reducedESRecHits')
    , genParticleCollection = cms.InputTag('genParticles')
    , genJetCollection = cms.InputTag('ak4GenJets')
    , trackCollection = cms.InputTag("generalTracks")
    , rhoLabel = cms.InputTag("fixedGridRhoFastjetAll")
    )

process.TFileService = cms.Service("TFileService",
    #fileName = cms.string('histo.root')
    fileName = cms.string(options.outputFile)
    )

process.hltFilter = cms.EDFilter("HLTHighLevel",
                                          eventSetupPathsKey = cms.string(''),
                                          TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                          HLTPaths = cms.vstring('HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v*'),
                                          andOr = cms.bool(True),
                                          throw = cms.bool(False)
                                          )

#process.p = cms.Path(process.fevt)
process.p = cms.Path(process.hltFilter*process.fevt)
