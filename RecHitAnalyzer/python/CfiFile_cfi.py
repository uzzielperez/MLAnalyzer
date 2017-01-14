import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('RecHitAnalyzer'
     ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
