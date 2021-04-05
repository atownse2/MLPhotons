import sys
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process("mlphotons")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v20'
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1) )

process.source = cms.Source("PoolSource", 
                            fileNames =
                            #cms.untracked.vstring('file:'+sys.argv[2]))
                            cms.untracked.vstring('file:/home/sclark/ml_photons/CMSSW_10_6_14/src/ML_Photons/ml_photons/python/test/aGun_flatMoE_barrel_MINIAODSIM_1.root'))

process.mlphotons = cms.EDProducer(
				'ml_photons',
        classifier_path = cms.FileInPath("ML_Photons/ml_photons/plugins/classifier.onnx"),
        regressor_path = cms.FileInPath("ML_Photons/ml_photons/plugins/regressor.onnx"),
				PhoInputTag = cms.InputTag('slimmedPhotons', '', sys.argv[3]),
				CluInputTag = cms.InputTag('reducedEgamma', 'reducedEBEEClusters', sys.argv[3]),
				HEEInputTag = cms.InputTag('reducedEgamma', 'reducedEERecHits', sys.argv[3]),
				HEBInputTag = cms.InputTag('reducedEgamma', 'reducedEBRecHits', sys.argv[3]),
				RHInputTag = cms.InputTag('reducedEgamma', 'reducedEBRecHits', sys.argv[3]),
        TriggerInputTag_HLT = cms.InputTag('TriggerResults', '', "HLT"),
        VtxInputTag = cms.InputTag('offlineSlimmedPrimaryVertices', '', sys.argv[3]),
        MATCH_DeltaR = cms.double(float(sys.argv[4])),
			)
process.p = cms.Path(process.mlphotons)