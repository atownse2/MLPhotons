import sys
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
import os

pwd = os.getcwd()
pwd = pwd[:pwd.rfind("/")]

process = cms.Process("mlphotons")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v20'
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100) )

process.source = cms.Source("PoolSource", 
                            fileNames =
                            #cms.untracked.vstring('file:'+sys.argv[2])) //Local File
                            #cms.untracked.vstring("root://cmsxrootd.fnal.gov//"+sys.argv[2])) #Official MC (or other nonlocal file)
                            cms.untracked.vstring('file:' + pwd + '/python/test/test_data18.root')) #Test file

process.mlphotons = cms.EDProducer(
				'ml_photons',
        classifier_path = cms.string(pwd + "/plugins/classifier.onnx"),
        regressor_path = cms.string(pwd + "/plugins/regressor.onnx"),
				PhoInputTag = cms.InputTag('slimmedPhotons', '', sys.argv[3]),
				CluInputTag = cms.InputTag('reducedEgamma', 'reducedEBEEClusters', sys.argv[3]),
				HEEInputTag = cms.InputTag('reducedEgamma', 'reducedEERecHits', sys.argv[3]),
				HEBInputTag = cms.InputTag('reducedEgamma', 'reducedEBRecHits', sys.argv[3]),
        genpartInputTag = cms.InputTag('prunedGenParticles', '', sys.argv[3]),
				RHInputTag = cms.InputTag('reducedEgamma', 'reducedEBRecHits', sys.argv[3]),
        TriggerInputTag_HLT = cms.InputTag('TriggerResults', '', "HLT"),
        VtxInputTag = cms.InputTag('offlineSlimmedPrimaryVertices', '', sys.argv[3]),
        cluster_name = cms.string("RUCLUs"),
			)

infname = sys.argv[2][(sys.argv[2].rfind("/") + 1 ) :]
outfname = sys.argv[4] + 'RUCLU_tree_' + infname
#outfname = 'file:' + pwd + '/python/test/output_test.root'

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(outfname),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    compressionAlgorithm = cms.untracked.string("LZMA"),
    compressionLevel = cms.untracked.int32(4),
    outputCommands = cms.untracked.vstring('keep *',
     "drop *_gtStage2Digis_*_*",
     "drop *_caloStage2Digis_*_*",
     "drop *_gmtStage2Digis_*_*",
     "drop *_hcalnoise_*_*",
     "drop *_gtDigis_*_*",
     "drop *_fixedGridRho*_*_*",
     "drop *_scalersRawToDigi_*_*",
     "drop *_l1extraParticles_*_*",
     "drop *_bunchSpacingProducer_*_*",
     "drop *_BeamHaloSummary_*_*",
     "drop *_CSCHaloData_*_*",
     )
    )

process.p = cms.Path(process.mlphotons)

process.e = cms.EndPath(process.out)
