import sys
import os

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

CMSSW_BASE = os.environ['CMSSW_BASE']

options = VarParsing('analysis')
options.register('year',
                 '2018',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 'Year to process')

options.parseArguments()

inputFiles = options.inputFiles
outputFile = options.outputFile
maxEvents = options.maxEvents


process = cms.Process("mlphotons")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v20'
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(maxEvents))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(*inputFiles))

process.mlphotons = cms.EDProducer(
    'ml_photons',
    classifier_path = cms.string(CMSSW_BASE + "/src/MLPhotons/MLPhotons/plugins/classifier.onnx"),
    regressor_path = cms.string(CMSSW_BASE + "/src/MLPhotons/MLPhotons/plugins/regressor.onnx"),
    CluInputTag = cms.InputTag('reducedEgamma', 'reducedEBEEClusters', 'PAT'),
    HEEInputTag = cms.InputTag('reducedEgamma', 'reducedEERecHits', 'PAT'),
    HEBInputTag = cms.InputTag('reducedEgamma', 'reducedEBRecHits', 'PAT'),
    pfcandInputTag = cms.InputTag('packedPFCandidates', '', 'PAT'),
    cluster_name = cms.string("RUCLUs"),
    # PhoInputTag = cms.InputTag('slimmedPhotons', '', 'PAT'),
    # VtxInputTag = cms.InputTag('offlineSlimmedPrimaryVertices', '', 'PAT'),
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(outputFile),
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