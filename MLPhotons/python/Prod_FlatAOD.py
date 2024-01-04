import sys
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018


# Get CMSSW_BASE path
import os
CMSSW_BASE = os.environ['CMSSW_BASE']


options = VarParsing('analysis')

options.register('isMC',
                 True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'Is this MC?')
options.register('year',
                 '2018',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 'Year to process')

options.parseArguments()

inputFiles = options.inputFiles
outputFile = options.outputFile
maxEvents = options.maxEvents
isMC = options.isMC
year = options.year

if year != '2018':
    raise ValueError('Only 2018 data is supported at this time')

# Store all trigger names in a text file, not sure if this is the best way to do this
triggers = '/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton/analysis/metadata/triggers/triggerNames_{}.txt'.format(year)

if not os.path.isfile(triggers):
    print("Trigger file does not exist: {}".format(triggers))
    sys.exit()

process = cms.Process("mlphotons")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v20'
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(maxEvents))

process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(*inputFiles))

process.TFileService = cms.Service("TFileService",
        fileName = cms.string(outputFile)
      )

process.mlphotons = cms.EDProducer(
				'ml_photons',
        classifier_path = cms.string(CMSSW_BASE + "/src/MLPhotons/MLPhotons/classifier.onnx"),
        regressor_path = cms.string(CMSSW_BASE + "/src/MLPhotons/MLPhotons/regressor.onnx"),
				PhoInputTag = cms.InputTag('slimmedPhotons', '', 'PAT'),
				CluInputTag = cms.InputTag('reducedEgamma', 'reducedEBEEClusters', 'PAT'),
				HEEInputTag = cms.InputTag('reducedEgamma', 'reducedEERecHits', 'PAT'),
				HEBInputTag = cms.InputTag('reducedEgamma', 'reducedEBRecHits', 'PAT'),
        genpartInputTag = cms.InputTag('prunedGenParticles', '', 'PAT'),
				RHInputTag = cms.InputTag('reducedEgamma', 'reducedEBRecHits', 'PAT'),
        TriggerInputTag_HLT = cms.InputTag('TriggerResults', '', "HLT"),
        VtxInputTag = cms.InputTag('offlineSlimmedPrimaryVertices', '', 'PAT'),
        cluster_name = cms.string("RUCLUs"),
			)
process.p = cms.Path(process.mlphotons)

process.flattener = cms.EDAnalyzer(
        'flattener',
        TriggerInputTag_HLT = cms.InputTag('TriggerResults', '', "HLT"),
        genpartInputTag = cms.InputTag('prunedGenParticles', '', 'PAT'),
        patjetInputTag = cms.InputTag('slimmedJets', '', 'PAT'),
        pfcandInputTag = cms.InputTag('packedPFCandidates', '', 'PAT'),
        metInputTag = cms.InputTag('slimmedMETs', '', 'PAT'),
        muonInputTag = cms.InputTag('slimmedMuons', '', 'PAT'),
        electronInputTag = cms.InputTag('slimmedElectrons', '', 'PAT'),
        patPhoInputTag = cms.InputTag('slimmedPhotons', '', 'PAT'),
        pvtxInputTag = cms.InputTag('offlineSlimmedPrimaryVertices', '', 'PAT'),
        svtxInputTag = cms.InputTag('slimmedSecondaryVertices', '', 'PAT'),

        ruclu_etaTag = cms.InputTag('mlphotons', 'RUCLUsEta', "mlphotons"),
        ruclu_phiTag = cms.InputTag('mlphotons', 'RUCLUsPhi', "mlphotons"),
        ruclu_energyTag = cms.InputTag('mlphotons', 'RUCLUsE', "mlphotons"),
        ruclu_r1Tag = cms.InputTag('mlphotons', 'RUCLUsR1', "mlphotons"),
        ruclu_r2Tag = cms.InputTag('mlphotons', 'RUCLUsR2', "mlphotons"),
        ruclu_r3Tag = cms.InputTag('mlphotons', 'RUCLUsR3', "mlphotons"),
        diphoInputTag = cms.InputTag('mlphotons', 'RUCLUsDipho', "mlphotons"),
        monophoInputTag = cms.InputTag('mlphotons', 'RUCLUsMonopho', "mlphotons"),
        hadronInputTag = cms.InputTag('mlphotons', 'RUCLUsHadron', "mlphotons"),
        moeInputTag = cms.InputTag('mlphotons', 'RUCLUsMoE', "mlphotons"),

        tfile_path = cms.string(triggers),
        btag_name = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"), # See: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation#Different_SF_campaigns
        isMC = cms.bool(isMC),

        weightInput = cms.double(1/1)

      )
process.f = cms.Path(process.flattener)

process.schedule = cms.Schedule(process.p, process.f)