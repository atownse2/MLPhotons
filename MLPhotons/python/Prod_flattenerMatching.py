import sys
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
import os

pwd = os.getcwd()
pwd = pwd[:pwd.rfind("/")]

process = cms.Process("flattenerMatching")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v20'
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

inpath = sys.argv[2]
midfilename = sys.argv[3]
useXRD = sys.argv[4]=='useXRD'

prefix = 'file:'

redir = "ndcms.crc.nd.edu"
#redir = "cmsxrootd.fnal.gov"
if useXRD:
  prefix = 'root://'+redir+'/' 

process.source = cms.Source("PoolSource", 
                            fileNames =
                            #cms.untracked.vstring(*inputMC))
                            cms.untracked.vstring(prefix+inpath))
                            #cms.untracked.vstring('file:' + pwd + '/python/test/test_qcd.root'))


process.TFileService = cms.Service("TFileService",
        fileName = cms.string(midfilename)
      )

process.flattenerMatching = cms.EDAnalyzer(
				'flattenerMatching',
        TriggerInputTag_HLT = cms.InputTag('TriggerResults', '', "HLT"),

        genpartInputTag = cms.InputTag('prunedGenParticles', '', sys.argv[3]),
        patjetInputTag = cms.InputTag('slimmedJets', '', sys.argv[3]),
        metInputTag = cms.InputTag('slimmedMETs', '', sys.argv[3]),
        muonInputTag = cms.InputTag('slimmedMuons', '', sys.argv[3]),
        electronInputTag = cms.InputTag('slimmedElectrons', '', sys.argv[3]),
        patPhoInputTag = cms.InputTag('slimmedPhotons', '', sys.argv[3]),
        pvtxInputTag = cms.InputTag('offlineSlimmedPrimaryVertices', '', sys.argv[3]),
        svtxInputTag = cms.InputTag('slimmedSecondaryVertices', '', sys.argv[3]),


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

        tfile_path = cms.string(pwd + "/plugins/tnames16.txt"),
        btag_name = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"), # See: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation#Different_SF_campaigns

        weightInput = cms.double(float(sys.argv[4]) / float(sys.argv[5]))

			)


process.p = cms.Path(process.flattenerMatching)
