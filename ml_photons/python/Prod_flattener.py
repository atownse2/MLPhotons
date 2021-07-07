import sys
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
import os

pwd = os.getcwd()
pwd = pwd[:pwd.rfind("/")]

process = cms.Process("flattener")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v20'
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1) )

process.source = cms.Source("PoolSource", 
                            fileNames =
                            cms.untracked.vstring('file:'+sys.argv[2])) #Local File
                            #cms.untracked.vstring("root://cmsxrootd.fnal.gov//"+sys.argv[2])) #Official MC (or other nonlocal file)
                            #cms.untracked.vstring('file:' + pwd + '/python/test/test_rucluAOD.root')) #Test file

infname = sys.argv[2]
fname_only = infname.split("/")[-1]
outfname = "/cms/sclark/RUCLU_Outputs/Data/2017/flat/" + sys.argv[2][infname.find("Run"):infname.rfind("/")+1] + fname_only.replace("RUCLU","flat")
print(outfname)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string(outfname)
      )

process.flattener = cms.EDAnalyzer(
				'flattener',
        TriggerInputTag_HLT = cms.InputTag('TriggerResults', '', "HLT"),

        patjetInputTag = cms.InputTag('slimmedJetsPuppi', '', "PAT"),#TODO: PAT might be different for data, qcd ,etc. Make arg
        metInputTag = cms.InputTag('slimmedMETsPuppi', '', "PAT"),#TODO: PAT might be different for data, qcd ,etc. Make arg
        muonInputTag = cms.InputTag('slimmedMuons', '', "PAT"),#TODO: PAT might be different for data, qcd ,etc. Make arg
        pvtxInputTag = cms.InputTag('offlineSlimmedPrimaryVertices', '', "PAT"),
        svtxInputTag = cms.InputTag('slimmedSecondaryVertices', '', "PAT"),


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

        weightInput = cms.double(float(sys.argv[3]) / float(sys.argv[4]))

			)


process.p = cms.Path(process.flattener)
