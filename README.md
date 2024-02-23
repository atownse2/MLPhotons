## Scripts for Clustering, Classifying and Regressing Diphoton Events

Runs over MiniAOD, creates clusters based on hits in ECal. Then applies classification NN to classify clusters as monophoton, diphoton, or hadronic, and then predicts the mass of the clusters.

Also includes scripts for producing flat trees with some objects relevant to an EGamma analysis.

To deploy:

```
cmsrel CMSSW_10_6_19_patch2
cd CMSSW_10_6_19_patch2/src
cmsenv
git clone https://github.com/atownse2/MLPhotons.git
git checkout nano-integration
scram b
```

Get a test file:
```bash
xrdcp root://ndcms.crc.nd.edu//store/mc/RunIISummer20UL18MiniAODv2/BkkToGRadionToGGG_M1-1000_R0-12p5_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/50000/128D50BB-33ED-E742-BAC1-A4F60467D2AA.root ./test.root
```

Test MLPhoton object production:
```bash
cmsRun Prod_MLPhotons.py inputFiles=file:test.root maxEvents=10
```
Test MLPhoton Table production (+NanoAOD steps)
```bash
cmsRun Prod_NanoAOD.py inputFiles=file:test.root maxEvents=10
```
