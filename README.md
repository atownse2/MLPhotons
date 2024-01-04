## Analysis Code for Clustering, Classifying and Regressing Diphoton Events

Runs over MiniAOD, creates clusters based on hits in ECal. Then applies classification NN to classify clusters as monophoton, diphoton, or hadronic, and then predicts the mass of the clusters.

Also includes scripts for producing flat trees with some objects relevant to an EGamma analysis.

To deploy:

```
cmsrel CMSSW_10_6_19_patch2
cd CMSSW_10_6_19_patch2/src
cmsenv
git clone https://github.com/atownse2/MLPhotons.git
scram b
```

You can run a test job by doing:
```
cd MLPhotons/MLPhotons/python
cmsRun Prod_FlatAOD.py --inputFiles /scratch365/atownse2/tmp/test_MiniAOD.root --outputFile test.root --maxEvents 10 --isMC True
```
