## Analysis Code for Clustering, Classifying and Regressing Diphoton Events
Runs over MiniAOD, creates clusters based on hits in ECal. Then applies classification NN to classify clusters as monophoton, diphoton, or hadronic, and then predicts the mass of the clusters.

To run, move this directory into CMSSW_10_6_14/src or later. Then build with 'scram b' and run 'cmsRun ml_photons/python/Prod_ml_photons.py'

Owner: Steven Clark, s.clark@rutgers.edu
