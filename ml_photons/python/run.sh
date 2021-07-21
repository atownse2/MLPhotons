rm test/out*
cmsRun Prod_ml_photons.py placeholder DQM test/
cmsRun Prod_flattener.py test/RUCLU_tree_placeholder DQM 1 1
