rm test/out*
#cmsRun Prod_ml_photons.py placeholder PAT test/
cmsRun Prod_ml_photons.py test/test_data16.root DQM test/
cmsRun Prod_flattener.py test/RUCLU_tree_test_data16.root DQM 1 1
