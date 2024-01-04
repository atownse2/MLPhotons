rm test/out*
#cmsRun Prod_ml_photons.py placeholder PAT test/
cmsRun Prod_ml_photons.py test/test_qcd.root PAT test/
cmsRun Prod_flattener.py test/RUCLU_tree_test_qcd.root PAT 1 1
