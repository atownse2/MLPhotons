rm test/out*
cmsRun Prod_ml_photons.py placeholder PAT test/
cmsRun Prod_flattener.py test/RUCLU_tree_placeholder 1 1
