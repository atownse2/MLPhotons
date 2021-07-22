rm test/out*
cmsRun Prod_ml_photons.py placeholder PAT test/
cmsRun Prod_flattenerMatching.py test/RUCLU_tree_placeholder PAT 1 1
