rm test/out*
cmsRun Prod_ml_photons.py placeholder RECO test/
cmsRun Prod_flattener.py test/RUCLU_tree_placeholder RECO 1 1
