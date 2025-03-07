#!/usr/bin/env python

from pyvela import SPNTA

datasets = ["NGC6440E", "sim1", "sim2", "sim_sw.wb"]

for dataset in datasets:
    print(f"Processing {dataset} ...")

    parfile = f"../pyvela/examples/{dataset}.par"
    timfile = f"../pyvela/examples/{dataset}.tim"
    prior_file = f"../pyvela/tests/datafiles/custom_priors.json"
    spnta = SPNTA(parfile, timfile, cheat_prior_scale=50, custom_priors=prior_file)

    spnta.save_jlso(f"datafiles/{dataset}.jlso")