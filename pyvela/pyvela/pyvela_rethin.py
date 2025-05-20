import json
from typing import Iterable
import os

import corner
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
import emcee
from pint import DMconst, dmu

from pyvela import SPNTA


def rethin_samples(result_dir: str, thin: int, burnin: int):

    reader = emcee.backends.HDFBackend(f"{result_dir}/chain.h5", read_only=True)
    samples_raw = reader.get_chain(flat=True, discard=burnin, thin=thin)

    with open(f"{result_dir}/summary.json") as summary_file:
        summary_info = json.load(summary_file)

    spnta = (
        SPNTA(
            f"{result_dir}/{summary_info['input']['par_file']}",
            f"{result_dir}/{summary_info['input']['tim_file']}",
            cheat_prior_scale=summary_info["input"]["cheat_prior_scale"],
            custom_priors=(
                f"{result_dir}/{summary_info['input']['custom_prior_file']}"
                if summary_info["input"]["custom_prior_file"] is not None
                else {}
            ),
            # marginalize_gp_noise=(not args.no_marg_gp_noise),
        )
        if summary_info["input"]["jlso_file"] is None
        else SPNTA.load_jlso(
            summary_info["input"]["jlso_file"],
            summary_info["input"]["par_file"],
            summary_info["input"]["tim_file"],
        )
    )

    spnta.save_results(
        result_dir,
        samples_raw,
        {
            "sampler": "emcee",
            "nwalkers": summary_info["sampler"]["nwalkers"],
            "nsteps": summary_info["sampler"]["nsteps"],
            "burnin": burnin,
            "thin": thin,
            "vectorized": True,
        },
        None,
    )
