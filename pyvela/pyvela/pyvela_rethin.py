import json

import emcee
from pyvela import SPNTA


def rethin_samples(result_dir: str, thin: int, burnin: int):
    """Change the burn-in length and thinning fraction for an MCMC chain created using the
    `pyvela` script."""
    reader = emcee.backends.HDFBackend(f"{result_dir}/chain.h5", read_only=True)
    samples_raw = reader.get_chain(flat=True, discard=burnin, thin=thin)

    with open(f"{result_dir}/summary.json") as summary_file:
        summary_info = json.load(summary_file)

    spnta = SPNTA.load_jlso(
        f"{result_dir}/{summary_info["input"]["jlso_file"]}",
        f"{result_dir}/{summary_info["input"]["par_file"]}",
        f"{result_dir}/{summary_info["input"]["tim_file"]}",
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
        f"{result_dir}/{summary_info["input"]["truth_par_file"]}",
    )
