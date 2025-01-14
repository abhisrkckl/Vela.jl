from copy import deepcopy
import datetime
import getpass
import json
import os
import platform
import shutil
import sys
from argparse import ArgumentParser

import emcee
import numpy as np
import pint

import pint.models
import pyvela
from pyvela import SPNTA
from pyvela import Vela as vl


def info_dict(args):
    info_dict = {
        "input": {
            "par_file": os.path.basename(args.par_file),
            "tim_file": os.path.basename(args.tim_file),
            "prior_file": os.path.basename(args.prior_file),
            "cheat_prior_scale": args.cheat_prior_scale,
        },
        "sampler": {
            "nsteps": args.nsteps,
            "burnin": args.burnin,
            "thin": args.thin,
        },
        "env": {
            "launch_time": datetime.datetime.now().isoformat(),
            "user": getpass.getuser(),
            "host": platform.node(),
            "os": platform.platform(),
            "julia_threads": vl.nthreads(),
            "python": sys.version,
            "julia": str(vl.VERSION),
            "pyvela": pyvela.__version__,
            "pint": pint.__version__,
            "emcee": emcee.__version__,
        },
    }

    if args.truth is not None:
        info_dict["input"]["truth_par_file"] = os.path.basename(args.truth)

    return info_dict


def parse_args(argv):
    parser = ArgumentParser(
        prog="pyvela",
        description="A command line interface for the Vela.jl pulsar timing &"
        " noise analysis package. Uses emcee for sampling. This may not be "
        "appropriate for more complex datasets. Write your own scripts for "
        "such cases.",
    )
    parser.add_argument(
        "par_file",
        help="The pulsar ephemeris file. Should be readable using PINT. The "
        "uncertainties listed in the file will be used for 'cheat' priors where applicable.",
    )
    parser.add_argument(
        "tim_file", help="The pulsar TOA file. Should be readable using PINT."
    )
    parser.add_argument(
        "-P",
        "--prior_file",
        help="A JSON file containing the prior distributions for each free parameter.",
    )
    parser.add_argument(
        "-T",
        "--truth",
        help="Pulsar ephemeris file containing the true timing and noise parameter values. "
        "Relevant for simulation studies.",
    )
    parser.add_argument(
        "-C",
        "--cheat_prior_scale",
        default=50,
        type=float,
        help="The scale factor by which the frequentist uncertainties are multiplied to "
        "get the 'cheat' prior distributions.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default="pyvela_results",
        help="The output directory. Will throw an error if it already exists.",
    )
    parser.add_argument(
        "-N",
        "--nsteps",
        default=6000,
        type=int,
        help="Number of ensemble MCMC iterations",
    )
    parser.add_argument(
        "-b",
        "--burnin",
        default=1500,
        type=int,
        help="Burn-in length for MCMC chains",
    )
    parser.add_argument(
        "-t",
        "--thin",
        default=100,
        type=int,
        help="Thinning factor for MCMC chains",
    )

    return parser.parse_args(argv)


def prepare_outdir(args):
    summary_info = info_dict(args)

    assert not os.path.isdir(args.outdir), "The output directory already exists!"

    os.mkdir(args.outdir)
    shutil.copy(args.par_file, args.outdir)
    shutil.copy(args.tim_file, args.outdir)
    if args.prior_file is not None:
        shutil.copy(args.prior_file, args.outdir)
    with open(f"{args.outdir}/summary.json", "w") as summary_file:
        json.dump(summary_info, summary_file, indent=4)

    if args.truth is not None:
        shutil.copy(args.truth, args.outdir)


def get_true_values(spnta: SPNTA, args):
    true_model = pint.models.get_model(args.truth)
    return (
        np.array(
            [
                true_model[par].value if par in true_model else np.nan
                for par in spnta.param_names
            ]
        )
        * spnta.scale_factors
    )


def save_spnta_attrs(spnta: SPNTA, args):
    np.savetxt(f"{args.outdir}/param_names.txt", spnta.param_names, fmt="%s")
    np.savetxt(f"{args.outdir}/param_prefixes.txt", spnta.param_prefixes, fmt="%s")
    np.savetxt(f"{args.outdir}/param_units.txt", spnta.param_units, fmt="%s")
    np.savetxt(f"{args.outdir}/param_scale_factors.txt", spnta.scale_factors)

    if args.truth is not None:
        np.savetxt(f"{args.outdir}/param_true_values.txt", get_true_values(spnta, args))


def save_new_parfile(
    spnta: SPNTA, params: np.ndarray, param_uncertainties: np.ndarray, filename: str
):
    param_vals = spnta.rescale_samples(params)
    param_errs = spnta.rescale_samples(param_uncertainties)

    model1 = deepcopy(spnta.model_pint)
    for pname, pval, perr in zip(spnta.param_names, param_vals, param_errs):
        model1[pname].value = (
            pval
            if pname != "F0"
            else (
                np.longdouble(spnta.model.param_handler._default_params_tuple.F_.x)
                + pval
            )
        )
        model1[pname].uncertainty_value = perr

    model1.write_parfile(filename)


def save_resids(spnta: SPNTA, params: np.ndarray, outdir: str) -> None:
    wb = spnta.is_wideband()

    ntoas = len(spnta.toas)
    mjds = spnta.get_mjds()
    tres = spnta.time_residuals(params)
    terr = spnta.scaled_toa_unceritainties(params)

    res_arr = np.zeros((ntoas, 2 * (1 + int(wb)) + 1))
    res_arr[:, 0] = mjds
    res_arr[:, 1] = tres
    res_arr[:, 2] = terr

    if wb:
        dres = spnta.dm_residuals(params)
        derr = spnta.scaled_dm_unceritainties(params)

        res_arr[:, 3] = dres
        res_arr[:, 4] = derr

    np.savetxt(f"{outdir}/residuals.txt", res_arr)


def main(argv=None):
    args = parse_args(argv)
    prepare_outdir(args)

    spnta = SPNTA(
        args.par_file,
        args.tim_file,
        cheat_prior_scale=args.cheat_prior_scale,
        custom_priors=(args.prior_file if args.prior_file is not None else {}),
    )

    save_spnta_attrs(spnta, args)

    nwalkers = spnta.ndim * 5
    p0 = np.array(
        [spnta.prior_transform(cube) for cube in np.random.rand(nwalkers, spnta.ndim)]
    )

    sampler = emcee.EnsembleSampler(
        nwalkers,
        spnta.ndim,
        spnta.lnpost_vectorized,
        # moves=[emcee.moves.StretchMove(), emcee.moves.DESnookerMove()],
        vectorize=True,
        backend=emcee.backends.HDFBackend(f"{args.outdir}/chain.h5"),
    )
    sampler.run_mcmc(p0, args.nsteps, progress=True, progress_kwargs={"mininterval": 1})
    samples_raw = sampler.get_chain(flat=True, discard=args.burnin, thin=args.thin)
    samples = spnta.rescale_samples(samples_raw)

    with open(f"{args.outdir}/samples_raw.npy", "wb") as f:
        np.save(f, samples_raw)
    with open(f"{args.outdir}/samples.npy", "wb") as f:
        np.save(f, samples)

    param_uncertainties = np.std(samples_raw, axis=0)

    params_maxpost = sampler.chain[
        *np.unravel_index(
            np.argmax(sampler.get_log_prob().T), sampler.get_log_prob().T.shape
        ),
        :,
    ]
    np.savetxt(f"{args.outdir}/params_maxpost.txt", params_maxpost)
    save_new_parfile(
        spnta,
        params_maxpost,
        param_uncertainties,
        f"{args.outdir}/{spnta.model.pulsar_name}.maxpost.par",
    )

    params_median = np.median(samples_raw, axis=0)
    np.savetxt(f"{args.outdir}/params_median.txt", params_median)
    save_new_parfile(
        spnta,
        params_median,
        param_uncertainties,
        f"{args.outdir}/{spnta.model.pulsar_name}.median.par",
    )

    save_resids(spnta, params_median, args.outdir)

    np.savetxt(f"{args.outdir}/param_default_values.txt", spnta.default_params)
