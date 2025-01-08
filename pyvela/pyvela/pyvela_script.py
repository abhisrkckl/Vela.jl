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

import pyvela
from pyvela import SPNTA
from pyvela import Vela as vl


def info_dict(args):
    return {
        "input": {
            "par_file": args.par_file,
            "tim_file": args.tim_file,
            "prior_file": args.prior_file,
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


def parse_args(argv):
    parser = ArgumentParser(
        prog="pyvela",
        description="A command line interface for the Vela.jl pulsar timing & noise analysis package",
    )
    parser.add_argument("par_file")
    parser.add_argument("tim_file")
    parser.add_argument("-P", "--prior_file")
    parser.add_argument("--cheat_prior_scale", default=10, type=float)
    parser.add_argument("-o", "--outdir", default="pyvela_results")
    parser.add_argument("-N", "--nsteps", default=6000, type=int)
    parser.add_argument("-b", "--burnin", default=1500, type=int)
    parser.add_argument("-t", "--thin", default=100, type=int)

    return parser.parse_args(argv)


def prepare_outdir(args):
    summary_info = info_dict(args)

    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    shutil.copy(args.par_file, args.outdir)
    shutil.copy(args.tim_file, args.outdir)
    if args.prior_file is not None:
        shutil.copy(args.prior_file, args.outdir)
    with open(f"{args.outdir}/summary.json", "w") as summary_file:
        json.dump(summary_info, summary_file, indent=4)


def save_spnta_attrs(spnta: SPNTA, args):
    np.savetxt(f"{args.outdir}/param_names.txt", spnta.param_names, fmt="%s")
    np.savetxt(f"{args.outdir}/param_units.txt", spnta.param_units, fmt="%s")
    np.savetxt(f"{args.outdir}/param_scale_factors.txt", spnta.scale_factors, fmt="%s")


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
