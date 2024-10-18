#!/usr/bin/env python

import os
import shutil
from argparse import ArgumentParser
import sys
from typing import Optional

import corner
import emcee
import numpy as np
from juliacall import Main as jl
from matplotlib import pyplot as plt
from pint.logging import setup as setup_log
from . import SPNTA, Vela as vl


def parse_args(argv):
    parser = ArgumentParser(
        description="Run single-pulsar timing & noise analysis using Vela.jl"
    )
    parser.add_argument("parfile", help="Pulsar ephemeris file", type=str)
    parser.add_argument("timfile", help="Pulsar TOA file", type=str)
    parser.add_argument(
        "-o", "--outdir", help="Output directory", default=".", type=str, required=True
    )
    parser.add_argument(
        "--cheat_prior_scale",
        help=(
            "Set the parameter prior to be a Uniform distribution centered at the value given in the par file "
            "with a width equal to `cheat_prior_scale * uncertainty_in_parfile`, unless a a default or custom "
            "prior is available"
        ),
        default=10,
        type=float,
    )
    parser.add_argument(
        "-N",
        "--nsteps",
        help="Number of steps for the MCMC sampler",
        default=6000,
        type=int,
    )
    parser.add_argument(
        "-b",
        "--burnin",
        help="Burn-in length for the MCMC chain",
        default=2000,
        type=int,
    )
    parser.add_argument(
        "-t", "--thin", help="Thinning length for the MCMC chain", default=100, type=int
    )
    parser.add_argument(
        "-r", "--resume", help="Resume from a previous run", action="store_true"
    )
    parser.add_argument(
        "-p",
        "--plot_only",
        help="Don't run the sampler, only plot previous results",
        action="store_true",
    )
    parser.add_argument(
        "--no_plot_params",
        help="Exclude these parameters from plotting",
        nargs="*",
        default=[],
    )
    parser.add_argument(
        "--plot_params",
        help="Include these parameters from plotting (supersedes --no_plot_params)",
        nargs="*",
        default=[],
    )

    args = parser.parse_args(argv)

    return args


def prepare_outdir(outdir: str, parfile: str, timfile: str):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    datadir = f"{outdir}/data"
    if not os.path.isdir(datadir):
        os.mkdir(datadir)

    shutil.copy(parfile, datadir)
    shutil.copy(timfile, datadir)


def run_emcee(
    spnta: SPNTA,
    nsteps,
    burnin,
    thin,
    outdir,
    resume,
    plot_only,
):
    nwalkers = spnta.ndim * 5
    p0 = (
        np.array([spnta.prior_transform(cube) for cube in np.random.rand(nwalkers, spnta.ndim)])
        if not resume
        else None
    )

    sampler = emcee.EnsembleSampler(
        nwalkers,
        spnta.ndim,
        spnta.lnpost,
        moves=[emcee.moves.StretchMove(), emcee.moves.DESnookerMove()],
        vectorize=True,
        # parameter_names=param_names,
        backend=emcee.backends.HDFBackend(f"{outdir}/chain.h5"),
    )

    if not plot_only:
        sampler.run_mcmc(p0, nsteps, progress=True, progress_kwargs={"mininterval": 1})

    chain = sampler.get_chain(flat=True, discard=burnin, thin=thin)
    lnprobs = sampler.get_log_prob(flat=True, discard=burnin, thin=thin)

    return chain, lnprobs


def print_results(scaled_samples, shifts, param_names, no_plot_params, plot_params):
    means = np.mean(scaled_samples + shifts, axis=0)
    stds = np.std(scaled_samples, axis=0)

    param_plot_mask = [
        idx
        for idx, par in enumerate(param_names)
        if (
            all(not par.startswith(pnp) for pnp in no_plot_params)
            or any(par.startswith(pnp) for pnp in plot_params)
        )
    ]

    for idx, (pname, mean, std) in enumerate(zip(param_names, means, stds)):
        if idx in param_plot_mask:
            print(f"{pname}\t\t{mean}\t\t{std}")


def plot_corner(
    scaled_samples,
    param_names,
    param_labels,
    pulsar_name,
    no_plot_params,
    plot_params,
    outdir,
):
    param_plot_mask = [
        idx
        for idx, par in enumerate(param_names)
        if (
            all(not par.startswith(pnp) for pnp in no_plot_params)
            or any(par.startswith(pnp) for pnp in plot_params)
        )
    ]

    param_labels = [
        f"\n\n{label}\n\n"
        for idx, label in enumerate(param_labels)
        if idx in param_plot_mask
    ]
    samples_for_plot = scaled_samples[:, param_plot_mask]
    fig = corner.corner(
        samples_for_plot,
        labels=param_labels,
        label_kwargs={"fontsize": 11},
        range=[0.999] * len(param_labels),
        plot_datapoints=False,
        hist_kwargs={"density": True},
    )

    plt.suptitle(pulsar_name)
    plt.tight_layout()
    plt.savefig(f"{outdir}/corner.pdf")
    plt.show()


def plot_postfit_resids(mv, tv, samples, outdir):
    params_median = vl.read_params(mv, np.median(samples, axis=0))
    rv = (
        list(map(vl.value, vl.form_residuals(mv, tv, params_median)))
        if not jl.isa(tv[1], vl.WidebandTOA)
        else [vl.value(wr[0]) for wr in vl.form_residuals(mv, tv, params_median)]
    )

    ctoas = [vl.correct_toa(mv, tvi, params_median) for tvi in tv]
    errs = np.sqrt(
        [vl.value(vl.scaled_toa_error_sqr(tvi, ctoa)) for (tvi, ctoa) in zip(tv, ctoas)]
    )

    mjds = [float(vl.value(tvi.value)) for tvi in tv]

    plt.errorbar(
        mjds,
        rv,
        errs,
        ls="",
        marker="+",
        color="blue",
    )
    plt.axhline(0, color="grey", ls="dotted")
    plt.xlabel("MJD")
    plt.ylabel("Residuals (s)")
    plt.suptitle(mv.pulsar_name)
    plt.tight_layout()
    plt.savefig(f"{outdir}/postfit.pdf")
    plt.show()


def main(argv: Optional[list] = None):
    setup_log(level="WARNING")

    args = parse_args(argv)

    prepare_outdir(args.outdir, args.parfile, args.timfile)

    spnta = SPNTA(args.parfile, args.timfile, cheat_prior_scale=args.cheat_prior_scale)

    samples, lnprs = run_emcee(
        spnta,
        args.nsteps,
        args.burnin,
        args.thin,
        args.outdir,
        args.resume,
        args.plot_only,
    )

    scaled_samples = spnta.rescale_samples(samples)
    shifts = [
        (
            float(spnta.model.param_handler._default_params_tuple.F_.x)
            if pname == "F0"
            else 0
        )
        for pname in spnta.param_names
    ]
    print_results(
        scaled_samples, shifts, spnta.param_names, args.no_plot_params, args.plot_params
    )

    plot_corner(
        scaled_samples,
        spnta.param_names,
        spnta.param_labels,
        spnta.model.pulsar_name,
        args.no_plot_params,
        args.plot_params,
        args.outdir,
    )

    plot_postfit_resids(spnta.model, spnta.toas, samples, args.outdir)


# if __name__ == "__main__":
#     main(sys.argv[1:])
