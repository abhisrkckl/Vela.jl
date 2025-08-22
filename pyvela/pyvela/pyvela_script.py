"""Script for running pulsar timing & noise analysis using Vela.jl with emcee."""

from copy import deepcopy
import json
import os
import shutil
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import emcee
import numpy as np
import astropy.units as u
from scipy.linalg import cholesky, cho_solve, solve_triangular, LinAlgError
from pint.residuals import Residuals, WidebandTOAResiduals

from pyvela import SPNTA
from pyvela.parameters import (
    analytic_marginalizable_names,
    analytic_marginalizable_prefixes,
)


def parse_args(argv):
    parser = ArgumentParser(
        prog="pyvela",
        description="A command line interface for the Vela.jl pulsar timing &"
        " noise analysis package. Uses emcee for sampling. This may not be "
        "appropriate for more complex datasets. Write your own scripts for "
        "such cases.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "par_file",
        help="The pulsar ephemeris file. Should be readable using PINT. The "
        "uncertainties listed in the file will be used for 'cheat' priors where applicable.",
    )
    parser.add_argument(
        "tim_file",
        help="The pulsar TOA file. Should be readable using PINT. Either this or a JLSO file (-J) should be provided.",
    )
    parser.add_argument(
        "-J",
        "--jlso_file",
        help="The JLSO file containing pulsar timing and noise model & TOAs created using "
        "`pyvela-jlso`. JLSO files may need to be recreated after updating `Vela.jl` since "
        "the data format may change. These files are faster to read and parse.",
    )
    parser.add_argument(
        "-P",
        "--prior_file",
        help="A JSON file containing the prior distributions for each free parameter. (Ignored if `-J` option is used.)",
    )
    parser.add_argument(
        "--no_marg_gp_noise",
        action="store_true",
        help="Don't analytically marginalize the correlated Gaussian noise amplitudes.",
    )
    parser.add_argument(
        "-A",
        "--analytic_marg",
        nargs="+",
        default=[],
        help="Parameters to analytically marginalze (only some parameters are allowed).",
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
        default=100,
        type=float,
        help="The scale factor by which the frequentist uncertainties are multiplied to "
        "get the 'cheat' prior distributions.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default="pyvela_results",
        help="The output directory. Will throw an error if it already exists (unless -f is given).",
    )
    parser.add_argument(
        "-f",
        "--force_rewrite",
        action="store_true",
        help="Force rewrite the output directory if it exists.",
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
    parser.add_argument(
        "-r",
        "--resume",
        default=False,
        action="store_true",
        help="Resume from an existing run",
    )
    parser.add_argument(
        "-s",
        "--initial_sample_spread",
        default=0.3,
        type=float,
        help="Spread of the starting samples around the default parameter values. "
        "Must be > 0 and <= 1. 0 represents no spread and 1 represents prior draws.",
    )

    return parser.parse_args(argv)


def validate_input(args):
    assert os.path.isfile(
        args.par_file
    ), f"Invalid par file {args.par_file}. Make sure that the path is correct."
    assert os.path.isfile(
        args.tim_file
    ), f"Invalid tim file {args.tim_file}. Make sure that the path is correct."

    if args.jlso_file is not None:
        assert os.path.isfile(
            args.jlso_file
        ), f"Invalid JLSO file {args.jlso_file}. Make sure that the path is correct."

    assert args.prior_file is None or os.path.isfile(
        args.prior_file
    ), f"Prior file {args.prior_file} not found. Make sure the path is correct."
    assert args.truth is None or os.path.isfile(
        args.truth
    ), f"Truth par file {args.truth} not found.  Make sure the path is correct."

    if args.resume:
        args.force_rewrite = True
    assert args.force_rewrite or not os.path.isdir(
        args.outdir
    ), f"The output directory {args.outdir} already exists! Use `-f` option to force overwrite."

    assert (
        args.initial_sample_spread > 0 and args.initial_sample_spread <= 1
    ), "initial_sample_spread must be > 0 and <= 1."


def prepare_outdir(args):
    if args.force_rewrite and os.path.isdir(args.outdir) and not args.resume:
        shutil.rmtree(args.outdir)

    if not args.resume and not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    if not os.path.exists(f"{args.outdir}/{os.path.basename(args.par_file)}"):
        shutil.copy(args.par_file, args.outdir)

    if not os.path.exists(f"{args.outdir}/{os.path.basename(args.tim_file)}"):
        shutil.copy(args.tim_file, args.outdir)

    if args.jlso_file is not None and not os.path.exists(
        f"{args.outdir}/{os.path.basename(args.jlso_file)}"
    ):
        shutil.copy(args.jlso_file, args.outdir)

    if args.prior_file is not None and not os.path.exists(
        f"{args.outdir}/{os.path.basename(args.prior_file)}"
    ):
        shutil.copy(args.prior_file, args.outdir)

    if args.truth is not None and not os.path.exists(
        f"{args.outdir}/{os.path.basename(args.truth)}"
    ):
        shutil.copy(args.truth, args.outdir)


def main(argv=None):
    args = parse_args(argv)

    if args.resume:
        # copy info from the prior run into the current arguments
        # to make sure they agree
        with open(f"{args.outdir}/summary.json") as summary_file:
            summary_info = json.load(summary_file)
        args.par_file = f'{args.outdir}/{summary_info["input"]["par_file"]}'
        args.tim_file = f'{args.outdir}/{summary_info["input"]["tim_file"]}'
        args.cheat_prior_scale = summary_info["input"]["cheat_prior_scale"]
        args.analytic_marg = summary_info["input"]["analytic_marginalized_params"]
        args.prior_fie = (
            f'{args.outdir}/{summary_info["input"]["custom_prior_file"]}'
            if summary_info["input"]["custom_prior_file"] is not None
            else None
        )
        args.jlso_file = f'{args.outdir}/{summary_info["input"]["jlso_file"]}'

    if "all" in args.analytic_marg:
        assert (
            len(args.analytic_marg) == 1
        ), "Other parameters cannot be specified when `-A all` is given."
        args.analytic_marg = (
            analytic_marginalizable_names + analytic_marginalizable_prefixes
        )

    validate_input(args)

    prepare_outdir(args)

    spnta = (
        SPNTA(
            args.par_file,
            args.tim_file,
            cheat_prior_scale=args.cheat_prior_scale,
            custom_priors=(args.prior_file if args.prior_file is not None else {}),
            marginalize_gp_noise=(not args.no_marg_gp_noise),
            analytic_marginalized_params=args.analytic_marg,
        )
        if args.jlso_file is None
        else SPNTA.load_jlso(
            args.jlso_file,
            args.par_file,
            args.tim_file,
            custom_prior_file=args.prior_file,
            cheat_prior_scale=args.cheat_prior_scale,
            analytic_marginalized_params=args.analytic_marg,
        )
    )

    if not args.resume and spnta.jlsofile is None:
        jlsofile = f"{args.outdir}/_{spnta.model.pulsar_name}.jlso"
        spnta.save_jlso(jlsofile)
        spnta.jlsofile = jlsofile

    nwalkers = spnta.ndim * 5

    spnta.save_pre_analysis_summary(
        args.outdir,
        {
            "sampler": "emcee",
            "nwalkers": nwalkers,
            "nsteps": args.nsteps,
            "burnin": args.burnin,
            "thin": args.thin,
            "vectorized": True,
        },
        args.truth,
    )

    p0 = get_start_samples(spnta, args.initial_sample_spread, nwalkers)

    sampler = emcee.EnsembleSampler(
        nwalkers,
        spnta.ndim,
        spnta.lnpost_vectorized,
        moves=[emcee.moves.StretchMove(), emcee.moves.DESnookerMove()],
        vectorize=True,
        backend=emcee.backends.HDFBackend(f"{args.outdir}/chain.h5"),
    )
    if not args.resume:
        sampler.run_mcmc(
            p0, args.nsteps, progress=True, progress_kwargs={"mininterval": 1}
        )
    else:
        sampler.run_mcmc(
            None, args.nsteps, progress=True, progress_kwargs={"mininterval": 1}
        )

    samples_raw = sampler.get_chain(flat=True, discard=args.burnin, thin=args.thin)

    spnta.save_results(
        args.outdir,
        samples_raw,
    )


def get_start_samples(spnta: SPNTA, s: float, nwalkers: int) -> np.ndarray:
    """Get starting samples for the MCMC. nwalkers is the number of samples
    to be returned.

    The samples are obtained as follows.

        1. Find the maximum-posterior point θ_max=(a_max, b_max).
        2. Compute the design matrix M assuming the timing parameters a_max.
        3. Draw noise parameters b from the prior distribution.
        4. Assuming the noise parameters b, do a linear GLS fit to find the
           timing parameters a.
        5. The sample is ((1-s)*θ_max + s*θ) where θ=(a,b) and 0<s<1.
    """
    m1 = deepcopy(spnta.model_pint)

    # Set timing parameters to the maximum posterior values
    pmax_pint = spnta.maxpost_params / spnta.scale_factors
    for pname, pval in zip(
        spnta.param_names[: spnta.ntmdim], pmax_pint[: spnta.ntmdim]
    ):
        if pname == "F0":
            m1[pname].value += pval
        else:
            m1[pname].value = pval

    # Unfreeze analytically marginalized timing model parameters
    for pname in spnta.analytic_marginalized_params:
        if pname in m1:
            m1[pname].frozen = False

    y = (
        WidebandTOAResiduals(spnta.toas_pint, m1).calc_wideband_resids()
        if spnta.wideband
        else Residuals(spnta.toas_pint, m1).calc_time_resids().to_value(u.s)
    ).astype(float)

    M, params_M = (
        m1.full_wideband_designmatrix(spnta.toas_pint)[:2]
        if spnta.wideband
        else m1.full_designmatrix(spnta.toas_pint)[:2]
    )
    M = M.astype(float)

    ii, iter = 0, 0
    result = np.empty((nwalkers, len(spnta.param_names)))
    while ii < nwalkers:
        print(ii, iter)
        iter += 1
        assert iter <= nwalkers * 10

        # Draw a sample from prior.
        cube = np.random.rand(spnta.ndim)
        p0 = spnta.prior_transform(cube) / spnta.scale_factors

        # Set noise parameters in using prior draws
        for pname, pval in zip(spnta.param_names[spnta.ntmdim :], p0[spnta.ntmdim :]):
            m1[pname].value = pval

        Phidiag = m1.full_basis_weight(spnta.toas_pint).astype(float)
        Ndiag = (
            m1.scaled_wideband_uncertainty(spnta.toas_pint)
            if spnta.wideband
            else m1.scaled_toa_uncertainty(spnta.toas_pint).to_value(u.s)
        ).astype(float)

        Ninv_M = M / Ndiag[:, None]
        MT_Ninv_y = Ninv_M.T @ y
        MT_Ninv_M = M.T @ Ninv_M
        Sigmainv = np.diag(1 / Phidiag) + MT_Ninv_M

        try:
            Sigmainv_cf = cholesky(Sigmainv, lower=False)
            ahat = cho_solve((Sigmainv_cf, False), MT_Ninv_y)
            a1 = ahat
        except LinAlgError:
            continue

        delta_a_dict = dict(zip(params_M, a1))

        delta_param = np.array(
            [
                (
                    delta_a_dict[parname] * scale_factor
                    if parname in delta_a_dict
                    else (p0i - pmaxpost)
                )
                for (parname, p0i, pmaxpost, scale_factor) in zip(
                    spnta.param_names,
                    (p0 * spnta.scale_factors),
                    spnta.maxpost_params,
                    spnta.scale_factors,
                )
            ]
        )

        sample = spnta.maxpost_params + s * delta_param

        if not np.isfinite(spnta.lnpost(sample)):
            continue

        result[ii, :] = sample
        ii += 1

    return result
