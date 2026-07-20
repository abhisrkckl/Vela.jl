"""Script for running pulsar timing & noise analysis using Vela.jl with emcee."""

from collections import namedtuple
import os
import shutil
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import pocomc as poco
import numpy as np

from pyvela import SPNTA
from pyvela.parameters import (
    analytic_marginalizable_names,
    analytic_marginalizable_prefixes,
)
from pyvela.results import SPNTAResults


def parse_args(argv):
    parser = ArgumentParser(
        prog="pyvela-poco",
        description="A command line interface for the Vela.jl pulsar timing & "
        "noise analysis package. Uses pocoMC for sampling. This may not be "
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
        "--nsample",
        default=4096,
        type=int,
        help="Effective sample size required before the sampling terminates.",
    )
    parser.add_argument(
        "-r",
        "--resume",
        default=False,
        action="store_true",
        help="Resume from an existing run",
    )
    parser.add_argument(
        "-c",
        "--center_epochs",
        default=False,
        action="store_true",
        help="Center the epochs of the pulsar timing model.",
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

        results = SPNTAResults(args.outdir)

        summary_info = results.summary
        args.par_file = results.input_par_file
        args.tim_file = results.input_tim_file
        args.cheat_prior_scale = summary_info["input"]["cheat_prior_scale"]
        args.analytic_marg = summary_info["input"]["analytic_marginalized_params"]
        args.prior_fie = (
            f'{args.outdir}/{summary_info["input"]["custom_prior_file"]}'
            if summary_info["input"]["custom_prior_file"] is not None
            else None
        )
        args.jlso_file = results.jlso_file
        args.center_epochs = summary_info["input"]["center_epochs"]

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
            marginalize_gp_noise=True,
            analytic_marginalized_params=args.analytic_marg,
            center_epochs=args.center_epochs,
        )
        if args.jlso_file is None
        else SPNTA.load_jlso(
            args.jlso_file,
            args.par_file,
            args.tim_file,
            custom_prior_file=args.prior_file,
            cheat_prior_scale=args.cheat_prior_scale,
            analytic_marginalized_params=args.analytic_marg,
            center_epochs=args.center_epochs,
        )
    )

    if not args.resume and spnta.jlsofile is None:
        jlsofile = f"{args.outdir}/_{spnta.model.pulsar_name}.jlso"
        spnta.save_jlso(jlsofile)
        spnta.jlsofile = jlsofile

    spnta.save_pre_analysis_summary(
        args.outdir,
        {
            "sampler": "pocoMC",
            "nsamples": args.nsamples,
            "vectorized": True,
        },
        args.truth,
    )

    prior = namedtuple("PocoPrior", "dim bounds logpdf rvs")(
        spnta.ndim,
        spnta.prior_bounds,
        spnta.lnprior_vectorized,
        spnta.draw_from_prior,
    )

    sampler = poco.Sampler(
        prior=prior,
        likelihood=spnta.lnlike_vectorized,
        vectorize=True,
        output_dir=args.outdir,
    )
    if not args.resume:
        sampler.run(
            n_total=args.nsamples,
            n_evidence=args.nsamples,
            save_every=10,
            progress=True,
        )
    else:
        pass

    samples_raw, logl, logp = sampler.posterior(resample=True)

    spnta.save_results(
        args.outdir,
        samples_raw,
    )
