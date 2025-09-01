"""Script to resample `pyvela` output with new burn-in length and thinning factor."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pint.logging import setup as setup_log

from pyvela import pyvela_rethin

setup_log(level="WARNING")


def parse_args(argv):
    parser = ArgumentParser(
        prog="pyvela-rethin",
        description="Resample `pyvela` output with new values of thinning and burn-in.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "result_dir", help="A directory containing the output of the `pyvela` script."
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


def main(argv=None):
    args = parse_args(argv)
    pyvela_rethin.rethin_samples(args.result_dir, burnin=args.burnin, thin=args.thin)
