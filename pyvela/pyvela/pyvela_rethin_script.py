from argparse import ArgumentParser

from pyvela import pyvela_rethin


def parse_args(argv):
    parser = ArgumentParser(
        prog="pyvela-rethin",
        description="Resample emcee output for new values of thinning and burning",
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
