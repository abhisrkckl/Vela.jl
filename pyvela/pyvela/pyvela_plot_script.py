from argparse import ArgumentParser

from pyvela import pyvela_plot


def parse_args(argv):
    parser = ArgumentParser(
        prog="pyvela-plot",
        description="Create a corner plot from pyvela results.",
    )
    parser.add_argument(
        "result_dir", help="A directory containing the output of the `pyvela` script."
    )
    parser.add_argument(
        "-I",
        "--ignore_params",
        nargs="+",
        default=[],
        help="Parameters to exclude from the corner plot.",
    )
    parser.add_argument(
        "-P",
        "--include_params",
        nargs="+",
        default=None,
        help="Parameters to include in the corner plot (will include only these).",
    )
    parser.add_argument(
        "--priors", action="store_true", default=False, help="Plot priors?"
    )
    parser.add_argument("-o", "--out", default=None, help="Output file for plot")
    parser.add_argument(
        "--labelpad", default=0.2, type=float, help="Padding for plot labels"
    )

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    pyvela_plot.plot(
        args.result_dir,
        ignore_params=args.ignore_params,
        include_params=args.include_params,
        out=args.out,
        plot_priors=args.priors,
        labelpad=args.labelpad,
    )
