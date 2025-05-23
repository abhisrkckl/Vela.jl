from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from pyvela import pyvela_plot


def parse_args(argv):
    parser = ArgumentParser(
        prog="pyvela-plotchains",
        description="Create plots for parameter chains",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "result_dir", help="A directory containing the output of the `pyvela` script."
    )
    parser.add_argument(
        "-o",
        "--out",
        default=None,
        help="Output directory (if different from result_dir)",
    )
    parser.add_argument(
        "-e",
        "--ext",
        default="png",
        choices=["png", "pdf", "gif"],
        help="Extension for plots",
    )

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    pyvela_plot.plot_chains(args.result_dir, outdir=args.out, extension=args.ext)
