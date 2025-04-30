import json
from argparse import ArgumentParser
from typing import Iterable

import corner
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from pint import DMconst, dmu
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
        out=args.out,
        labelpad=args.labelpad,
    )
