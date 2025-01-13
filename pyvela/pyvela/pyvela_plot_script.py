from typing import Iterable
import numpy as np
import matplotlib.pyplot as plt
import corner
from argparse import ArgumentParser


def parse_args(argv):
    parser = ArgumentParser(
        prog="pyvela-plot",
        description="Plot pyvela results.",
    )
    parser.add_argument("result_dir")
    parser.add_argument("-I", "--ignore_params", nargs="+", default=[])

    return parser.parse_args(argv)


def get_param_plot_mask(
    param_names: Iterable[str], param_prefixes: Iterable[str], args: ArgumentParser
) -> np.ndarray:
    ignore_params_default = {
        "PHOFF",
        "PLREDSIN_",
        "PLREDCOS_",
        "PLDMSIN_",
        "PLDMCOS_",
        "PLCHROMSIN_",
        "PLCHROMCOS_",
        "DMX_",
        "WXSIN_",
        "WXCOS_",
        "DMWXSIN_",
        "DMWXCOS_",
        "CMWXSIN_",
        "CMWXCOS_",
    }
    ignore_params = ignore_params_default.union(args.ignore_params)

    return [
        idx
        for idx, (pname, pprefix) in enumerate(zip(param_names, param_prefixes))
        if (pname not in ignore_params and pprefix not in ignore_params)
    ]


def main(argv=None):
    args = parse_args(argv)

    samples = np.load(f"{args.result_dir}/samples.npy")
    param_names = np.genfromtxt(f"{args.result_dir}/param_names.txt", dtype=str)
    param_prefixes = np.genfromtxt(f"{args.result_dir}/param_prefixes.txt", dtype=str)
    with open(f"{args.result_dir}/param_units.txt", "r") as f:
        param_units = np.array([s.strip() for s in f.readlines()])

    param_plot_mask = get_param_plot_mask(param_names, param_prefixes, args)

    plot_labels = [
        f"{pname}\n{punit}"
        for pname, punit in zip(
            param_names[param_plot_mask], param_units[param_plot_mask]
        )
    ]

    residuals_data = np.genfromtxt(f"{args.result_dir}/residuals.txt")
    wb = residuals_data.shape[1] == 5
    if wb:
        mjds, tres, terr, dres, derr = residuals_data.T
    else:
        mjds, tres, terr = residuals_data.T

    fig = corner.corner(
        samples[:, param_plot_mask],
        labels=plot_labels,
        label_kwargs={"fontsize": 12},
        labelpad=0.02 * len(param_plot_mask),
        max_n_ticks=3,
        plot_datapoints=False,
        hist_kwargs={"density": True},
        range=[0.999] * len(param_plot_mask),
    )

    plt.subplot(5, 4, 4)
    plt.errorbar(mjds, tres, terr, marker="+", ls="")
    plt.axhline(0, ls="dotted")
    plt.ylabel("Time res (s)")
    if wb:
        plt.xticks([])
        plt.subplot(5, 4, 8)
        plt.errorbar(mjds, dres, derr, marker="+", ls="")
        plt.axhline(0, ls="dotted")
        plt.ylabel("DM res (dmu)")
    plt.xlabel("MJD")
    plt.show()
