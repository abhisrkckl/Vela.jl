from typing import Iterable
import numpy as np
import matplotlib.pyplot as plt
import corner
import json
from argparse import ArgumentParser
from astropy import units as u
from pint import dmu, DMconst


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


def read_true_values(args):
    with open(f"{args.result_dir}/summary.json", "r") as summary_file:
        summary = json.load(summary_file)

    if "truth_par_file" not in summary["input"]:
        return None

    true_values_raw = np.genfromtxt(f"{args.result_dir}/param_true_values.txt")
    scale_factors = np.genfromtxt(f"{args.result_dir}/param_scale_factors.txt")

    return true_values_raw / scale_factors


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
        dres = (dres * u.Hz / DMconst).to_value(dmu)
        derr = (derr * u.Hz / DMconst).to_value(dmu)
    else:
        mjds, tres, terr = residuals_data.T

    true_values = read_true_values(args)[param_plot_mask]

    fig = corner.corner(
        samples[:, param_plot_mask],
        labels=plot_labels,
        label_kwargs={"fontsize": 9},
        labelpad=0.04 * len(param_plot_mask),
        max_n_ticks=3,
        plot_datapoints=False,
        hist_kwargs={"density": True},
        range=[0.999] * len(param_plot_mask),
        truths=true_values,
    )

    for ax in fig.get_axes():
        ax.tick_params(axis="both", labelsize=8)
        ax.yaxis.get_offset_text().set_fontsize(8)
        ax.xaxis.get_offset_text().set_fontsize(8)

    plt.subplot(5, 3, 3)
    plt.errorbar(mjds, tres, terr, marker="+", ls="", alpha=0.9)
    plt.axhline(0, ls="dotted", color="k")
    plt.ylabel("Time res (s)")
    if wb:
        plt.xticks([])
        plt.subplot(5, 3, 6)
        plt.errorbar(mjds, dres, derr, marker="+", ls="", alpha=0.9)
        plt.axhline(0, ls="dotted", color="k")
        plt.ylabel("DM res (dmu)")
    plt.xlabel("MJD - PEPOCH")
    plt.show()
