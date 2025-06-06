import json
from typing import Iterable, Literal, Optional

try:
    from tqdm import tqdm
except ImportError:
    # allow for systems that don't have tqdm
    def tqdm(iterable, **kwargs):
        return iterable


import corner
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from pint import DMconst, dmu
from pint.models import get_model


def get_param_plot_mask(
    param_names: Iterable[str],
    param_prefixes: Iterable[str],
    ignore_params: Iterable[str] = [],
    include_params: Iterable[str] = None,
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
    ignore_params = ignore_params_default.union(ignore_params)
    if (include_params is None) or (include_params == []):
        return [
            idx
            for idx, (pname, pprefix) in enumerate(zip(param_names, param_prefixes))
            if (pname not in ignore_params and pprefix not in ignore_params)
        ]
    return [
        idx
        for idx, (pname, pprefix) in enumerate(zip(param_names, param_prefixes))
        if (pname in include_params or pprefix in include_params)
    ]


def read_true_values(result_dir):
    with open(f"{result_dir}/summary.json", "r") as summary_file:
        summary = json.load(summary_file)

    if (
        "truth_par_file" not in summary["input"]
        or summary["input"]["truth_par_file"] is None
    ):
        return None

    true_values_raw = np.genfromtxt(f"{result_dir}/param_true_values.txt")
    scale_factors = np.genfromtxt(f"{result_dir}/param_scale_factors.txt")

    return true_values_raw / scale_factors


def get_psrname(result_dir: str) -> Optional[str]:
    try:
        with open(f"{result_dir}/summary.json") as summary_file:
            summary = json.load(summary_file)

        parfile = summary["input"]["par_file"]
        model = get_model(f"{result_dir}/{parfile}")
        return model["PSR"].value
    except:
        return ""


def plot(
    result_dir: str,
    ignore_params: Iterable[str] = [],
    include_params: Iterable[str] = None,
    plot_priors: bool = False,
    outfile: str = None,
    labelpad: float = 0.2,
):
    """Plot `pyvela` output and optionally save it to a file. The output includes a corner plot of the
    posterior samples and the post-fit whitened residuals."""

    samples = np.load(f"{result_dir}/samples.npy")
    param_names = np.genfromtxt(f"{result_dir}/param_names.txt", dtype=str)
    param_prefixes = np.genfromtxt(f"{result_dir}/param_prefixes.txt", dtype=str)
    with open(f"{result_dir}/param_units.txt", "r") as f:
        param_units = np.array([s.strip() for s in f.readlines()])
    scale_factors = np.genfromtxt(f"{result_dir}/param_scale_factors.txt")

    param_plot_mask = get_param_plot_mask(
        param_names, param_prefixes, ignore_params, include_params=include_params
    )

    plot_labels = [
        f"{pname}\n{punit if punit != '1' else ''}"
        for pname, punit in zip(
            param_names[param_plot_mask], param_units[param_plot_mask]
        )
    ]

    residuals_data = np.genfromtxt(f"{result_dir}/residuals.txt")
    wb = residuals_data.shape[1] == 7
    if wb:
        mjds, tres, tres_w, terr, dres, dres_w, derr = residuals_data.T
        dres = (dres * u.Hz / DMconst).to_value(dmu)
        dres_w = (dres_w * u.Hz / DMconst).to_value(dmu)
        derr = (derr * u.Hz / DMconst).to_value(dmu)
    else:
        mjds, tres, tres_w, terr = residuals_data.T

    true_values_all = read_true_values(result_dir)
    true_values = (
        true_values_all[param_plot_mask] if true_values_all is not None else None
    )

    fig = corner.corner(
        samples[:, param_plot_mask],
        labels=plot_labels,
        label_kwargs={"fontsize": 9},
        labelpad=labelpad,
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

    if plot_priors:
        # Plot the pre-evaluated priors
        prior_evals = np.load(f"{result_dir}/prior_evals.npy")
        nplots = len(param_plot_mask)
        for jj, ii in enumerate(param_plot_mask):
            plt.subplot(nplots, nplots, jj * (nplots + 1) + 1)
            xs = prior_evals[:, 2 * ii] / scale_factors[ii]
            ys = prior_evals[:, 2 * ii + 1]
            # normalize the prior to match the plotted histogram
            plt.plot(xs, ys / np.trapz(ys, xs))

    ax = plt.subplot(5, 3, 3)
    ax.errorbar(
        mjds, tres, terr, marker="+", ls="", alpha=1, color="orange", label="Pre-fit"
    )
    ax.set_ylabel("Time res (pre) (s)")
    # ax.legend()

    ax1 = ax.twinx()
    ax1.errorbar([], [], [], ls="", marker="+", color="orange", label="Pre-fit")
    ax1.errorbar(
        mjds,
        tres_w,
        terr,
        marker="+",
        ls="",
        alpha=0.4,
        color="blue",
        label="Post fit whitened",
    )
    ax1.legend()
    ax1.set_ylabel("Time res (post) (s)")
    ax1.axhline(0, ls="dotted", color="k")

    if wb:
        plt.xticks([])
        ax = plt.subplot(5, 3, 6)
        ax.errorbar(
            mjds,
            dres,
            derr,
            marker="+",
            ls="",
            alpha=1,
            color="orange",
            label="Pre-fit",
        )

        ax.set_ylabel("DM res (pre) (dmu)")

        ax1 = ax.twinx()
        ax1.errorbar(
            mjds,
            dres_w,
            derr,
            marker="+",
            ls="",
            alpha=0.4,
            color="blue",
            label="Post fit whitened",
        )
        ax1.set_ylabel("DM res (post) (dmu)")
        ax1.axhline(0, ls="dotted", color="k")

    ax.set_xlabel("MJD - PEPOCH")

    psrname = get_psrname(result_dir)
    if psrname is not None:
        plt.suptitle(psrname, y=0.98, x=0.4)

    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile)


def plot_chains(
    result_dir: str, outdir: str = None, extension: Literal["png", "pdf", "gif"] = "png"
):
    if outdir is None:
        outdir = result_dir
    params = np.loadtxt(f"{result_dir}/param_names.txt", dtype=str)
    d = np.load(f"{result_dir}/samples.npy")
    with open(f"{result_dir}/summary.json") as summary_file:
        summary_info = json.load(summary_file)
    nwalkers = summary_info["sampler"]["nwalkers"]
    thin = summary_info["sampler"]["thin"]

    for i in tqdm(range(len(params))):
        plt.clf()
        plt.plot(np.arange(d.shape[0]), d[:, i], ",")
        ax = plt.gca()
        ax.set_ylabel(params[i])
        ax.set_xlabel("Raw Samples (thinned steps * walkers)")
        ax2 = ax.secondary_xaxis(
            "top",
            functions=(lambda x: x * thin / nwalkers, lambda x: x * nwalkers / thin),
        )
        ax2.set_xlabel("MCMC Steps After Burnin")
        plt.savefig(f"{outdir}/chain_{params[i]}.{extension}")
