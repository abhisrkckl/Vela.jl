#!/usr/bin/env python

# %%
from argparse import ArgumentParser
from timeit import timeit

import corner
import emcee
import numpy as np
from matplotlib import pyplot as plt

from pyvela import SPNTA
from pyvela import Vela as vl

parser = ArgumentParser()
parser.add_argument("par_file")
parser.add_argument("tim_file")
parser.add_argument("-P", "--priors", default={})
args = parser.parse_args()


spnta = SPNTA(
    args.par_file, args.tim_file, cheat_prior_scale=10, custom_priors=args.priors
)


shifts = [
    (spnta.model.param_handler._default_params_tuple.F_.x if pname == "F0" else 0)
    for pname in spnta.param_names
]


maxlike_params_v = np.array([spnta.maxlike_params])


print(spnta.lnpost_vectorized(maxlike_params_v))
print(
    timeit("spnta.lnpost_vectorized(maxlike_params_v)", globals=globals(), number=1000)
)


nwalkers = spnta.ndim * 5
p0 = np.array(
    [spnta.prior_transform(cube) for cube in np.random.rand(nwalkers, spnta.ndim)]
)

sampler = emcee.EnsembleSampler(
    nwalkers,
    spnta.ndim,
    spnta.lnpost_vectorized,
    moves=[emcee.moves.StretchMove(), emcee.moves.DESnookerMove()],
    vectorize=True,
)
sampler.run_mcmc(p0, 6000, progress=True, progress_kwargs={"mininterval": 1})


samples_v_0 = sampler.get_chain(flat=True, discard=2500, thin=40)
samples_v = spnta.rescale_samples(samples_v_0)


params_no_plot = [
    "PLREDSIN_",
    "PLREDCOS_",
    "PLDMSIN_",
    "PLDMCOS_",
    "PLCHROMSIN_",
    "PLCHROMCOS_",
]
param_plot_mask = [
    idx
    for idx, par in enumerate(spnta.param_names)
    if all(not par.startswith(pnp) for pnp in params_no_plot)
]


means = (np.mean(samples_v_0, axis=0) + shifts) / spnta.scale_factors
stds = np.std(samples_v, axis=0)
for idx, (pname, mean, std) in enumerate(zip(spnta.param_names, means, stds)):
    if idx in param_plot_mask:
        print(f"{pname}\t\t{mean:.18e}\t\t{std:.4e}")

params_median = np.median(samples_v_0, axis=0)
rv = spnta.time_residuals(params_median)
errs = spnta.scaled_toa_unceritainties(params_median)

samples_for_plot = samples_v[:, param_plot_mask]
fig = corner.corner(
    samples_for_plot,
    labels=np.array(spnta.param_labels)[param_plot_mask],
    label_kwargs={"fontsize": 11},
    range=[0.999] * spnta.ndim,
    truths=(maxlike_params_v[0] / spnta.scale_factors)[param_plot_mask],
    plot_datapoints=False,
    hist_kwargs={"density": True},
    labelpad=0.3,
    max_n_ticks=3,
)

plt.suptitle(spnta.model.pulsar_name)

plt.subplot(6, 3, 3)
plt.errorbar(
    spnta.get_mjds(),
    rv,
    errs,
    ls="",
    marker="+",
    color="blue",
)
plt.axhline(0, color="grey", ls="dotted")
plt.xlabel("MJD")
plt.ylabel("Residuals (s)")
plt.suptitle(spnta.model.pulsar_name)

plt.show()
