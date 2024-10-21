#!/usr/bin/env python

# %%
from timeit import timeit

import corner
import emcee
import numpy as np
from matplotlib import pyplot as plt
from argparse import ArgumentParser

from pyvela import SPNTA, Vela as vl
from pyvela.priors import parse_custom_prior_file


parser = ArgumentParser()
parser.add_argument("par_file")
parser.add_argument("tim_file")
parser.add_argument("-P", "--priors", default={})
args = parser.parse_args()

# %%
spnta = SPNTA(
    args.par_file, args.tim_file, cheat_prior_scale=10, custom_priors=args.priors
)

# %%
shifts = [
    (spnta.model.param_handler._default_params_tuple.F_.x if pname == "F0" else 0)
    for pname in spnta.param_names
]

# %%
maxlike_params_v = np.array([spnta.maxlike_params])

# %%
print(spnta.lnpost_vectorized(maxlike_params_v))
print(
    timeit("spnta.lnpost_vectorized(maxlike_params_v)", globals=globals(), number=1000)
)

# %%
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

# %%
samples_v_0 = sampler.get_chain(flat=True, discard=2500, thin=40)
samples_v = spnta.rescale_samples(samples_v_0)

# %%
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

# %%
means = (np.mean(samples_v_0, axis=0) + shifts) / spnta.scale_factors
stds = np.std(samples_v, axis=0)
for idx, (pname, mean, std) in enumerate(zip(spnta.param_names, means, stds)):
    if idx in param_plot_mask:
        print(f"{pname}\t\t{mean:.18e}\t\t{std:.4e}")

# %%
# param_labels = [f"\n\n{pname}\n({m[pname].units})\n" for pname in param_names]
param_labels = [
    f"\n\n{label}\n\n"
    for idx, label in enumerate(spnta.param_labels)
    if idx in param_plot_mask
]
samples_for_plot = samples_v[:, param_plot_mask]
fig = corner.corner(
    samples_for_plot,
    labels=param_labels,
    label_kwargs={"fontsize": 11},
    range=[0.999] * len(param_labels),
    truths=(maxlike_params_v[0] / spnta.scale_factors)[param_plot_mask],
    plot_datapoints=False,
    hist_kwargs={"density": True},
)

plt.suptitle(spnta.model.pulsar_name)
plt.tight_layout()
plt.show()

# %%
params_median = vl.read_params(spnta.model, np.median(samples_v_0, axis=0))
rv = (
    list(map(vl.value, vl.form_residuals(spnta.model, spnta.toas, params_median)))
    if not spnta.is_wideband()
    else [
        vl.value(wr[0])
        for wr in vl.form_residuals(spnta.model, spnta.toas, params_median)
    ]
)

ctoas = [vl.correct_toa(spnta.model, tvi, params_median) for tvi in spnta.toas]
errs = np.sqrt(
    [
        vl.value(vl.scaled_toa_error_sqr(tvi, ctoa))
        for (tvi, ctoa) in zip(spnta.toas, ctoas)
    ]
)

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
plt.tight_layout()
plt.show()
