#!/usr/bin/env python

# %%
import sys
from timeit import timeit

import corner
import emcee
import numpy as np
from matplotlib import pyplot as plt
from pint.logging import setup as setup_log
from pint.models import get_model_and_toas

from pyvela import Vela as vl
from pyvela import read_model_and_toas

# %%
setup_log(level="WARNING")

# %%
parfile, timfile = sys.argv[1], sys.argv[2]
m, t = get_model_and_toas(parfile, timfile)

# %%
mv, tv = read_model_and_toas(parfile, timfile, cheat_prior_scale=10)
lnpost = vl.get_lnpost_func(mv, tv, True)
prior_transform = vl.get_prior_transform_func(mv)

# %%
param_names = vl.get_free_param_names(mv.param_handler)
scale_factors = vl.get_scale_factors(mv.param_handler)
shifts = [(float(m.F0.value) if pname == "F0" else 0) for pname in param_names]

# %%
maxlike_params_v = np.array(
    [
        vl.read_param_values_to_vector(
            mv.param_handler, mv.param_handler._default_params_tuple
        )
    ]
)

# %%
print(lnpost(maxlike_params_v))
print(timeit("lnpost(maxlike_params_v)", globals=globals(), number=1000))

# %%
ndim = len(param_names)
nwalkers = ndim * 5
p0 = np.array([prior_transform(cube) for cube in np.random.rand(nwalkers, ndim)])

sampler = emcee.EnsembleSampler(
    nwalkers,
    ndim,
    lnpost,
    moves=[emcee.moves.StretchMove(), emcee.moves.DESnookerMove()],
    vectorize=True,
)
sampler.run_mcmc(p0, 6000, progress=True, progress_kwargs={"mininterval": 1})

# %%
samples_v_0 = sampler.get_chain(flat=True, discard=2500, thin=40)
samples_v = samples_v_0 / scale_factors

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
    for idx, par in enumerate(param_names)
    if all(not par.startswith(pnp) for pnp in params_no_plot)
]

# %%
means = (np.mean(samples_v_0, axis=0) + shifts) / scale_factors
stds = np.std(samples_v, axis=0)
for idx, (pname, mean, std) in enumerate(zip(param_names, means, stds)):
    if idx in param_plot_mask:
        print(f"{pname}\t\t{mean}\t\t{std}")

# %%
# param_labels = [f"\n\n{pname}\n({m[pname].units})\n" for pname in param_names]
param_labels = [
    f"\n\n{label}\n\n"
    for idx, label in enumerate(vl.get_free_param_labels(mv))
    if idx in param_plot_mask
]
samples_for_plot = samples_v[:, param_plot_mask]
fig = corner.corner(
    samples_for_plot,
    labels=param_labels,
    label_kwargs={"fontsize": 11},
    range=[0.999] * len(param_labels),
    truths=(maxlike_params_v[0] / scale_factors)[param_plot_mask],
    plot_datapoints=False,
    hist_kwargs={"density": True},
)

plt.suptitle(m.PSR.value)
plt.tight_layout()
plt.show()

# %%
params_median = vl.read_params(mv, np.median(samples_v_0, axis=0))
rv = (
    list(map(vl.value, vl.form_residuals(mv, tv, params_median)))
    if not t.is_wideband()
    else [vl.value(wr[0]) for wr in vl.form_residuals(mv, tv, params_median)]
)

ctoas = [vl.correct_toa(mv, tvi, params_median) for tvi in tv]
errs = np.sqrt(
    [vl.value(vl.scaled_toa_error_sqr(tvi, ctoa)) for (tvi, ctoa) in zip(tv, ctoas)]
)

plt.errorbar(
    t.get_mjds(),
    rv,
    errs,
    ls="",
    marker="+",
    color="blue",
)
plt.axhline(0, color="grey", ls="dotted")
plt.xlabel("MJD")
plt.ylabel("Residuals (s)")
plt.suptitle(m.PSR.value)
plt.tight_layout()
plt.show()
