#!/usr/bin/env python

# %%
from pint.models import get_model_and_toas
from pint.models.priors import Prior
from pint.logging import setup as setup_log
from pint.bayesian import BayesianTiming

import numpy as np
import nestle
import corner
from scipy.stats import uniform
from matplotlib import pyplot as plt
import sys
import time

from pint2vela import read_model_and_toas, Vela as vl

# %%
setup_log(level="WARNING")

# %%
parfile, timfile = sys.argv[1], sys.argv[2]

m, t = get_model_and_toas(parfile, timfile)
t.compute_pulse_numbers(m)

# %%
for par in m.free_params:
    param = getattr(m, par)
    param_min = float(param.value - 10 * param.uncertainty_value)
    param_span = float(20 * param.uncertainty_value)
    param.prior = Prior(uniform(param_min, param_span))

# %%
bt = BayesianTiming(m, t, use_pulse_numbers=True)

# %%
mv, tv = read_model_and_toas(parfile, timfile)
lnlike = vl.get_lnlike_func(mv, tv)
prior_transform = vl.get_prior_transform_func(mv)

# %%
param_names = vl.get_free_param_names(mv.param_handler)
param_idxs = [bt.param_labels.index(p) for p in param_names]
scale_factors = vl.get_scale_factors(mv.param_handler)
shifts = [(float(m.F0.value) if pname == "F0" else 0) for pname in param_names]

# %%
maxlike_params_p = np.array([param.value for param in bt.params], dtype=float)
maxlike_params_v = np.array(
    vl.read_param_values_to_vector(
        mv.param_handler, mv.param_handler._default_params_tuple
    )
)

# %%
# Make sure that the parameter order is OK.
assert np.allclose(
    maxlike_params_p[param_idxs] * scale_factors - shifts, maxlike_params_v
)

# %%
print(lnlike(maxlike_params_v), bt.lnlikelihood(maxlike_params_p))

# %%
# %timeit lnlike_vela(maxlike_params_p)
# %timeit bt.lnlikelihood(maxlike_params_p)

# %%
begin = time.time()
result_nestle_v = nestle.sample(
    lnlike,
    prior_transform,
    bt.nparams,
    method="multi",
    npoints=150,
    dlogz=0.01,
    callback=nestle.print_progress,
)
end = time.time()
print(f"\nTime elapsed = {end-begin} s")

# %%
begin = time.time()
result_nestle_p = nestle.sample(
    bt.lnlikelihood,
    bt.prior_transform,
    bt.nparams,
    method="multi",
    npoints=150,
    dlogz=0.01,
    callback=nestle.print_progress,
)
end = time.time()
print(f"\nTime elapsed = {end-begin} s")

samples_v = (result_nestle_v.samples + shifts) / scale_factors

# %%
fig = corner.corner(
    result_nestle_p.samples[:, param_idxs],
    weights=result_nestle_p.weights,
    labels=np.array(bt.param_labels)[param_idxs],
    label_kwargs={"fontsize": 15},
    range=[0.999] * bt.nparams,
)
corner.corner(
    samples_v,
    weights=result_nestle_v.weights,
    labels=bt.param_labels,
    label_kwargs={"fontsize": 15},
    range=[0.999] * bt.nparams,
    fig=fig,
    color="red",
)
plt.suptitle(m.PSR.value)
plt.tight_layout()
plt.show()
