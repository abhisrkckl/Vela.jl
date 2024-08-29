#!/usr/bin/env python

# %%
from pint.models import get_model_and_toas
from pint.logging import setup as setup_log

import numpy as np
import emcee
import corner
from matplotlib import pyplot as plt
import sys
from timeit import timeit

from pint2vela import read_model_and_toas, Vela as vl

# %%
setup_log(level="WARNING")

# %%
parfile, timfile = sys.argv[1], sys.argv[2]
m, t = get_model_and_toas(parfile, timfile)

# %%
mv, tv = read_model_and_toas(parfile, timfile)
lnpost = vl.get_lnpost_func(mv, tv)
prior_transform = vl.get_prior_transform_func(mv)

# %%
param_names = vl.get_free_param_names(mv.param_handler)
scale_factors = vl.get_scale_factors(mv.param_handler)
shifts = [(float(m.F0.value) if pname == "F0" else 0) for pname in param_names]

# %%
maxlike_params_v = np.array(
    vl.read_param_values_to_vector(
        mv.param_handler, mv.param_handler._default_params_tuple
    )
)

# %%
print(lnpost(maxlike_params_v))
print(timeit("lnpost(maxlike_params_v)", globals=globals(), number=1000))

# %%
ndim = len(param_names)
nwalkers = 100
p0 = np.array([prior_transform(cube) for cube in np.random.rand(nwalkers, ndim)])

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnpost)
sampler.run_mcmc(p0, 2500, progress=True)

# %%
samples_v_0 = sampler.get_chain(flat=True, discard=500, thin=10)
samples_v = samples_v_0 / scale_factors

# %%
means = (np.mean(samples_v_0, axis=0) + shifts) / scale_factors
stds = np.std(samples_v, axis=0)
for pname, mean, std in zip(param_names, means, stds):
    print(f"{pname}\t\t{mean}\t\t{std}")

# %%
# param_labels = [f"\n\n{pname}\n({m[pname].units})\n" for pname in param_names]
param_labels = vl.get_free_param_labels(mv)
fig = corner.corner(
    samples_v,
    labels=param_labels,
    label_kwargs={"fontsize": 11},
    range=[0.999] * ndim,
    truths=maxlike_params_v / scale_factors,
    plot_datapoints=False,
    hist_kwargs={"density": True},
)

# for ii, prior in enumerate(mv.priors):
#     a = min(samples_v[:, ii])
#     b = max(samples_v[:, ii])
#     xs = np.linspace(a, b, 100)
#     distr = vl.distr(prior, mv.param_handler._default_params_tuple)
#     ys = np.array([vl.logpdf(distr, x) for x in xs])

#     plt_num = 1 + ii * (ndim + 1)
#     plt.subplot(ndim, ndim, plt_num)
#     plt.plot(xs, 100 * ys)

# if os.path.isfile(f"{m.PSR.value}_chain_emcee_jl.txt"):
#     samples_v1 = np.genfromtxt(f"{m.PSR.value}_chain_emcee_jl.txt", skip_header=1)
#     corner.corner(
#         samples_v1,
#         range=[0.999999] * ndim,
#         plot_datapoints=False,
#         fig=fig,
#         color="blue",
#         hist_kwargs={"density": True},
#     )

plt.suptitle(m.PSR.value)
plt.tight_layout()
plt.show()
