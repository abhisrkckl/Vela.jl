#!/usr/bin/env python

# %%
import time
from argparse import ArgumentParser

import corner
import nestle
import numpy as np
import wquantiles as wq
from matplotlib import pyplot as plt

from pyvela import SPNTA
from pyvela import Vela as vl

np.product = np.prod


parser = ArgumentParser()
parser.add_argument("par_file")
parser.add_argument("tim_file")
parser.add_argument("-P", "--priors", default={})
args = parser.parse_args()


spnta = SPNTA(
    args.par_file, args.tim_file, cheat_prior_scale=10, custom_priors=args.priors
)


# shifts = [
#     (spnta.model.param_handler._default_params_tuple.F_.x if pname == "F0" else 0)
#     for pname in spnta.param_names
# ]


begin = time.time()
result_nestle_v = nestle.sample(
    spnta.lnlike,
    spnta.prior_transform,
    spnta.ndim,
    method="multi",
    npoints=500,
    dlogz=0.001,
    callback=nestle.print_progress,
)
end = time.time()
print(f"\nTime elapsed = {end-begin} s")


samples_v = spnta.rescale_samples(result_nestle_v.samples)


means, cov = nestle.mean_and_cov(result_nestle_v.samples, result_nestle_v.weights)
means /= spnta.scale_factors
stds = np.sqrt(np.diag(cov)) / spnta.scale_factors
for pname, mean, std in zip(spnta.param_names, means, stds):
    if pname == "F0":
        F0_ = np.longdouble(spnta.model.param_handler._default_params_tuple.F_.x)
        print(f"{pname}\t\t{(F0_ + mean):.18f}\t\t{std}")
    else:
        print(f"{pname}\t\t{mean}\t\t{std}")


corner.corner(
    samples_v,
    weights=result_nestle_v.weights,
    labels=spnta.param_names,
    label_kwargs={"fontsize": 15},
    range=[0.999999] * spnta.ndim,
    color="red",
    hist_kwargs={"density": True},
)
plt.suptitle(spnta.model.pulsar_name)
plt.tight_layout()
plt.show()


params_median = vl.read_params(
    spnta.model, wq.median(result_nestle_v.samples.T, weights=result_nestle_v.weights)
)
rv = (
    [
        vl.value(wr[0])
        for wr in vl.form_residuals(spnta.model, spnta.toas, params_median)
    ]
    if spnta.is_wideband()
    else list(map(vl.value, vl.form_residuals(spnta.model, spnta.toas, params_median)))
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
