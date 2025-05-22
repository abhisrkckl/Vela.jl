# Getting started: The `pyvela` interface

`Vela.jl` must interact with Python for three reasons: (1) most pulsar astronomers are more 
familiar with Python than Julia (2) Python has many samplers that have no counterpart in 
Julia, and most importantly, (3) it is a real pain to implement `par` and `tim` file readers.
The `pyvela` interface allows one to access `Vela.jl` from Python. This is the simplest way
to get started with `Vela.jl`.

The `pyvela` interface is demonstrated below using an example with the `emcee` sampler.

## Reading `par` and `tim` files using the `SPNTA` class

A pulsar timing dataset usually contains a `par` file and a `tim` file. The `tim` file contains the
pulse time of arrival (TOA) measurements and related metadata, and the `par` file contains a timing &
noise model that fits the TOAs along with its parameter values. `pyvela` reads these files with the 
help of the [`PINT`](https://nanograv-pint.readthedocs.io/) package. `PINT` also does the clock 
corrections and solar system ephemeris calculations. See [this page](https://nanograv-pint.readthedocs.io/en/latest/explanation.html)
for detailed explanations on these topics.

Here is how we read read a pair of `par` and `tim` files in `pyvela`:
```
from pyvela import SPNTA, Vela, pyvela_plot
import numpy as np

parfile, timfile = "NGC6440E.par", "NGC6440E.tim"
spnta = SPNTA(
    parfile, 
    timfile,
    cheat_prior_scale=100,
    custom_priors={},
    analytic_marginalized_params=["PHOFF", "F"],
)
```
Here, the `cheat_prior_scale` and `custom_priors` arguments are used for specifying the prior distributions for
model parameters. See [Representing priors in a `JSON` file](@ref) for more details.
`analytic_marginalized_params` represents parameters that are analytically marginalized assuming infinitely 
wide improper priors. Only parameters whose effect on the timing residuals is approximately linear are allowed
to be analytically marginalized. This is different from `ENTERPRISE`, where the timing model is linearized in 
all timing parameters.

An `SPNTA` object stores a pulsar timing & noise model along with the measured TOAs.
Here, `spnta.model` is a `Vela.TimingModel` object and `spnta.toas` is a `Vector{Vela.TOA}` object. 
The `SPNTA` class reads the `par` and `tim` files using the [`pint.models.get_model_and_toas()`](https://nanograv-pint.readthedocs.io/en/latest/_autosummary/pint.models.model_builder.get_model_and_toas.html)
function under the hood and converts the resulting `pint.models.TimingModel` and `pint.toa.TOAs` objects
into the corresponding `Vela.jl` objects. Note that this conversion does not preserve all the information
present in the `PINT` objects, and is not reversible.

The `SPNTA` class internally uses the `Vela.Pulsar` type to store the timing model and the TOAs together.
Currently it is assumed that all the TOAs are of the same paradigm (narrowband or wideband). Inhomogeneous
datasets are not yet supported.
```@docs
Pulsar
```

Everything defined in `Vela.jl` will be available through the `Vela` namespace above, e.g., `Vela.TimingModel` 
and `Vela.TOA`. Things that were explicitly imported into `Vela.jl` are also available, e.g., `Vela.GQ` is the `GQ` type 
from `GeometricUnits.jl` (see [Quantities](@ref)). Other things from Julia are accessed using the 
[`juliacall`](https://juliapy.github.io/PythonCall.jl/stable/juliacall/) package.

The `SPNTA` object created above provides functions like `lnlike()`, `lnprior()`, `prior_transform()`, `lnpost()`, and 
`lnpost_vectorized()`. These provide the log-likelihood, log-prior, prior transform, and log-posterior 
functions which call `Vela.jl` under the hood. The difference between `spnta.lnpost` and 
`spnta.lnpost_vectorized` is that if multiple threads are allowed (by setting the `PYTHON_JULIACALL_THREADS`
environment variable), `spnta.lnpost` parallelizes a single log-posterior computation across TOAs, whereas
`spnta.lnpost_vectorized` computes the log-posterior at multiple points in the parameter space parallelly.
The latter can be used with samplers such as `emcee` and `zeus`. Make sure not to set `PYTHON_JULIACALL_THREADS` 
to a value greater than the number of available CPU cores.

The `SPNTA` object also provides several useful attributes. These are:
1. `ndim` - Number of free parameters.
2. `ntmdim` - Number of timing model parameters (excludes noise parameters and hyperparameters).
3. `param_names` - Free paramater names in `PINT` convention.
4. `param_units` - Free paramater unit strings in `astropy.units` convention as used in `PINT`.
5. `param_labels` - Free parameter labels containing names and units.
6. `param_prefixes` - Free parameter prefixes in `PINT` convention.
7. `marginalized_param_names` - Names of the analytically marginalized parameters. This may include correlated noise amplitudes as well as linear timing model paramaters as defined by the `analytic_marginalized_params` option. This follows `PINT` conventions as much as possible.
8. `scale_factors` - Scale factors used to convert parameters from normal pulsar timing units to `Vela.jl`'s internal units.
9. `default_params` - Default free parameter values taken from the `par` file (in `Vela.jl`'s internal units).
10. `maxpost_params` - Maximum-posterior parameters estimated using the Nelder-Mead method (in `Vela.jl`'s internal units).
11. `has_marginalized_gp_noise` - Whether the timing & noise model contains marginalized Gaussian process components (Boolean).
12. `has_ecorr_noise` - Whether the timing & noise model contains ECORR noise (Boolean).
13. `wideband` - Whether the TOAs are wideband (Boolean).
14. `mjds` - Observing epochs in MJD.

The methods available in `SPNTA` other than the ones mentioned above are:
1. `rescale_samples(samples_raw)` - Rescales samples from `Vela.jl` internal units to the usual units used in pulsar astronomy (see [this page](https://nanograv-pint.readthedocs.io/en/latest/timingmodels.html#supported-parameters)).
2. `get_marginalized_gp_noise_realization(params)` - Given free parameter values, compute the Gaussian process realization of the analytically marginalized parameters.
3. `time_residuals(params)` - Compute the time residuals given a set of paramater values.
4. `whitened_time_residuals(params)` - Compute the whitened time residuals given a set of paramater values.
5. `dm_residuals(params)` - Compute the DM residuals given a set of paramater values (wideband only).
6. `whitened_dm_residuals(params)` - Compute the whitened DM residuals given a set of paramater values (wideband only).
7. `scaled_toa_unceritainties(params)` - Compute the scaled TOA uncertainties given a set of paramater values.
8. `scaled_dm_unceritainties(params)` - Compute the scaled DM uncertainties given a set of paramater values (wideband only).
9. `model_dm(params)` - Compute the DM predicted by the model given a set of paramater values.

## Setting up the sampler

We show below an example using the [`emcee`](https://emcee.readthedocs.io/) sampler.

```
import emcee

nwalkers = spnta.ndim * 5
p0 = np.array([spnta.prior_transform(cube) for cube in np.random.rand(nwalkers, spnta.ndim)])

sampler = emcee.EnsembleSampler(
    nwalkers,
    spnta.ndim,
    spnta.lnpost_vectorized,
    moves=[emcee.moves.StretchMove(), emcee.moves.DESnookerMove()],
    vectorize=True,
)
```
`p0` contains random draws from the prior distribution (see [this page](https://en.wikipedia.org/wiki/Inverse_transform_sampling)
for an explanation on how this works). Note the `vectorize=True` while creating the sampler. This must be given if we are using the 
vectorized log-posterior. In practice, the `moves` should be optimized based on the problem at hand.
 
## Running MCMC and getting the samples
```
sampler.run_mcmc(p0, 6000, progress=True)
```
This will only take a few seconds because the dataset is very small. Larger datasets may take
minutes or hours depending on the size and the computing power available.

```
samples_raw = sampler.get_chain(flat=True, discard=1000, thin=50)
```
This flattens the multiple chains `emcee` was running. Also note the burn-in (`discard`) 
and the thinning. These samples are in `Vela.jl`'s internal units (see [Quantities](@ref)).

```
samples = spnta.rescale_samples(samples_raw)
```
This converts the samples into the usual units used in pulsar astronomy.

## Printing out the results
```
means = np.mean(samples_v, axis=0)
stds = np.std(samples_v, axis=0)
for idx, (pname, mean, std) in enumerate(zip(spnta.param_names, means, stds)):
    if pname == "F0":
        F0_ = np.longdouble(spnta.model.param_handler._default_params_tuple.F_.x)
        print(f"{pname}\t\t{mean + F0_}\t\t{std}")    
    else:
        print(f"{pname}\t\t{mean}\t\t{std}")
```
The special treatment for F0 is explained in [Precision](@ref).

## Plotting
```
import corner
import matplotlib.pyplot as plt

fig = corner.corner(
    samples_v,
    labels=spnta.param_labels,
    label_kwargs={"fontsize": 11},
    range=[0.999] * spnta.ndim,
    plot_datapoints=False,
    hist_kwargs={"density": True},
    labelpad=0.3,
)
plt.suptitle(spnta.model.pulsar_name)
plt.tight_layout()
plt.show()
```

The output looks something like this:
![NGC6440E_posterior](NGC6440E_posterior.png)

Note that the F0 plot is centered around 0. Refer to [Precision](@ref) for why this is.