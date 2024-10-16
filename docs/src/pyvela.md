# The `pyvela` interface

`Vela.jl` must interact with Python for three reasons: (1) most pulsar astronomers are more 
familiar with Python than Julia (2) Python has many samplers that have no counterpart in 
Julia, and most importantly, (3) it is a real pain to implement `par` and `tim` file readers.

The `pyvela` interface is demonstrated below using an example with the `emcee` sampler.

## Reading `par` and `tim` files
While the `par` file format is outwardly simple, it is rife with special cases and arbitrary units,
and reconstructing a timing & noise model from a `par` file is pretty hard. The `tim` files are
kind of easier to handle, but there are multiple `tim` file formats and some datasets can have
multiple format TOAs in the same file. So I am not re-inventing the wheel and instead using
`PINT` to do this. This also has the advantage of not having to implement clock corrections and
solar system ephemeris calculations.

Here is how we read read a pair of `par` and `tim` files in `pyvela`:
```
from pyvela import read_model_and_toas, Vela as vl
from juliacall import Main as jl
import numpy as np

jl.seval("using Distributions")

parfile, timfile = "NGC6440E.par", "NGC6440E.tim"
mv, tv = read_model_and_toas(
    parfile, 
    timfile,
    cheat_prior_scale=5,
    custom_prior_dicts={
        "PHOFF": jl.Uniform(-0.5, 0.5)
    }
)
```
Here, `mv` is a `Vela.TimingModel` object and `tv` is a `Vector{Vela.TOA}` object. The `read_model_and_toas`
function reads the `par` and `tim` files using `pint.models.get_model_and_toas` function under the hood and
converts the resulting `pint.models.TimingModel` and `pint.toa.TOAs` objects. Note that this conversion does 
not conserve all the information, and the reverse is not possible. 

`cheat_prior_scale` defines the scale factor by which the frequentist uncertainties are multiplied to obtain the "cheat" prior widths. 
The `custom_prior_dicts` argument contains the user-defined prior distributions, which should be instances of the 
`Distributions.UnivariateDistribution` type. It supports both parameter names and prefixes as dict keys. For example,
if an entry for "EFAC" is present, it will set the prior for all EFAC parameters. If "EFAC1" is present, it will 
set the prior for EFAC1 specifically. If both "EFAC" and "EFAC1" are present, the latter sets the prior for EFAC1, whereas
the former sets the priors for all other EFACs. 

In the above example, the priors for all parameters except "PHOFF" are set using "cheat" priors.

Everything defined in `Vela.jl` will be available through the `vl` namespace above, e.g., `vl.TimingModel` and `vl.TOA`.
Things that were explicitly imported into `Vela.jl` are also available, e.g., `vl.GQ` is the `GQ` type from `GeometricUnits.jl`.
Other things from Julia are accessed using the `juliacall` package.

## Getting the log-posterior function etc.
```
lnpost = vl.get_lnpost_func(mv, tv, True)
prior_transform = vl.get_prior_transform_func(mv)
```
The last argument is `vectorize`, which tells the function to output a callable that executes 
the log-posterior for multiple data points at the same time. The prior transform is useful here
to easily draw samples from the prior distribution even if we are not using a nested sampler. 

`vectorize` improves periformance only when multiple threads are available to the process. The number of threads 
can be controlled using the `PYTHON_JULIACALL_THREADS` environment variable. Don't set it to a value greater
than the number of available CPU cores.

```
param_names = vl.get_free_param_names(mv)
scale_factors = vl.get_scale_factors(mv)
```
`param_names` contains the names of the free parameters in the model. `scale_factors` contain the
scale factors that convert the different parameters from `Vela.jl` internal representation to their
usual units. 

```
maxlike_params = np.array(
    [
        vl.read_param_values_to_vector(mv)
    ]
)
```
The `maxlike_params` variable contains the free parameter values read from the par file.
Note that we are enclosing this within a 1 x Ndim matrix. This is because we set `vectorize` 
to true, so that `lnpost` expects a collection of points in the parameter space. 

Now, let us see if `lnpost` actually works.
```
print(lnpost(maxlike_params))
```

## Setting up the sampler
```
import emcee

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
```
`p0` contains a number of random draws from the prior distribution. Note the `vectorize=True`
while creating the sampler. This should match what we passed into `get_lnpost_func`.

## Running MCMC
```
sampler.run_mcmc(p0, 6000, progress=True, progress_kwargs={"mininterval": 1})
```
This will only take a few seconds because the dataset is very small. Larger datasets may take
minutes or hours depending on the size and the computing power available.

```
samples_v_0 = sampler.get_chain(flat=True, discard=1000, thin=40)
```
This flattens the multiple chains `emcee` was running. Also note the burn-in (`discard`) 
and the thinning.

## Printing out the results
```
shifts = [(mv.param_handler._default_params_tuple.F_.x if pname == "F0" else 0) for pname in param_names]
means = (np.mean(samples_v_0, axis=0) + shifts) / scale_factors
stds = np.std(samples_v, axis=0)
for idx, (pname, mean, std) in enumerate(zip(param_names, means, stds)):
    print(f"{pname}\t\t{mean}\t\t{std}")
```
Recall that the F0 value is stored as a sum of two `Float64`s. `shifts` is used to put these two parts 
back together for printing. 

## Plotting
```
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
```

The output looks something like this:
![NGC6440E_posterior](NGC6440E_posterior.png)
The blue lines represent the maximum-likelihood estimates.
