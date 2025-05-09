# Prior and Posterior distributions

The prior distributions of each free parameter is represented as a `Vela.Prior` object.  
These use the `Distribution`s defined by the [`Distributions.jl`](https://juliastats.org/Distributions.jl/) package
under the hood. 
```@docs
Prior
```

The priors corresponding to `Parameter`s and `MultiParameter`s are represented by two subtypes
of `Prior`.
```@eval
using InteractiveUtils
using AbstractTrees
using Vela
using Markdown

AbstractTrees.children(d::DataType) = subtypes(d)
Markdown.MD(Markdown.Code(repr_tree(Prior)))
```

```@docs
SimplePrior
SimplePriorMulti
```

Each `Prior` has `lnprior` and `prior_transform` methods which compute  the log-prior distribution 
and the [prior transform function](http://kylebarbary.com/nestle/prior.html) respectively
for that parameter. The former is necessary for MCMC samplers and the latter for nested samplers.
Please note that these functions act on `Float64`s rather than `GQ`s because the samplers 
only provide `Float64`s. A `TimingModel` also has the `lnprior` and `prior_transform` methods; they evaluate 
the joint log-prior and the prior transform over all free parameters. The `get_lnprior_func` and `get_prior_transform_func`
functions return callables that can be passed on to samplers.
```@docs
lnprior
prior_transform
get_lnprior_func
get_prior_transform_func
```

The log-posterior is the sum of the log-likelihood and the log-prior up to an additive 
constant. The `get_lnpost_func` function returns a callable that evaluates the log-posterior 
can be passed on to samplers. Note that the expensive log-likelihood is evaluated only if the 
log-prior is finite. 
```@docs
get_lnpost_func
```

## Priors for different parameters

In principle, the prior for each parameter has to be set based on our prior knowledge. Indeed,
we may have prior information on some of the parameters from previous timing experiments, VLBI 
campaigns, detection of counterparts in other parts of the electromagnetic spectrum (e.g., using GAIA),
etc. Or priors may be estimated from population statistics using something like 
[`psrcat`](https://www.atnf.csiro.au/research/pulsar/psrcat/).

However, for many parameters, pulsar timing provides so much signal-to-noise ratio that the effect of the prior
on the posterior distrubution is entirely negligible. This is the case for parameters like `F0`, `F1`, `RAJ`, `DECJ`,
etc. In such a case, it may be enough to use "cheat" priors that are based on the frequentist uncertainties for
these parameters (given in the `par` file). Specifically, we use uniform distributions centered around the frequentist
estimate whose width is several times (e.g., 10x) the frequentist uncertainty. Care must be taken to ensure that the
data provides enough S/N for the parameter for the "cheat" prior to be valid, otherwise we will be effectively 
[double dipping](https://en.wikipedia.org/wiki/Circular_analysis).

For some parameters, e.g., the orbital inclination, we have physically motivated default prior distributions.
```@docs
KINPriorDistribution
SINIPriorDistribution
STIGMAPriorDistribution
SHAPMAXPriorDistribution
```

The preference for the prior distributions is user defined distribution > default distribution > "cheat" distribution.
**Please take care to ensure that the wrong parameter doesn't end up with a "cheat" distribution.**

See the documentation for `Distrubutions.jl` to see what distributions are available.

## Specifying priors in `pyvela`

While creating a `pyvela.SPNTA` object, the priors are specified using the `cheat_prior_scale` and `custom_priors`
arguments. 

The `cheat_prior_scale` argument defines the scale factor by which the frequentist uncertainties are multiplied to 
obtain the "cheat" prior widths. 

The `custom_priors` argument contains the user-defined prior distributions as a dictionary. It can be a Python `dict` 
or a filename (`str`) / `IO` object containing a `JSON` representation of the dictionary (see below). It supports both parameter 
names and prefixes as dict keys.  For example, if an entry for "EFAC" is present, it will set the prior for all EFAC 
parameters. If "EFAC1" is present, it will set the prior for EFAC1 specifically. If both "EFAC" and "EFAC1" are present, 
the latter sets the prior for EFAC1, whereas the former sets the priors for all other EFACs. 

## Representing priors in a `JSON` file
Prior distributions available in `Distributions.jl` can be represented as a `JSON` file like so:
```
{
    "EFAC": {
        "distribution": "LogNormal",
        "args": [0.0, 0.5]
    },
    "EQUAD": {
        "distribution": "LogUniform",
        "args": [1e-2, 2.0]
    },
    "M2": {
        "distribution": "Normal",
        "args": [0.1, 0.02],
        "lower": 0.0
    }
}
```
The `distribution` attribute for each parameter corresponds to a `UnivariateDistribution` available in `Distributions.jl`
(see [here](https://juliastats.org/Distributions.jl/stable/univariate/)). `args` are arguments to the 
`UnivariateDistribution` type's constructor. The `lower` and `upper` attributes represent the lower and 
upper bounds for truncating the distribution (see [`truncated`](https://juliastats.org/Distributions.jl/stable/truncate/)).

Note that the values above should be given in their "normal" units as they appear in the par files. 
Specifically, the prior on M2 corresponds  to 0.1 ± 0.02 Msun, where the normal distribution is truncated 
at a lower bound of 0. 

## Default priors

The following table lists the parameters for which default prior distributions are available.
More parameters will be added to this table in the future.

| **Parameter** | **Description & Unit**                                                      | **Prior distribution**                                          | **Remarks**      |
|:-------------:|:----------------------------------------------------------------------------|:----------------------------------------------------------------|:-----------------|
| DMEFAC        | Wideband DM uncertainty correction (multiplicative)                         | `LogNormal[0, 0.25]`                                            |                  |
| DMEQUAD       | Wideband DM uncertainty correction (added in quadrature) (pc/cm^3)          | `LogUniform[1e-8, 1e-2]`                                        |                  |
| DMX_          | Piecewise-constant DM parameters (pc / cm^3)                                | `Normal[0, 1e-2]`                                               |                  |
| EFAC          | TOA uncertainty correction (multiplicative)                                 | `LogNormal[0, 0.25]`                                            |                  |
| EQUAD         | TOA uncertainty correction (added in quadrature) (μs)                       | `LogUniform[1e-3, 1e+2]`                                        |                  |
| ECORR         | Correlation between narrowband TOAs measured from the same observation (μs) | `LogUniform[1e-3, 1e+2]`                                        |                  |
| H3            | Orthometric Shapiro delay range (s)                                         | `LogNormal[1e-9*G*Msun/c^3, 100*G*Msun/c^3]`                    |                  |
| JUMP          | Instrumental time jump                                                      | `Uniform[-1/(2*F0), 1/(2*F0)]`                                  |                  |
| KIN           | Orbital inclination (deg)                                                   | `P(KIN) = sin(KIN)`                                             | 0 <= KIN <= 90.  |
| KOM           | Longitude of ascending node (deg)                                           | `LogNormal[0, 360]`                                             |                  |
| M2            | Companion mass (MSun)                                                       | `LogUniform[1e-9, 100]`                                         |                  |
| PHOFF         | Overall phase offset                                                        | `Uniform[-0.5, 0.5]`                                            |                  |
| PX            | Parallax (mas)                                                              | `P(PX) = 3 / Rmax^3 / PX^4`                                     | 0 <= PX <= Rmax  |
| SHAPMAX       | Log-scaled inclination parameter                                            | `P(SHAPMAX) = (1 - exp(-SHAPMAX)) / (sqrt(2 * exp(SHAPMAX) - 1))`| 0 <= SHAPMAX    |
| SINI          | Orbital sine-inclination                                                    | `P(SINI) = SINI / sqrt(1 - SINI^2)`                             | 0 <= SINI <= 1   |
| STIGMA        | Orthometric Shapiro delay shape                                             | `P(STIGMA) = 4 * STIGMA / (1 + STIGMA^2)^2`                     | 0 <= STIGMA <= 1 |
| T0            | Epoch of periapsis passage                                                  | `Uniform[T0_1 - PB/2, T0_1 + PB/2]`                             | T0_1 is the default value of T0. |
| TASC          | Epoch of acending node passage                                              | `Uniform[TASC_1 - PB/2, TASC_1 + PB/2]`                         | TASC_1 is the default value of TASC. |
| TNCHROMAMP    | Chromatic noise log-amplitude                                               | `Uniform[-18, -9]`                                              |                  |
| TNCHROMGAM    | Chromatic noise log-amplitude                                               | `Uniform[0, 7]`                                                 |                  |
| TNDMAMP       | DM noise log-amplitude                                                      | `Uniform[-18, -9]`                                              |                  |
| TNDMGAM       | DM noise spectral index                                                     | `Uniform[0, 7]`                                                 |                  |
| TNREDAMP      | Spin noise log-amplitude                                                    | `Uniform[-18, -9]`                                              |                  |
| TNREDGAM      | Spin noise spectral index                                                   | `Uniform[0, 7]`                                                 |                  |

