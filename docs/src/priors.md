# Prior and Posterior distributions

The prior distributions of each free parameter is represented as a `Prior` object.  
These use the `Distribution`s defined by [`Distributions.jl`](https://juliastats.org/Distributions.jl/) 
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
    }
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
at a lower bound 0.