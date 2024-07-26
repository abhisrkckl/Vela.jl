export Prior,
    SimplePrior, lnprior, prior_transform, get_lnprior_func, get_prior_transform_func, distr

"""Abstract base class of all priors."""
abstract type Prior end

"""Extract the value of a parameter for evaluating its prior."""
param_value(prior::Prior, params::NamedTuple) = param_value(prior, params[prior.name])
param_value(::Prior, param::GQ) = value(param)
param_value(prior::Prior, params::NTuple) = value(params[prior.index])

"""Evaluate the log-prior."""
lnprior(prior::Prior, params::NamedTuple) =
    logpdf(distr(prior, params), param_value(prior, params))
lnprior(priors, params::NamedTuple) = sum(prior -> lnprior(prior, params), priors)
lnprior(model::TimingModel, params::NamedTuple) = lnprior(model.priors, params)
get_lnprior_func(model::TimingModel) = params -> lnprior(model, params)

"""Evaluate the prior transform function."""
prior_transform(priors, cube) = map(prior_transform, priors, cube)
prior_transform(model::TimingModel, cube) = prior_transform(model.priors, cube)
get_prior_transform_func(model::TimingModel) = cube -> prior_transform(model, cube)

struct SimplePrior{D<:Distribution} <: Prior
    name::Symbol
    index::UInt
    distribution::D
end

SimplePrior(name::Symbol, distr::Distribution) = SimplePrior(name, UInt(0), distr)

distr(sp::SimplePrior, ::NamedTuple)::Distribution = sp.distribution

prior_transform(prior::SimplePrior, q) = quantile(prior.distribution, q)
