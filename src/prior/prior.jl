abstract type Prior end

lnprior(prior::Prior, params::NamedTuple) = logpdf(distr(prior, params), value(params[prior.name]))
prior_transform(prior::Prior, q) = GQ(quantile(distr(prior, params), q), prior.dim)

lnprior(priors, params::NamedTuple) = sum(prior -> lnprior(prior, params), priors)
prior_transform(priors, cube) = map(prior_transform, priors, cube)

lnprior(model::TimingModel, params::NamedTuple) = lnprior(model.priors, params)
prior_transform(model::TimingModel, cube) = prior_transform(model.priors, cube)

get_lnprior_func(model::TimingModel) = params -> lnprior(model, params)
get_prior_transform_func(model::TimingModel) = cube -> prior_transform(model, cube)

struct SimplePrior{D<:Distribution} <: Prior
    name::Symbol
    dim::Int
    distribution::D
end

distr(sp::SimplePrior, ::NamedTuple)::Distribution = sp.distribution
