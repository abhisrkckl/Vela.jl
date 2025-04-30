export Prior,
    lnprior,
    calc_lnprior,
    prior_transform,
    get_lnprior_func,
    get_prior_transform_func,
    distr

"""Abstract base class of all priors."""
abstract type Prior end

"""Evaluate the log-prior."""
lnprior(prior::Prior, params::NamedTuple) =
    logpdf(distr(prior, params), param_value(prior, params))
lnprior(priors, params::NamedTuple) = sum(prior -> lnprior(prior, params), priors)
lnprior(priors, params) = mapreduce(lnprior, +, priors, params)

calc_lnprior(model::TimingModel, params::NamedTuple) = lnprior(model.priors, params)
calc_lnprior(model::TimingModel, params) = lnprior(model.priors, params)

"""
    get_lnprior_func(::TimingModel)::Function

Returns a callable that evaluates the log-prior given a collection of parameter values.
"""
get_lnprior_func(model::TimingModel) = params -> calc_lnprior(model, params)

"""Evaluate the prior transform function."""
prior_transform(priors, cube) = map(prior_transform, priors, cube)
prior_transform(model::TimingModel, cube) = prior_transform(model.priors, cube)

"""
    get_prior_transform_func(::TimingModel)::Function

Returns a callable that evaluates the prior transform given a point within a unit hypercube.
"""
get_prior_transform_func(model::TimingModel) = cube -> prior_transform(model, cube)

distr_args(d::Distribution) = ntuple(i -> getfield(d, i), fieldcount(typeof(d)))
distr_args(d::Truncated) = distr_args(d.untruncated)
distr_args(p::Prior) = distr_args(p.distribution)

distr_name(d::Distribution) = typeof(d).name.name
distr_name(d::Truncated) = distr_name(d.untruncated)
distr_name(p::Prior) = distr_name(p.distribution)
