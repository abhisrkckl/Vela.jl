export SimplePrior, SimplePriorMulti, PriorSourceType

@enum PriorSourceType DEFAULT_PRIOR CHEAT_PRIOR USER_DEFINED_PRIOR

abstract type SimplePriorBase <: Prior end

distr(sp::SimplePriorBase, ::NamedTuple)::Distribution = sp.distribution
prior_transform(prior::SimplePriorBase, q) = quantile(prior.distribution, q)

lnprior(prior::SimplePriorBase, param::Float64) = logpdf(prior.distribution, param)

"""
    SimplePrior{name,D<:Distribution}

A univariate prior for a single `Parameter`.
"""
struct SimplePrior{name,D<:Distribution} <: SimplePriorBase
    distribution::D
    source_type::PriorSourceType
end

function SimplePrior{name}(distr::Distribution, source_type::PriorSourceType) where {name}
    @assert name isa Symbol
    return SimplePrior{name,typeof(distr)}(distr, source_type)
end

param_name(::SimplePrior{name,D}) where {name,D<:Distribution} = name
param_value(::SimplePrior{name,D}, params::NamedTuple) where {name,D<:Distribution} =
    value(params[name])

"""
    SimplePriorMulti{name,index,D<:Distribution}

A univariate prior for a single parameter belonging to a `MultiParameter`.
"""
struct SimplePriorMulti{name,index,D<:Distribution} <: SimplePriorBase
    distribution::D
    source_type::PriorSourceType
end

function SimplePriorMulti{name,index}(
    distr::Distribution,
    source_type::PriorSourceType,
) where {name,index}
    @assert name isa Symbol
    @assert index isa Integer && index >= 0
    return SimplePriorMulti{name,index,typeof(distr)}(distr, source_type)
end

param_name(::SimplePriorMulti{name,index,D}) where {name,index,D<:Distribution} = name
param_index(::SimplePriorMulti{name,index,D}) where {name,index,D<:Distribution} = index
param_value(
    ::SimplePriorMulti{name,index,D},
    params::NamedTuple,
) where {name,index,D<:Distribution} = value(params[name][index])
