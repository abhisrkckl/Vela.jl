export SimplePrior, SimplePriorMulti

abstract type SimplePriorBase <: Prior end

distr(sp::SimplePriorBase, ::NamedTuple)::Distribution = sp.distribution
prior_transform(prior::SimplePriorBase, q) = quantile(prior.distribution, q)

lnprior(prior::SimplePriorBase, param::Float64) = logpdf(prior.distribution, param)

struct SimplePrior{name,D<:Distribution} <: SimplePriorBase
    distribution::D
end

function SimplePrior{name}(distr::Distribution) where {name}
    @assert name isa Symbol
    return SimplePrior{name,typeof(distr)}(distr)
end

param_name(::SimplePrior{name,D}) where {name,D<:Distribution} = name
param_value(::SimplePrior{name,D}, params::NamedTuple) where {name,D<:Distribution} =
    value(params[name])

struct SimplePriorMulti{name,index,D<:Distribution} <: SimplePriorBase
    distribution::D
end

function SimplePriorMulti{name,index}(distr::Distribution) where {name,index}
    @assert name isa Symbol
    @assert index isa Integer && index >= 0
    return SimplePriorMulti{name,index,typeof(distr)}(distr)
end

param_name(::SimplePriorMulti{name,index,D}) where {name,index,D<:Distribution} = name
param_index(::SimplePriorMulti{name,index,D}) where {name,index,D<:Distribution} = index
param_value(
    ::SimplePriorMulti{name,index,D},
    params::NamedTuple,
) where {name,index,D<:Distribution} = value(params[name][index])