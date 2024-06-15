using Distributions

abstract type Prior end

struct SimplePrior{Distr<:Distribution} <: Prior
    param::Symbol
    distribution::Distr
end

struct MultiPrior{DistrTuple<:Tuple} <: Prior
    param::Symbol
    distributions::DistrTuple
end

logpdf(sp::SimplePrior, params::NamedTuple) = logpdf(sp.distribution, value(params[sp.param]))

function logpdf(mp::MultiPrior, params::NamedTuple)
    vals = value.(params[mp.param])
    @assert length(vals) <= length(mp.distributions) 

end

logpdf(priors::Tuple, params::NamedTuple) = sum((logpdf(prior, params) for prior in priors))
