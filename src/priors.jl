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
    pars = params[mp.param]
    @assert length(pars) <= length(mp.distributions) 

    result = 0.0
    for (par, dist) in zip(pars, mp.distributions[:length(pars)])
        result += logpdf(dist, value(par))
    end

    return sum(
        (logpdf(dist, value(par)) for (par, dist) in zip(pars, mp.distributions[:length(pars)]))
    )
end

logpdf(priors::Tuple, params::NamedTuple) = sum((logpdf(prior, params) for prior in priors))
