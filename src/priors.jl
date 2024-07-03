function _lnprior(priors_PAR::Tuple, params::NamedTuple, _PAR::Symbol)
    vals = map(value, params[_PAR])
    dists = priors_PAR[1:length(vals)]
    return sum(map(logpdf, dists, vals))
end
