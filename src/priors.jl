export lnprior

function _lnprior(priors_PAR::Tuple, params::NamedTuple, _PAR::Symbol)
    vals = map(value, params[_PAR])
    dists = priors_PAR[1:length(vals)]
    return sum(map(logpdf, dists, vals))
end

function lnprior(model::TimingModel, params::NamedTuple)
    _lnpr = component -> lnprior(component, params)
    return sum(map(_lnpr, model.components))
end
