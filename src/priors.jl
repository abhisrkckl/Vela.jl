export calc_lnprior

function _lnprior(priors_PAR::Tuple, params::NamedTuple, _PAR::Symbol)
    vals = map(value, params[_PAR])
    dists = priors_PAR[1:length(vals)]
    return sum(map(logpdf, dists, vals))
end

@unroll function lnprior(components::Tuple, params::NamedTuple)
    result = 0.0
    @unroll for component in components
        result += lnprior(component, params)
    end
    return result
end

calc_lnprior(model::TimingModel, params::NamedTuple) = lnprior(model.components, params)
