export calc_lnpost

function calc_lnpost(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    lnpr = calc_lnprior(model, params)
    return isfinite(lnpr) ? (lnpr + calc_lnlike(model, toas, params)) : -Inf
end
