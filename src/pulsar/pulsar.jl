export Pulsar

struct Pulsar{M<:TimingModel,T<:TOABase}
    model::M
    toas::Vector{T}
end

calc_lnlike(psr::Pulsar, params) = calc_lnlike(psr.model, psr.toas, params)
calc_lnlike_serial(psr::Pulsar, params) = calc_lnlike_serial(psr.model, psr.toas, params)

calc_lnprior(psr::Pulsar, params) = calc_lnprior(psr.model, params)
prior_transform(psr::Pulsar, cube) = prior_transform(psr.model, cube)

calc_lnpost(psr::Pulsar, params) = calc_lnpost(psr.model, psr.toas, params)
calc_lnpost_serial(psr::Pulsar, params) = calc_lnpost_serial(psr.model, psr.toas, params)
calc_lnpost_vectorized(psr::Pulsar, paramss) =
    calc_lnpost_vectorized(psr.model, psr.toas, paramss)
