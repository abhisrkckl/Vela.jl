function form_residual(model::TimingModel, toa::TOA, params::NamedTuple, tzrphase::GQ)::GQ
    ctoa = correct_toa(model, toa, params)
    dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
    return dphase / doppler_shifted_spin_frequency(ctoa)
end

function form_residuals(
    model::TimingModel,
    toas::Vector{TOA},
    params::NamedTuple,
)::Vector{GQ}
    tzrphase = calc_tzr_phase(model, params)
    return [form_residual(model, toa, params, tzrphase) for toa in toas]
end

function calc_tzr_phase(model::TimingModel, params::NamedTuple)
    ctzrtoa = correct_toa(model, model.tzr_toa, params)
    return ctzrtoa.phase
end