export calc_tzr_phase, form_residual, form_residuals

"""
    form_residual(::TimingModel, ::TOA, params::NamedTuple, tzrphase::GQ)::GQ
    
Compute the timing residual corresponding to a single narrowband TOA.
"""
function form_residual(model::TimingModel, toa::TOA, params::NamedTuple, tzrphase::GQ)::GQ
    ctoa = correct_toa(model, toa, params)
    dphase = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)
    return dphase / doppler_shifted_spin_frequency(ctoa)
end

"""
    form_residuals(::TimingModel, ::Vector{TOA}, params::NamedTuple)::Vector{GQ}

Compute the timing residuals corresponding to a collection of narrowband TOAs.
"""
function form_residuals(
    model::TimingModel,
    toas::Vector{TOA},
    params::NamedTuple,
)::Vector{GQ}
    tzrphase = calc_tzr_phase(model, params)
    return map(toa -> form_residual(model, toa, params, tzrphase), toas)
end

"""Compute the phase corresponding to the TZR TOA."""
function calc_tzr_phase(model::TimingModel, params::NamedTuple)
    ctzrtoa = correct_toa(model, model.tzr_toa, params)
    return ctzrtoa.phase
end
