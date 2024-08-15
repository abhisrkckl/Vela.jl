export calc_tzr_phase, form_residual, form_residuals

"""Compute the timing residual corresponding to a single TOA."""
function form_residual(model::TimingModel, toa::TOA, params::NamedTuple, tzrphase::GQ)::GQ
    ctoa = correct_toa(model, toa, params)
    dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
    return dphase / doppler_shifted_spin_frequency(ctoa)
end

"""Compute the timing residuals corresponding to a collection of TOAs.

Roughly corresponds to `Residuals.calc_time_resids` in `PINT`."""
function form_residuals(
    model::TimingModel,
    toas::Vector{TOA},
    params::NamedTuple,
)::Vector{GQ}
    tzrphase = calc_tzr_phase(model, params)
    return [form_residual(model, toa, params, tzrphase) for toa in toas]
end

"""Compute the phase corresponding to the TZR TOA."""
function calc_tzr_phase(model::TimingModel, params::NamedTuple)
    ctzrtoa = correct_toa(model, model.tzr_toa, params)
    return ctzrtoa.phase
end
