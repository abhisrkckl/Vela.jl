"""Compute the timing residual corresponding to a single TOA."""
function form_wideband_residual(
    model::TimingModel,
    wtoa::WidebandTOA,
    params::NamedTuple,
    tzrphase::GQ,
)
    cwtoa = correct_toa(model, wtoa, params)
    dphase = GQ{Float64}(phase_residual(cwtoa.corrected_toa) - tzrphase)
    time_residual = dphase / doppler_shifted_spin_frequency(cwtoa.corrected_toa)
    dm_residual = dm_residual(cwtoa.corrected_dminfo)
    return time_residual, dm_residual
end

function form_wideband_residuals(
    model::TimingModel,
    wtoas::Vector{WidebandTOA},
    params::NamedTuple,
)
    tzrphase = calc_tzr_phase(model, params)
    return map(toa -> form_wideband_residual(model, wtoa, params, tzrphase), wtoas)
end
