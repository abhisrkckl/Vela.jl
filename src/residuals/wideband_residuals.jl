"""Compute the timing residual corresponding to a single TOA."""
function form_residual(
    model::TimingModel,
    wtoa::WidebandTOA,
    params::NamedTuple,
    tzrphase::GQ,
)
    cwtoa = correct_toa(model, wtoa, params)
    dphase = GQ{Float64}(phase_residual(cwtoa.corrected_toa) - tzrphase)
    tres = dphase / doppler_shifted_spin_frequency(cwtoa.corrected_toa)
    dmres = dm_residual(cwtoa.corrected_dminfo)
    return tres, dmres
end

function form_residuals(model::TimingModel, wtoas::Vector{WidebandTOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)
    return map(wtoa -> form_residual(model, wtoa, params, tzrphase), wtoas)
end
