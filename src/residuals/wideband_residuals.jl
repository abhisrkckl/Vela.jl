function form_residual(
    model::TimingModel,
    wtoa::WidebandTOA,
    params::NamedTuple,
    tzrphase::GQ,
)
    cwtoa = correct_toa(model, wtoa, params)
    dphase = GQ{Float64}(phase_residual(wtoa.toa, cwtoa.toa_correction) - tzrphase)
    tres = dphase / doppler_shifted_spin_frequency(cwtoa.toa_correction)
    dmres = dm_residual(wtoa.dminfo, cwtoa.dm_correction)
    return tres, dmres
end

function form_residuals(model::TimingModel, wtoas::Vector{WidebandTOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)
    return map(wtoa -> form_residual(model, wtoa, params, tzrphase), wtoas)
end
