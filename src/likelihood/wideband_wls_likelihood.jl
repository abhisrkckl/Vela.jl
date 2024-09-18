function _wls_lnlike_term(
    model::TimingModel,
    wtoa::WidebandTOA,
    params::NamedTuple,
    tzrphase,
)
    cwtoa = correct_toa(model, wtoa, params)

    dphase = GQ{Float64}(phase_residual(wtoa.toa, cwtoa.toa_correction) - tzrphase)
    tres = dphase / doppler_shifted_spin_frequency(cwtoa.toa_correction)
    err2 = scaled_toa_error_sqr(wtoa.toa, cwtoa.toa_correction)
    norm = log(value(err2))
    lnlike_term_toa = value(tres * tres / err2) + norm

    dmres = dm_residual(wtoa.dminfo, cwtoa.dm_correction)
    dmerr2 = scaled_dm_error_sqr(wtoa.dminfo, cwtoa.dm_correction)
    norm = log(value(dmerr2))
    lnlike_term_dm = value(dmres * dmres / dmerr2) + norm

    return lnlike_term_toa + lnlike_term_dm
end
