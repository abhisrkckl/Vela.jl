function _wls_chi2_term(model::TimingModel, wtoa::WidebandTOA, params::NamedTuple, tzrphase)
    cwtoa::WidebandTOACorrection = correct_toa(model, wtoa, params)

    dphase = GQ{Float64}(phase_residual(wtoa.toa, cwtoa.toa_correction) - tzrphase)
    tres = dphase / doppler_shifted_spin_frequency(cwtoa.toa_correction)
    err2 = scaled_toa_error_sqr(wtoa.toa, cwtoa.toa_correction)
    chi2_term_toa = value(tres * tres / err2)

    dmres = dm_residual(wtoa.dminfo, cwtoa.dm_correction)
    dmerr2 = scaled_dm_error_sqr(wtoa.dminfo, cwtoa.dm_correction)
    chi2_term_dm = value(dmres * dmres / dmerr2)

    return chi2_term_toa + chi2_term_dm
end
