function _chi2_term(model::TimingModel, wtoa::WidebandTOA, params::NamedTuple, tzrphase)
    cwtoa::CorrectedWidebandTOA = correct_toa(model, wtoa, params)

    dphase = GQ{Float64}(phase_residual(cwtoa.corrected_toa) - tzrphase)
    tres = dphase / doppler_shifted_spin_frequency(cwtoa.corrected_toa)
    err2 = scaled_toa_error_sqr(cwtoa.corrected_toa)
    chi2_term_toa = value(tres * tres / err2)

    dmres = dm_residual(cwtoa.corrected_dminfo)
    dmerr2 = scaled_dm_error_sqr(cwtoa.corrected_dminfo)
    chi2_term_dm = value(dmres * dmres / dmerr2)

    return chi2_term_toa + chi2_term_dm
end
