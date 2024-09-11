function _wls_lnlike_term(
    model::TimingModel,
    wtoa::WidebandTOA,
    params::NamedTuple,
    tzrphase,
)
    cwtoa = correct_toa(model, wtoa, params)

    dphase = GQ{Float64}(phase_residual(cwtoa.corrected_toa) - tzrphase)
    tres = dphase / doppler_shifted_spin_frequency(cwtoa.corrected_toa)
    err2 = scaled_toa_error_sqr(cwtoa.corrected_toa)
    norm = log(value(err2))
    lnlike_term_toa = value(tres * tres / err2) + norm

    dmres = dm_residual(cwtoa.corrected_dminfo)
    dmerr2 = scaled_dm_error_sqr(cwtoa.corrected_dminfo)
    norm = log(value(dmerr2))
    lnlike_term_dm = value(dmres * dmres / dmerr2) + norm

    return lnlike_term_toa + lnlike_term_dm
end
