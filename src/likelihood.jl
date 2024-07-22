function _chi2_term(model::TimingModel, toa::TOA, params::NamedTuple, tzrphase)
    ctoa = correct_toa(model, toa, params)
    dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
    tres = dphase / doppler_shifted_spin_frequency(ctoa)
    err2 = scaled_toa_error_sqr(ctoa)
    return value(tres * tres / err2)
end

"""Compute the χ^2 value for a given timing model and collection of TOAs (parallel execution)."""
function calc_chi2(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    chisq = Atomic{Float64}(0.0)
    @threads for toa in toas
        atomic_add!(chisq, _chi2_term(model, toa, params, tzrphase))
    end

    return chisq[]
end

calc_chi2(model::TimingModel, toas::Vector{TOA}, params) =
    calc_chi2(model, toas, read_params(model, params))

"""Compute the χ^2 value for a given timing model and collection of TOAs (serial execution)."""
function calc_chi2_serial(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)
    return sum(toa -> _chi2_term(model, toa, params, tzrphase), toas)
end

calc_chi2_serial(model::TimingModel, toas::Vector{TOA}, params) =
    calc_chi2_serial(model, toas, read_params(model, params))

function _lnlike_term(model::TimingModel, toa::TOA, params::NamedTuple, tzrphase)
    ctoa = correct_toa(model, toa, params)
    dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
    tres = dphase / doppler_shifted_spin_frequency(ctoa)
    err2 = scaled_toa_error_sqr(ctoa)
    norm = log(value(err2))
    return value(tres * tres / err2) + norm
end

"""Compute the log-likelihood value for a given timing model and collection of TOAs (parallel execution)."""
function calc_lnlike(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    result = Atomic{Float64}(0.0)
    @threads for toa in toas
        atomic_add!(result, _lnlike_term(model, toa, params, tzrphase))
    end

    return -result[] / 2
end

calc_lnlike(model::TimingModel, toas::Vector{TOA}, params) =
    calc_lnlike(model, toas, read_params(model, params))

"""Compute the log-likelihood value for a given timing model and collection of TOAs (serial execution)."""
function calc_lnlike_serial(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)
    return -sum(toa -> _lnlike_term(model, toa, params, tzrphase), toas) / 2
end

calc_lnlike_serial(model::TimingModel, toas::Vector{TOA}, params) =
    calc_lnlike_serial(model, toas, read_params(model, params))
