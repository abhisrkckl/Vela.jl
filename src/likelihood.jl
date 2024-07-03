export calc_chi2, calc_lnlike, get_lnlike_func

function calc_chi2(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    chisq = Atomic{Float64}(0.0)
    @threads for toa in toas
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
        tres = dphase / doppler_shifted_spin_frequency(ctoa)
        err2 = scaled_toa_error_sqr(ctoa)

        atomic_add!(chisq, Float64((tres * tres / err2).x))
    end

    return chisq[]
end

function calc_chi2_serial(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    chisq = 0.0
    for toa in toas
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
        tres = dphase / doppler_shifted_spin_frequency(ctoa)
        err2 = scaled_toa_error_sqr(ctoa)

        chisq += (tres * tres / err2).x
    end

    return chisq
end

read_params(model::TimingModel, values::Vector{Float64}) =
    return read_params(model.param_handler, values)

calc_chi2(model::TimingModel, toas::Vector{TOA}, params::Vector{Float64}) =
    calc_chi2(model, toas, read_params(model, params))

calc_chi2_serial(model::TimingModel, toas::Vector{TOA}, params::Vector{Float64}) =
    calc_chi2_serial(model, toas, read_params(model, params))

function calc_lnlike(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    result = Atomic{Float64}(0.0)
    @threads for toa in toas
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
        tres = dphase / doppler_shifted_spin_frequency(ctoa)
        err2 = scaled_toa_error_sqr(ctoa)
        norm = log(err2.x) / 2

        atomic_add!(result, (tres * tres / err2).x + norm)
    end

    return -result[] / 2
end

function calc_lnlike_serial(model::TimingModel, toas::Vector{TOA}, params::NamedTuple)
    tzrphase = calc_tzr_phase(model, params)

    result = 0.0
    for toa in toas
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
        tres = dphase / doppler_shifted_spin_frequency(ctoa)
        err2 = scaled_toa_error_sqr(ctoa)
        norm = log(err2.x) / 2

        result += (tres * tres / err2).x + norm
    end

    return -result / 2
end

calc_lnlike(model::TimingModel, toas::Vector{TOA}, params::Vector{Float64}) =
    calc_lnlike(model, toas, read_params(model, params))

calc_lnlike_serial(model::TimingModel, toas::Vector{TOA}, params::Vector{Float64}) =
    calc_lnlike_serial(model, toas, read_params(model, params))

get_lnlike_func(model::TimingModel, toas::Vector{TOA}) =
    params -> calc_lnlike(model, toas, params)
get_lnlike_serial_func(model::TimingModel, toas::Vector{TOA}) =
    params -> calc_lnlike_serial(model, toas, params)
