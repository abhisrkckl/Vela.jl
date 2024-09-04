export calc_chi2, calc_lnlike, get_chi2_serial_func, get_chi2_parallel_func, get_chi2_func

function _chi2_term(model::TimingModel, toa::TOA, params::NamedTuple, tzrphase)
    ctoa = correct_toa(model, toa, params)
    dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
    tres = dphase / doppler_shifted_spin_frequency(ctoa)
    err2 = scaled_toa_error_sqr(ctoa)
    return value(tres * tres / err2)
end

_chi2_chunk(
    model::TimingModel,
    toas::Vector{T},
    params,
    tzrphase,
    chunk,
) where {T<:Union{TOA,WidebandTOA}} =
    sum(ii -> _chi2_term(model, toas[ii], params, tzrphase), chunk)

"""Compute the χ^2 value for a given timing model and collection of TOAs (parallel execution)."""
function calc_chi2(
    model::TimingModel,
    toas::Vector{T},
    params::NamedTuple,
) where {T<:Union{TOA,WidebandTOA}}
    tzrphase = calc_tzr_phase(model, params)
    chunks = Iterators.partition(eachindex(toas), length(toas) ÷ nthreads())
    spawn_chunk(chunk) = @spawn _chi2_chunk(model, toas, params, tzrphase, chunk)
    tasks = map(spawn_chunk, chunks)
    result = sum(fetch, tasks)
    return result
end

calc_chi2(model::TimingModel, toas::Vector{T}, params) where {T<:Union{TOA,WidebandTOA}} =
    calc_chi2(model, toas, read_params(model, params))

"""Compute the χ^2 value for a given timing model and collection of TOAs (serial execution)."""
function calc_chi2_serial(
    model::TimingModel,
    toas::Vector{T},
    params::NamedTuple,
) where {T<:Union{TOA,WidebandTOA}}
    tzrphase = calc_tzr_phase(model, params)
    return sum(toa -> _chi2_term(model, toa, params, tzrphase), toas)
end

calc_chi2_serial(
    model::TimingModel,
    toas::Vector{T},
    params,
) where {T<:Union{TOA,WidebandTOA}} =
    calc_chi2_serial(model, toas, read_params(model, params))

function get_chi2_serial_func(
    model::TimingModel,
    toas::Vector{T},
) where {T<:Union{TOA,WidebandTOA}}
    toas_ = copy(toas)
    params -> calc_chi2_serial(model, toas_, params)
end

function get_chi2_parallel_func(
    model::TimingModel,
    toas::Vector{T},
) where {T<:Union{TOA,WidebandTOA}}
    toas_ = copy(toas)
    params -> calc_chi2(model, toas_, params)
end

"""Get the χ^2(params) function for a given timing model and collection of TOAs.
Serial or parallel execution is decided based on the number of available threads."""
get_chi2_func(model, toas) =
    (nthreads() == 1) ? get_chi2_serial_func(model, toas) :
    get_chi2_parallel_func(model, toas)
