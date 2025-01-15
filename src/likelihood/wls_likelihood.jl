export calc_lnlike,
    calc_lnlike_serial, get_lnlike_serial_func, get_lnlike_parallel_func, get_lnlike_func

const ln_2π = log(2π)

"""A single term in the pulsar timing log-likelihood expression (white noise-only).

A factor of 1/2 is excluded here."""
function _wls_lnlike_term(model::TimingModel, toa::TOA, params::NamedTuple, tzrphase)
    ctoa = correct_toa(model, toa, params)
    dphase = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)
    tres = dphase / doppler_shifted_spin_frequency(ctoa)
    err2 = scaled_toa_error_sqr(toa, ctoa)
    norm = log(value(err2))
    return value(tres * tres / err2) + norm
end

_wls_lnlike_chunk(
    model::TimingModel,
    toas::Vector{T},
    params,
    tzrphase,
    chunk,
) where {T<:TOABase} = sum(ii -> _wls_lnlike_term(model, toas[ii], params, tzrphase), chunk)

data_length(toas::Vector{TOA}) = length(toas)
data_length(toas::Vector{WidebandTOA}) = 2 * length(toas)

function calc_lnlike(
    model::TimingModel{ComponentsTuple,WhiteNoiseKernel,PriorsTuple},
    toas::Vector{T},
    params::NamedTuple,
) where {ComponentsTuple<:Tuple,PriorsTuple<:Tuple,T<:TOABase}
    tzrphase = calc_tzr_phase(model, params)
    chunks = Iterators.partition(eachindex(toas), length(toas) ÷ nthreads())
    spawn_chunk(chunk) = @spawn _wls_lnlike_chunk(model, toas, params, tzrphase, chunk)
    tasks = map(spawn_chunk, chunks)
    result = sum(fetch, tasks)
    return -result / 2 - (data_length(toas) * ln_2π / 2)
end

"""
    calc_lnlike(::TimingModel, ::Vector{T}, params)::Float64 where {T<:TOABase}

Compute the log-likelihood value for a given timing model and collection of TOAs 
(parallel execution).

Reference:
    [Lentati+ 2014](https://doi.org/10.1093/mnras/stt2122)
"""
calc_lnlike(model::TimingModel, toas::Vector{T}, params) where {T<:TOABase} =
    calc_lnlike(model, toas, read_params(model, params))

function calc_lnlike_serial(
    model::TimingModel{ComponentsTuple,WhiteNoiseKernel,PriorsTuple},
    toas::Vector{T},
    params::NamedTuple,
) where {ComponentsTuple<:Tuple,PriorsTuple<:Tuple,T<:TOABase}
    tzrphase = calc_tzr_phase(model, params)
    return -sum(toa -> _wls_lnlike_term(model, toa, params, tzrphase), toas) / 2 -
           (data_length(toas) * ln_2π / 2)
end

"""
    calc_lnlike_serial(::TimingModel, ::Vector{T}, params)::Float64 where {T<:TOABase}

Compute the log-likelihood value for a given timing model and collection of TOAs 
(serial execution).

Reference:
    [Lentati+ 2014](https://doi.org/10.1093/mnras/stt2122)
"""
calc_lnlike_serial(model::TimingModel, toas::Vector{T}, params) where {T<:TOABase} =
    calc_lnlike_serial(model, toas, read_params(model, params))

"""
    get_lnlike_serial_func(model, toas)::Function

Version of `get_lnlike_func` that always does serial execution.
"""
function get_lnlike_serial_func(model::TimingModel, toas::Vector{T}) where {T<:TOABase}
    toas_ = copy(toas)
    params -> calc_lnlike_serial(model, toas_, params)
end

function get_lnlike_parallel_func(model::TimingModel, toas::Vector{T}) where {T<:TOABase}
    toas_ = copy(toas)
    params -> calc_lnlike(model, toas_, params)
end

"""
    get_lnlike_func(model, toas)::Function

Get the log_likelihood(params) function for a given timing model and collection of TOAs.
Serial or parallel execution is decided based on the number of available threads.

Supports both narrowband and wideband TOAs.

Use `get_lnlike_serial_func(model, toas)` to force serial execution of the likelihood. 
The serial version should be used if parallelization is to be implemented at a different level 
(e.g., within the sampling method).

Reference:
    [Lentati+ 2014](https://doi.org/10.1093/mnras/stt2122),
    [Alam+ 2021](http://doi.org/10.3847/1538-4365/abc6a1),
    [Johnson+ 2024](https://doi.org/10.1103/PhysRevD.109.103012)
"""
get_lnlike_func(model, toas) =
    (nthreads() == 1) ? get_lnlike_serial_func(model, toas) :
    get_lnlike_parallel_func(model, toas)
