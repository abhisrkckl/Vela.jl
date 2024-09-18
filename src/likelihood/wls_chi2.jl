export calc_chi2,
    calc_lnlike,
    get_chi2_serial_func,
    get_chi2_parallel_func,
    get_chi2_func,
    degrees_of_freedom,
    calc_chi2_reduced

"""A single term in the pulsar timing χ^2 expression (white noise-only)."""
function _wls_chi2_term(model::TimingModel, toa::TOA, params::NamedTuple, tzrphase)
    ctoa = correct_toa(model, toa, params)
    dphase = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)
    tres = dphase / doppler_shifted_spin_frequency(ctoa)
    err2 = scaled_toa_error_sqr(toa, ctoa)
    return value(tres * tres / err2)
end

_wls_chi2_chunk(
    model::TimingModel,
    toas::Vector{T},
    params,
    tzrphase,
    chunk,
) where {T<:TOABase} = sum(ii -> _wls_chi2_term(model, toas[ii], params, tzrphase), chunk)

"""Compute the χ^2 value for a given timing model and collection of TOAs (parallel execution)."""
function calc_chi2(
    model::TimingModel{ComponentsTuple,WhiteNoiseKernel,PriorsTuple},
    toas::Vector{T},
    params::NamedTuple,
) where {ComponentsTuple<:Tuple,PriorsTuple<:Tuple,T<:TOABase}
    tzrphase = calc_tzr_phase(model, params)
    chunks = Iterators.partition(eachindex(toas), length(toas) ÷ nthreads())
    spawn_chunk(chunk) = @spawn _wls_chi2_chunk(model, toas, params, tzrphase, chunk)
    tasks = map(spawn_chunk, chunks)
    result = sum(fetch, tasks)
    return result
end

calc_chi2(model::TimingModel, toas::Vector{T}, params) where {T<:TOABase} =
    calc_chi2(model, toas, read_params(model, params))

"""Compute the χ^2 value for a given timing model and collection of TOAs (serial execution)."""
function calc_chi2_serial(
    model::TimingModel{ComponentsTuple,WhiteNoiseKernel,PriorsTuple},
    toas::Vector{T},
    params::NamedTuple,
) where {ComponentsTuple<:Tuple,PriorsTuple<:Tuple,T<:TOABase}
    tzrphase = calc_tzr_phase(model, params)
    return sum(toa -> _wls_chi2_term(model, toa, params, tzrphase), toas)
end

calc_chi2_serial(model::TimingModel, toas::Vector{T}, params) where {T<:TOABase} =
    calc_chi2_serial(model, toas, read_params(model, params))

function get_chi2_serial_func(model::TimingModel, toas::Vector{T}) where {T<:TOABase}
    toas_ = copy(toas)
    params -> calc_chi2_serial(model, toas_, params)
end

function get_chi2_parallel_func(model::TimingModel, toas::Vector{T}) where {T<:TOABase}
    toas_ = copy(toas)
    params -> calc_chi2(model, toas_, params)
end

"""Get the χ^2(params) function for a given timing model and collection of TOAs.
Serial or parallel execution is decided based on the number of available threads.

Supports both narrowband and wideband TOAs.

Reference:
    [Hobbs+ 2006](http://doi.org/10.1111/j.1365-2966.2006.10302.x),
    [Alam+ 2021](http://doi.org/10.3847/1538-4365/abc6a1),
    [Johnson+ 2024](https://doi.org/10.1103/PhysRevD.109.103012)
"""
get_chi2_func(model, toas) =
    (nthreads() == 1) ? get_chi2_serial_func(model, toas) :
    get_chi2_parallel_func(model, toas)

degrees_of_freedom(model::TimingModel, toas::Vector{TOA}) =
    length(toas) - model.param_handler._nfree
degrees_of_freedom(model::TimingModel, wtoas::Vector{WidebandTOA}) =
    2 * length(wtoas) - model.param_handler._nfree

calc_chi2_reduced(model::TimingModel, toas::Vector{T}, params) where {T<:TOABase} =
    calc_chi2(model, toas, params) / degrees_of_freedom(model, toas)
