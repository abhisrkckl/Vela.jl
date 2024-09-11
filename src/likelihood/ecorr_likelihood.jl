function _ecorr_lnlike_group(
    model::TimingModel,
    toas::Vector{TOA},
    params::NamedTuple,
    tzrphase::GQ,
    group::EcorrGroup,
)
    w = (group.index == 0) ? time(0.0) : params.ECORR[group.index]^Val(2)

    r_r = dimensionless(0.0)
    r_u = GQ{-1}(0.0)
    u_u = GQ{-2}(0.0)
    norm = 0.0

    for ii = group.start:group.stop
        toa = toas[ii]
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(ctoa) - tzrphase)
        tres = dphase / doppler_shifted_spin_frequency(ctoa)
        err2 = scaled_toa_error_sqr(ctoa)

        r_r += tres * tres / err2
        r_u += tres / err2
        u_u += 1 / err2
        norm += log(value(err2))
    end

    return r_r - w * r_u * r_u / (1 + w * u_u) + norm + 2 * log(1 + w * u_u)
end

_ecorr_lnlike_chunk(
    model::TimingModel,
    toas::Vector{T},
    params,
    tzrphase,
    chunk,
) where {T<:TOABase} = sum(
    jj -> _ecorr_lnlike_group(model, toas, params, tzrphase, model.kernel.ecorr_groups[jj]),
    chunk,
)

"""Compute the log-likelihood value for a given timing model and collection of TOAs 
(parallel execution).

Reference:
    [Lentati+ 2014](https://doi.org/10.1093/mnras/stt2122)
"""
function calc_lnlike(
    model::TimingModel{ComponentsTuple,EcorrKernel,PriorsTuple},
    toas::Vector{T},
    params::NamedTuple,
) where {ComponentsTuple<:Tuple,PriorsTuple<:Tuple,T<:TOABase}
    tzrphase = calc_tzr_phase(model, params)
    groups = model.kernel.ecorr_groups
    chunks = Iterators.partition(eachindex(groups), length(groups) ÷ nthreads())
    spawn_chunk(chunk) = @spawn _ecorr_lnlike_chunk(model, toas, params, tzrphase, chunk)
    tasks = map(spawn_chunk, chunks)
    result = sum(fetch, tasks)
    return -result / 2
end

function calc_lnlike_serial(
    model::TimingModel{ComponentsTuple,EcorrKernel,PriorsTuple},
    toas::Vector{T},
    params::NamedTuple,
) where {ComponentsTuple<:Tuple,PriorsTuple<:Tuple,T<:TOABase}
    tzrphase = calc_tzr_phase(model, params)
    return -sum(
        group -> _ecorr_lnlike_group(model, toas, params, tzrphase, group),
        model.kernel.ecorr_groups,
    ) / 2
end
