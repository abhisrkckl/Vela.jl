"""A single block in the pulsar timing log-likelihood expression when white noise and ECORR are present.

A factor of 1/2 is excluded here.

Reference:
    [Johnson+ 2024](https://doi.org/10.1103/PhysRevD.109.103012)
"""
function _ecorr_lnlike_group(
    model::TimingModel,
    toas::Vector{TOA},
    params::NamedTuple,
    tzrphase::GQ,
    group::EcorrGroup,
)
    ecorr = (group.index == 0) ? time(0.0) : params.ECORR[group.index]
    w = ecorr * ecorr

    r_r = dimensionless(0.0)
    r_u = GQ{-1}(0.0)
    u_u = GQ{-2}(0.0)
    norm = 0.0

    for ii in group.start:group.stop
        toa = toas[ii]
        ctoa = correct_toa(model, toa, params)
        dphase = GQ{Float64}(phase_residual(toa, ctoa) - tzrphase)
        tres = dphase / doppler_shifted_spin_frequency(ctoa)
        err2 = scaled_toa_error_sqr(toa, ctoa)

        r_r += tres * tres / err2
        r_u += tres / err2
        u_u += 1 / err2
        norm += log(value(err2))
    end

    return r_r - w * r_u * r_u / (1 + w * u_u) + norm + log(1 + w * u_u)
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
    return -value(result) / 2
end

function calc_lnlike_serial(
    model::TimingModel{ComponentsTuple,EcorrKernel,PriorsTuple},
    toas::Vector{T},
    params::NamedTuple,
) where {ComponentsTuple<:Tuple,PriorsTuple<:Tuple,T<:TOABase}
    tzrphase = calc_tzr_phase(model, params)
    return -value(
        sum(
            group -> _ecorr_lnlike_group(model, toas, params, tzrphase, group),
            model.kernel.ecorr_groups,
        ),
    ) / 2
end
