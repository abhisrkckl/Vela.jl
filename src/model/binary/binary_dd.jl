struct BinaryDD <: BinaryDDBase
    use_fbx::Bool
end

function DDState(dd::BinaryDD, ctoa::CorrectedTOA, params::NamedTuple)
    Δt = corrected_toa_value(ctoa) - params.T0
    n = mean_motion(Δt, params, dd.use_fbx)
    l = mean_anomaly(Δt, params, dd.use_fbx)

    et = params.ECC + Δt * params.EDOT
    er = et * (1 + params.DR)
    eϕ = et * (1 + params.DTH)

    u = mikkola(l, et)
    v = true_anomaly(u, eϕ)

    k = params.OMDOT / n
    ω = params.OM + k * v

    a1 = params.A1 + Δt * params.A1DOT

    return DDState(sincos(u), sincos(ω), et, er, eϕ, a1)
end
