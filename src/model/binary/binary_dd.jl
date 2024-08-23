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
    sinω, cosω = sincos(ω)

    a1 = params.A1 + Δt * params.A1DOT

    α = a1 * sinω
    β = a1 * sqrt(1 - eϕ * eϕ) * cosω
    γ = params.GAMMA

    m2 = params.M2
    sini = params.SINI

    return DDState((α, β, γ), sincos(u), et, er, a1, n, m2, sini)
end
