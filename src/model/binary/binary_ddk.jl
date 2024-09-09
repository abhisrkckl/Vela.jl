export BinaryDDK

struct BinaryDDK <: BinaryDDBase
    use_fbx::Bool
    ecliptic_coordinates::Bool
end

function DDState(ddk::BinaryDDK, ctoa::CorrectedTOA, params::NamedTuple)
    Δt = corrected_toa_value(ctoa) - params.T0
    n = mean_motion(Δt, params, ddk.use_fbx)
    l = mean_anomaly(Δt, params, ddk.use_fbx)

    et = params.ECC + Δt * params.EDOT

    er = et * (1 + params.DR)
    eϕ = et * (1 + params.DTH)

    u = mikkola(l, et)
    sinu, cosu = sincos(u)

    ηϕ = sqrt(1 - eϕ^2)

    βφ = (1 - ηϕ) / eϕ
    v_u = 2 * atan(βφ * sinu, 1 - βφ * cosu)
    v = v_u + u

    μα, μδ = proper_motion(params, ddk.ecliptic_coordinates)
    px = params.PX

    ssb_obs_pos = ctoa.toa.ephem.ssb_obs_pos
    ssb_psr_pos = ctoa.ssb_psr_pos

    a1 = params.A1 + Δt * params.A1DOT

    ι = params.KIN
    Ω = params.KOM

    δx, δω, δι = kopeikin_corrections(Δt, a1, ι, Ω, μα, μδ, px, ssb_obs_pos, ssb_psr_pos)

    a1 += δx
    ι += δι

    k = params.OMDOT / n
    ω = params.OM + k * v + δω
    sinω, cosω = sincos(ω)

    α = a1 * sinω
    β = a1 * ηϕ * cosω
    γ = params.GAMMA

    m2 = params.M2
    sinι = sin(ι)

    return DDState((α, β, γ), (sinu, cosu), et, er, a1, n, m2, sinι)
end

proper_motion(params::NamedTuple, ecliptic_coordinates::Bool)::NTuple{2,GQ{-1,Float64}} =
    ecliptic_coordinates ? (params.PMELONG, params.PMELAT) : (params.PMRA, params.PMDEC)

function kopeikin_I0_J0(ssb_psr_pos)
    sinδ = ssb_psr_pos[3]
    cosδ = sqrt(1 - sinδ * sinδ)
    cosα = ssb_psr_pos[1] / sinδ
    sinα = ssb_psr_pos[2] / sinδ

    I0 = (-sinα, cosα, zero(sinα))
    J0 = (-cosα * sinδ, -sinα * sinδ, cosδ)

    return I0, J0
end

function kopeikin_corrections(dt, x, ι, Ω, μα, μδ, px, ssb_obs_pos, ssb_psr_pos)
    sinι, cosι = sincos(ι)
    sinΩ, cosΩ = sincos(Ω)

    cotι = cosι / sinι
    cscι = 1 / sinι

    δι_pm = (-μα * sinΩ + μδ * cosΩ) * dt
    δx_pm = x * cotι * δι_pm
    δω_pm = cscι * (μα * cosΩ + μδ * sinΩ) * dt

    I0, J0 = kopeikin_I0_J0(ssb_psr_pos)

    ΔI0 = dot(ssb_obs_pos, I0)
    ΔJ0 = dot(ssb_obs_pos, J0)

    δx_px = x * cotι * px * (ΔI0 * sinΩ - ΔJ0 * cosΩ)
    δω_px = -cscι * px * (ΔI0 * cosΩ + ΔJ0 * sinΩ)

    δx = δx_pm + δx_px
    δω = δω_pm + δω_px
    δι = δι_pm

    return δx, δω, δι
end
