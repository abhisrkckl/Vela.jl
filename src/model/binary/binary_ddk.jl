export BinaryDDK

"""
    BinaryDDK

The Damour & Deruelle model for eccentric binaries with Kopeikin corrections included,
which account for the apparent changes in the orbital elements due to proper motion and 
parallax.

References:
    [Damour & Deruelle 1986](https://ui.adsabs.harvard.edu/abs/1986AIHPA..44..263D/abstract),
    [Kopeikin 1995](http://doi.org/10.1086/187731),
    [Kopeikin 1996](http://doi.org/10.1086/310201)
"""
struct BinaryDDK <: BinaryDDBase
    use_fbx::Bool
    ecliptic_coordinates::Bool
end

function DDState(ddk::BinaryDDK, toa::TOA, toacorr::TOACorrection, params::NamedTuple)
    Δt = corrected_toa_value(toa, toacorr, Float64) - params.T0
    n = mean_motion(Δt, params, ddk.use_fbx)
    l = mean_anomaly(Δt, params, ddk.use_fbx)

    et = params.ECC + Δt * params.EDOT

    er = et * (1 + params.DR)
    eϕ = et * (1 + params.DTH)

    u = mikkola(l, et)
    sinu, cosu = sincos(u)

    ηϕ = sqrt(1 - eϕ * eϕ)

    βϕ = (1 - ηϕ) / eϕ
    v_u = 2 * atan(βϕ * sinu, 1 - βϕ * cosu)
    v = v_u + u

    μα, μδ = proper_motion(params, ddk.ecliptic_coordinates)
    px = params.PX

    ssb_obs_pos = toa.ephem.ssb_obs_pos
    ssb_psr_pos = toacorr.ssb_psr_pos

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

"""Compute the vectors I0 and J0 that appear in the Kopeikin parallax corrections.

References:
    [Kopeikin 1995](http://doi.org/10.1086/187731)
"""
function kopeikin_I0_J0(ssb_psr_pos)
    sinδ = ssb_psr_pos[3]
    cosδ = sqrt(1 - sinδ * sinδ)
    cosα = ssb_psr_pos[1] / cosδ
    sinα = ssb_psr_pos[2] / cosδ

    I0 = (-sinα, cosα, zero(sinα))
    J0 = (-cosα * sinδ, -sinα * sinδ, cosδ)

    return I0, J0
end

"""Evaluate the Kopeikin corrections to the projected semi-major axis, argument of
periapsis, and the inclination due to proper motion and parallax.

References:
    [Kopeikin 1995](http://doi.org/10.1086/187731),
    [Kopeikin 1996](http://doi.org/10.1086/310201)
"""
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
