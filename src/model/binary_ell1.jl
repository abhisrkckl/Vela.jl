export BinaryELL1

struct BinaryELL1 <: BinaryComponent
    use_fbx::Bool
end

orbital_phase(Δt::GQ, params::NamedTuple, use_fbx::Bool) =
    use_fbx ? taylor_horner_integral(Δt, params.FB, dimensionless(0.0)) :
    (Δt / params.PB - 0.5 * (params.PBDOT + params.XPBDOT) * (Δt / params.PB)^2)

orbital_frequency(Δt::GQ, params::NamedTuple, use_fbx::Bool) =
    use_fbx ? taylor_horner(Δt, params.FB) :
    1 / (params.PB + (params.PBDOT + params.XPBDOT) * Δt)

function rømer_delay(Φ_trigs, a1, ϵ1, ϵ2)
    (sinΦ, cosΦ), (sin2Φ, cos2Φ), (sin3Φ, cos3Φ), (sin4Φ, cos4Φ) = Φ_trigs

    return a1 * (
        sinΦ + 0.5 * (ϵ2 * sin2Φ - ϵ1 * cos2Φ) -
        (1 / 8) * (
            (5 * ϵ2^2 + 3 * ϵ1^2) * sinΦ - 2 * ϵ2 * ϵ1 * cosΦ +
            (-3 * ϵ2^2 + 3 * ϵ1^2) * sin3Φ +
            6 * ϵ2 * ϵ1 * cos3Φ
        ) -
        (1 / 12) * (
            (5 * ϵ2^3 + 3 * ϵ1^2 * ϵ2) * sin2Φ +
            (-6 * ϵ1 * ϵ2^2 - 4 * ϵ1^3) * cos2Φ +
            (-4 * ϵ2^3 + 12 * ϵ1^2 * ϵ2) * sin4Φ +
            (12 * ϵ1 * ϵ2^2 - 4 * ϵ1^3) * cos4Φ
        )
    )
end

function d_rømer_delay_d_Φ(Φ_trigs, a1, ϵ1, ϵ2)
    (sinΦ, cosΦ), (sin2Φ, cos2Φ), (sin3Φ, cos3Φ), (sin4Φ, cos4Φ) = Φ_trigs

    return a1 * (
        cosΦ + ϵ1 * sin2Φ + ϵ2 * cos2Φ -
        (1.0 / 8) * (
            (5 * ϵ2^2 + 3 * ϵ1^2) * cosΦ +
            2 * ϵ1 * ϵ2 * sinΦ +
            (-9 * ϵ2^2 + 9 * ϵ1^2) * cos3Φ - 18 * ϵ1 * ϵ2 * sin3Φ
        ) -
        (1.0 / 12) * (
            (10 * ϵ2^3 + 6 * ϵ1^2 * ϵ2) * cos2Φ +
            (12 * ϵ1 * ϵ2^2 + 8 * ϵ1^3) * sin2Φ +
            (-16 * ϵ2^3 + 48 * ϵ1^2 * ϵ2) * cos4Φ +
            (-48 * ϵ1 * ϵ2^2 + 16 * ϵ1^3) * sin4Φ
        )
    )
end

function d2_rømer_delay_d_Φ2(Φ_trigs, a1, ϵ1, ϵ2)
    (sinΦ, cosΦ), (sin2Φ, cos2Φ), (sin3Φ, cos3Φ), (sin4Φ, cos4Φ) = Φ_trigs

    return a1 * (
        -sinΦ + 2 * ϵ1 * cos2Φ - 2 * ϵ2 * sin2Φ -
        (1.0 / 8) * (
            (-5 * ϵ2^2 - 3 * ϵ1^2) * sinΦ +
            2 * ϵ1 * ϵ2 * cosΦ +
            (27 * ϵ2^2 - 27 * ϵ1^2) * sin3Φ - 54 * ϵ1 * ϵ2 * cos3Φ
        ) -
        (1.0 / 12) * (
            (-20 * ϵ2^3 - 12 * ϵ1^2 * ϵ2) * sin2Φ +
            (24 * ϵ1 * ϵ2^2 + 16 * ϵ1^3) * cos2Φ +
            (64 * ϵ2^3 - 192 * ϵ1^2 * ϵ2) * sin4Φ +
            (-192 * ϵ1 * ϵ2^2 + 64 * ϵ1^3) * cos4Φ
        )
    )
end

function shapiro_delay(Φ_trigs, m2, sini)
    sinΦ = Φ_trigs[1][1]
    return -2 * m2 * log(1 - sini * sinΦ)
end

function delay(ell1::BinaryELL1, ctoa::CorrectedTOA, params::NamedTuple)
    Δt = corrected_toa_value(ctoa) - params.TASC
    a1 = params.A1 + Δt * params.A1DOT

    ϵ1 = params.EPS1 + Δt * params.EPS1DOT
    ϵ2 = params.EPS2 + Δt * params.EPS2DOT

    Φ = orbital_phase(Δt, params, ell1.use_fbx)

    Φ_trigs = sincos(Φ), sincos(2Φ), sincos(3Φ), sincos(4Φ)

    m2 = params.M2
    sini = params.SINI

    ΔR = rømer_delay(Φ_trigs, a1, ϵ1, ϵ2)
    ΔS = shapiro_delay(Φ_trigs, m2, sini)

    ΔRp = d_rømer_delay_d_Φ(Φ_trigs, a1, ϵ1, ϵ2)
    ΔRp2 = d2_rømer_delay_d_Φ2(Φ_trigs, a1, ϵ1, ϵ2)
    nhat = 2 * π * orbital_frequency(Δt, params, ell1.use_fbx)

    ΔR_inv = ΔR * (1 - nhat * ΔRp + (nhat * ΔRp)^2 + 0.5 * nhat^2 * ΔR * ΔRp2)

    return ΔR_inv + ΔS
end
