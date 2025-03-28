"""
    BinaryELL1Base
    
The abstract base type for all binary models representing nearly circular orbits.

Reference:
    [Lange+ 2001](http://doi.org/10.1046/j.1365-8711.2001.04606.x)
"""
abstract type BinaryELL1Base <: BinaryComponent end

ELL1ΦTrigs = NTuple{4,NTuple{2,GQ{0,Float64}}}

"""
    ELL1State

The instantaneous state of a nearly circular binary.

`Φ_trigs` contains the various trigonometric functions of the mean anomaly 
required for evaluating the Rømer delay and the Shapiro delay for such
binaries. `n` is the mean motion. `a1` is the projected semi-major axis of
the pulsar orbit. `ϵ1` and `ϵ2` are the Laplace-Lagrange parameters.

Reference:
    [Lange+ 2001](http://doi.org/10.1046/j.1365-8711.2001.04606.x)
"""
struct ELL1State
    Φ_trigs::ELL1ΦTrigs
    n::GQ{-1,Float64}
    a1::GQ{1,Float64}
    ϵ1::GQ{0,Float64}
    ϵ2::GQ{0,Float64}
    m2::GQ{1,Float64}
    sini::GQ{0,Float64}
end

function ELL1State(
    ell1::BinaryELL1Base,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
)
    Δt = corrected_toa_value(toa, toacorr, Float64) - params.TASC
    a1 = params.A1 + Δt * params.A1DOT
    ϵ1 = params.EPS1 + Δt * params.EPS1DOT
    ϵ2 = params.EPS2 + Δt * params.EPS2DOT

    Φ = mean_anomaly(Δt, params, ell1.use_fbx)
    Φ_trigs = sincos(Φ), sincos(2Φ), sincos(3Φ), sincos(4Φ)

    n = mean_motion(Δt, params, ell1.use_fbx)

    m2, sini = shapiro_delay_params(ell1, params)

    return ELL1State(Φ_trigs, n, a1, ϵ1, ϵ2, m2, sini)
end

"""Rømer delay due to a nearly circular binary. Includes terms up to the
cubic order in eccentricity.

Reference:
    [Lange+ 2001](http://doi.org/10.1046/j.1365-8711.2001.04606.x)
"""
function rømer_delay(::BinaryELL1Base, state::ELL1State)::GQ
    (sinΦ, cosΦ), (sin2Φ, cos2Φ), (sin3Φ, cos3Φ), (sin4Φ, cos4Φ) = state.Φ_trigs
    a1 = state.a1
    ϵ1 = state.ϵ1
    ϵ2 = state.ϵ2

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

"""Derivative of the Rømer delay w.r.t. the mean anomaly. Used for evaluating the inverse
timing formula.

Reference:
    [Lange+ 2001](http://doi.org/10.1046/j.1365-8711.2001.04606.x)
"""
function d_rømer_delay_d_Φ(::BinaryELL1Base, state::ELL1State)::GQ
    (sinΦ, cosΦ), (sin2Φ, cos2Φ), (sin3Φ, cos3Φ), (sin4Φ, cos4Φ) = state.Φ_trigs
    a1 = state.a1
    ϵ1 = state.ϵ1
    ϵ2 = state.ϵ2

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

"""Second derivative of the Rømer delay w.r.t. the mean anomaly. Used for evaluating the inverse
timing formula.

Reference:
    [Lange+ 2001](http://doi.org/10.1046/j.1365-8711.2001.04606.x)
"""
function d2_rømer_delay_d_Φ2(::BinaryELL1Base, state::ELL1State)::GQ
    (sinΦ, cosΦ), (sin2Φ, cos2Φ), (sin3Φ, cos3Φ), (sin4Φ, cos4Φ) = state.Φ_trigs
    a1 = state.a1
    ϵ1 = state.ϵ1
    ϵ2 = state.ϵ2

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

"""Shapiro delay due to a nearly circular binary.

Reference:
    [Lange+ 2001](http://doi.org/10.1046/j.1365-8711.2001.04606.x)
"""
function shapiro_delay(::BinaryELL1Base, state::ELL1State)::GQ
    sinΦ = state.Φ_trigs[1][1]
    m2 = state.m2
    sini = state.sini
    return -2 * m2 * log(1 - sini * sinΦ)
end

"""Update the `CorrectedTOA` object with delays and Doppler factor due to a nearly
circular binary.

Reference:
    [Lange+ 2001](http://doi.org/10.1046/j.1365-8711.2001.04606.x)
"""
function correct_toa(
    ell1::BinaryELL1Base,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
)
    state = ELL1State(ell1, toa, toacorr, params)

    ΔR = rømer_delay(ell1, state)
    ΔRp = d_rømer_delay_d_Φ(ell1, state)
    ΔRp2 = d2_rømer_delay_d_Φ2(ell1, state)
    nhat = state.n
    ΔR_inv = ΔR * (1 - nhat * ΔRp + nhat * nhat * ΔRp * ΔRp + 0.5 * nhat * nhat * ΔR * ΔRp2)

    ΔS = shapiro_delay(ell1, state)

    delay = ΔR_inv + ΔS

    doppler = -ΔRp * nhat

    return correct_toa_delay(toacorr; delay = delay, doppler = doppler)
end
