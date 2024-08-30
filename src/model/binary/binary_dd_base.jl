"""The abstract base type for all binary models representing eccentric orbits."""
abstract type BinaryDDBase <: BinaryComponent end

struct DDState
    rømer_einstein_coeffs::NTuple{3,GQ{1,Float64}}
    sincosu::SinCos
    et::GQ{0,Float64}
    er::GQ{0,Float64}
    a1::GQ{1,Float64}
    n::GQ{-1,Float64}
    m2::GQ{1,Float64}
    sini::GQ{0,Float64}
end

function rømer_einstein_delay(::BinaryDDBase, state::DDState)::GQ
    sinu, cosu = state.sincosu
    er = state.er
    α, β, γ = state.rømer_einstein_coeffs

    return α * (cosu - er) + (β + γ) * sinu
end

function d_rømer_einstein_delay_d_u(::BinaryDDBase, state::DDState)::GQ
    sinu, cosu = state.sincosu
    er = state.er
    α, β, γ = state.rømer_einstein_coeffs

    return -α * sinu + (β + γ) * cosu
end

function d2_rømer_einstein_delay_d_u2(::BinaryDDBase, state::DDState)::GQ
    sinu, cosu = state.sincosu
    er = state.er
    α, β, γ = state.rømer_einstein_coeffs

    return -α * cosu - (β + γ) * sinu
end

function shapiro_delay(::BinaryDDBase, state::DDState)
    m2 = state.m2
    sini = state.sini
    sinu, cosu = state.sincosu

    er = state.er
    et = state.et
    a1 = state.a1
    α, β, _ = state.rømer_einstein_coeffs

    return -2 * m2 * (1 - et * cosu - (sini / a1) * (α * (cosu - er) + β * sinu))
end

"""Total delay due to a nearly circular binary."""
function correct_toa(dd::BinaryDDBase, ctoa::CorrectedTOA, params::NamedTuple)
    state = DDState(dd, ctoa, params)

    ΔRE = rømer_einstein_delay(dd, state)
    ΔREp = d_rømer_einstein_delay_d_u(dd, state)
    ΔREp2 = d2_rømer_einstein_delay_d_u2(dd, state)

    et = state.et
    n = state.n
    sinu, cosu = state.sincosu
    nhat = n / (1 - et * cosu)

    ΔRE_inv =
        ΔRE * (
            1 - nhat * ΔREp +
            (nhat * nhat * ΔREp * ΔREp) +
            0.5 * nhat * nhat * ΔRE * ΔREp2 -
            0.5 * et * sinu / (1 - et * cosu) * nhat * nhat * ΔRE * ΔREp
        )

    ΔS = shapiro_delay(dd, state)

    delay = ΔRE_inv + ΔS

    # Is this accurate enough?
    doppler = ΔREp * nhat

    return correct_toa(ctoa; delay = delay, doppler = doppler)
end
