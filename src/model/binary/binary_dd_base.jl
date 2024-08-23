"""The abstract base type for all binary models representing eccentric orbits."""
abstract type BinaryDDBase <: BinaryComponent end

struct DDState
    rømer_einstein_coeffs::NTuple{3,GQ{Float64}}
    sincosu::SinCos
    et::GQ{Float64}
    er::GQ{Float64}
    a1::GQ{Float64}
    n::GQ{Float32}
    m2::GQ{Float32}
    sini::GQ{Float32}
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

    et = state.et
    a1 = state.a1
    α, β, _ = state.rømer_einstein_coeffs

    return -2 * m2 * (1 - et * cosu - (sini / a1) * (α * (cosu - er) + β * sinu))
end

"""Total delay due to a nearly circular binary."""
function delay(dd::BinaryDDBase, ctoa::CorrectedTOA, params::NamedTuple)::GQ
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
            1 - nhat * ΔREp + (nhat * ΔREp)^2 + 0.5 * nhat^2 * ΔRE * ΔREp2 -
            0.5 * et * sinu / (1 - et * cosu) * nhat^2 * ΔRE * ΔREp
        )

    ΔS = shapiro_delay(dd, state, params)

    return ΔRE_inv + ΔS
end
