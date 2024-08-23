"""The abstract base type for all binary models representing eccentric orbits."""
abstract type BinaryDDBase <: BinaryComponent end

struct DDState
    sincosu::SinCos
    sincosω::SinCos
    et::GQ{Float64}
    er::GQ{Float64}
    eϕ::GQ{Float64}
    a1::GQ{Float64}
end

"""Total delay due to a nearly circular binary."""
function delay(dd::BinaryDDBase, ctoa::CorrectedTOA, params::NamedTuple)
    state = DDState(dd, ctoa, params)

    ΔR = rømer_delay(dd, state)
    ΔRp = d_rømer_delay_d_Φ(dd, state)
    ΔRp2 = d2_rømer_delay_d_Φ2(dd, state)
    nhat = compute_nhat(state)
    ΔR_inv = ΔR * (1 - nhat * ΔRp + (nhat * ΔRp)^2 + 0.5 * nhat^2 * ΔR * ΔRp2)

    ΔS = shapiro_delay(dd, state, params)
    ΔE = einstein_delay(dd, state, params)

    return ΔR_inv + ΔS + ΔE
end

function ()

end
