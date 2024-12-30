export BinaryELL1k

"""A binary model representing a nearly circular orbit with large advance of 
periapsis. The difference between this model and ELL1 is that this model features
an exact treatment of advance of periapsis whereas ELL1 uses a linear approximation. 

Reference:
    [Susobhanan+ 2018](https://doi.org/10.1093/mnras/sty2177),
"""
struct BinaryELL1k <: BinaryELL1Base
    use_fbx::Bool
end

shapiro_delay_params(::BinaryELL1k, params::NamedTuple) = params.M2, params.SINI

function ELL1State(ell1::BinaryELL1k, toa::TOA, toacorr::TOACorrection, params::NamedTuple)
    Δt = corrected_toa_value(toa, toacorr, Float64) - params.TASC
    a1 = params.A1 + Δt * params.A1DOT

    ωdot = params.OMDOT
    lnedot = params.LNEDOT
    ϵ10 = params.EPS1
    ϵ20 = params.EPS2
    sin_ωdot_Δt, cos_ωdot_Δt = sincos(ωdot * Δt)
    ϵ1 = (1 + lnedot * Δt) * (ϵ10 * cos_ωdot_Δt + ϵ20 * sin_ωdot_Δt)
    ϵ2 = (1 + lnedot * Δt) * (ϵ20 * cos_ωdot_Δt - ϵ10 * sin_ωdot_Δt)

    Φ = mean_anomaly(Δt, params, ell1.use_fbx)
    Φ_trigs = sincos(Φ), sincos(2Φ), sincos(3Φ), sincos(4Φ)

    n = mean_motion(Δt, params, ell1.use_fbx)

    m2, sini = shapiro_delay_params(ell1, params)

    return ELL1State(Φ_trigs, n, a1, ϵ1, ϵ2, m2, sini)
end

function rømer_delay(ell1k::BinaryELL1k, state::ELL1State)::GQ
    return invoke(rømer_delay, Tuple{BinaryELL1Base,ELL1State}, ell1k, state) -
           (3 // 2) * state.a1 * state.ϵ1
end
