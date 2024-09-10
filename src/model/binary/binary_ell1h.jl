export BinaryELL1H

"""A binary model representing a nearly circular orbit with orthometric 
parametrization of the Shapiro delay. The Shapiro delay computed in this model
does not include the Fourier components that are fully covariant with the Rømer
delay.

Reference:
    [Lange+ 2001](http://doi.org/10.1046/j.1365-8711.2001.04606.x)
    [Freire & Wex 2010](http://doi.org/10.1111/j.1365-2966.2010.17319.x)
"""
struct BinaryELL1H <: BinaryELL1Base
    use_fbx::Bool
end

function shapiro_delay_params(::BinaryELL1H, params::NamedTuple)
    h3 = params.H3
    Ϛ = params.STIGMA
    m2 = h3 / Ϛ^3
    sini = 2 * Ϛ / (1 + Ϛ^2)
    return m2, sini
end

function shapiro_delay(ell1h::BinaryELL1H, state::ELL1State)
    ΔS_total = invoke(shapiro_delay, Tuple{BinaryELL1Base,ELL1State}, ell1h, state)

    (sinΦ, _), (_, cos2Φ), (_, _), (_, _) = state.Φ_trigs

    m2 = state.m2
    s = state.sini
    c̄ = sqrt(1 - s * s)
    Ϛ = s / (1 + c̄)
    a0 = -log(1 + Ϛ * Ϛ)
    b1 = -2 * Ϛ
    a2 = Ϛ * Ϛ
    ΔS_012 = -2 * m2 * (a0 + b1 * sinΦ + a2 * cos2Φ)

    return ΔS_total - ΔS_012
end
