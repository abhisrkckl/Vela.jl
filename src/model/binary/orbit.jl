"""Compute the mean anomaly given the orbital period and its first derivative
OR the orbital frequency and its derivatives."""
mean_anomaly(Δt::GQ, params::NamedTuple, use_fbx::Bool) =
    2π * (
        use_fbx ? taylor_horner_integral(Δt, params.FB, dimensionless(0.0)) :
        (Δt / params.PB - 0.5 * params.PBDOT * (Δt / params.PB)^2)
    )

"""Compute the mean motion given the orbital period and its first derivative
OR the orbital frequency and its derivatives."""
mean_motion(Δt::GQ, params::NamedTuple, use_fbx::Bool) =
    2π * (use_fbx ? taylor_horner(Δt, params.FB) : 1 / (params.PB + params.PBDOT * Δt))

SinCos = NTuple{2,GQ{Float64}}

"""Solves the Kepler equation using Mikkola's method."""
function mikkola(l, e)
    if iszero(e) || iszero(l)
        return l
    end

    @assert (zero(e) < e < oneunit(e)) "The eccentricity of an ellipse must lie within [0,1)."

    sgn = sign(l)
    l *= sgn # l > 0

    ncycles = floor(l / (2π))
    l -= 2π * ncycles # 0 ≤ l < 2π

    flag = l > π
    if flag
        l = 2π - l # 0 ≤ l ≤ π
    end

    α = (1 - e) / (4 * e + 0.5)
    α3 = α * α * α
    β = (l / 2) / (4 * e + 0.5)
    β2 = β * β

    z = (β > 0) ? cbrt(β + sqrt(α3 + β2)) : cbrt(β - sqrt(α3 + β2))

    s = (z - α / z)
    s5 = s * s * s * s * s
    w = (s - (0.078 * s5) / (1 + e))
    w3 = w * w * w

    E0 = (l + e * (3 * w - 4 * w3))
    u = E0

    esu, ecu = e .* sincos(u)

    fu = u - esu - l
    f1u = 1 - ecu
    f2u = esu
    f3u = ecu
    f4u = -esu

    u1 = -fu / f1u
    u2 = -fu / (f1u + f2u * u1 / 2)
    u3 = -fu / (f1u + f2u * u2 / 2 + f3u * (u2 * u2) / 6.0)
    u4 = -fu / (f1u + f2u * u3 / 2 + f3u * (u3 * u3) / 6.0 + f4u * (u3 * u3 * u3) / 24.0)
    xi = (E0 + u4)

    sol = flag ? (2π - xi) : xi

    u = sgn * (sol + ncycles * 2π)

    return u
end

true_anomaly(u, e) = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(u / 2))
