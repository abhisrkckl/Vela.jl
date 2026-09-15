export BinaryDDR

const DDR_R0 = distance(8.417380286943844e11)
const DDR_THETA0_OVER_C = dimensionless(0.0007338410094359345)
const DDR_RHO0 = GQ{-2}(4.517103049894966e-31)
const DDR_Z0 = distance(1.852688250978102e10)
const DDR_ZSUN = distance(2.05854250108678e9)
const ICRS_TO_GAL = (
    (-0.054875657712591633, -0.87343705195561583, -0.48383507361671546),
    (0.49410943719272682, -0.44482972122329512, 0.7469821839866676),
    (-0.8676661375596576, -0.19807633727300059, 0.45598381368730162),
)

"""Damour-Deruelle-Regular binary model.

This model uses Laplace-Lagrange coordinates at `TASC`, a regular
``q = \nu - M`` periastron advance, and an optional Cartesian viewing
projector at a fixed `TGEO`.

Reference: PINT `DDRmodel` / `ddr_kepler` (van Haasteren in prep.).
"""
struct BinaryDDR <: BinaryComponent
    use_fbx::Bool
    ecliptic_coordinates::Bool
    use_pk::Bool
    pbdot_kinematic::Bool
    use_geo::Bool
    use_kine::Bool

    function BinaryDDR(
        use_fbx::Bool,
        ecliptic_coordinates::Bool,
        use_pk::Bool,
        pbdot_kinematic::Bool,
        use_geo::Bool,
        use_kine::Bool,
    )
        if use_fbx
            @assert !pbdot_kinematic && !use_kine
        elseif use_geo
            @assert use_kine
        end
        new(use_fbx, ecliptic_coordinates, use_pk, pbdot_kinematic, use_geo, use_kine)
    end
end

struct DDRState
    x::GQ{1,Float64}
    n::GQ{-1,Float64}
    c::GQ{0,Float64}
    s::GQ{0,Float64}
    c_e::GQ{0,Float64}
    s_e::GQ{0,Float64}
    g_gamma::GQ{1,Float64}
    X::GQ{0,Float64}
    Y::GQ{0,Float64}
    dX::GQ{0,Float64}
    dY::GQ{0,Float64}
    d2X::GQ{0,Float64}
    d2Y::GQ{0,Float64}
    I::GQ{0,Float64}
    J::GQ{0,Float64}
    B_S::GQ{0,Float64}
    m2::GQ{1,Float64}
    valid::Bool
end

function _invalid_ddr_state()
    nan0 = dimensionless(NaN)
    return DDRState(
        time(NaN),
        frequency(NaN),
        nan0,
        nan0,
        nan0,
        nan0,
        time(NaN),
        nan0,
        nan0,
        nan0,
        nan0,
        nan0,
        nan0,
        nan0,
        nan0,
        nan0,
        mass(NaN),
        false,
    )
end

_ddr_isfinite(x::GQ) = isfinite(value(x))
_ddr_allfinite(xs::Tuple) = all(_ddr_isfinite, xs)

"""Reduce a secular longitude to the interval selected by nearest-orbit rounding."""
function reduce_longitude(λ::GQ{0,Float64})
    n_orb = round(value(λ) / (2π))
    return λ - n_orb * 2π, n_orb
end

"""Solve the regular Laplace-Lagrange Kepler equation."""
function solve_F(λ::GQ{0,Float64}, h::GQ{0,Float64}, k::GQ{0,Float64})
    λ_red, _ = reduce_longitude(λ)
    λv = value(λ_red)
    hv = value(h)
    kv = value(k)
    lo = λv - 1.0
    hi = λv + 1.0
    F = clamp(λv + kv * sin(λv) - hv * cos(λv), lo, hi)
    atol = 4 * eps(Float64) * max(1.0, abs(λv))
    converged = false

    for _ = 1:64
        sF, cF = sincos(F)
        R = F - kv * sF + hv * cF - λv
        if abs(R) <= atol
            converged = true
            break
        end
        D = 1.0 - kv * cF - hv * sF
        F_try = F - R / D
        good = isfinite(F_try) && lo <= F_try <= hi
        slo, clo = sincos(lo)
        R_lo = lo - kv * slo + hv * clo - λv
        go_hi = R_lo * R > 0
        lo_n = go_hi ? F : lo
        hi_n = go_hi ? hi : F
        F = good ? F_try : 0.5 * (lo_n + hi_n)
        lo = good ? lo : lo_n
        hi = good ? hi : hi_n
    end

    if !converged
        sF, cF = sincos(F)
        converged = abs(F - kv * sF + hv * cF - λv) <= atol
    end
    if !converged
        nan0 = dimensionless(NaN)
        return nan0, nan0, nan0, nan0
    end

    sF, cF = sincos(F)
    D = 1.0 - kv * cF - hv * sF
    c_e = kv * cF + hv * sF
    s_e = kv * sF - hv * cF
    return dimensionless(F), dimensionless(D), dimensionless(c_e), dimensionless(s_e)
end

"""Regular static orbital projections and their frozen-element derivatives."""
function static_XY(F::GQ{0,Float64}, h::GQ{0,Float64}, k::GQ{0,Float64})
    E2 = h * h + k * k
    η = sqrt(1 - E2)
    b = 1 / (1 + η)
    sF, cF = sincos(F)
    Y0 = (1 - b * k * k) * sF + b * h * k * cF - h
    X0 = (1 - b * h * h) * cF + b * h * k * sF - k
    dY0 = (1 - b * k * k) * cF - b * h * k * sF
    d2Y0 = -(1 - b * k * k) * sF - b * h * k * cF
    dX0 = -(1 - b * h * h) * sF + b * h * k * cF
    d2X0 = -(1 - b * h * h) * cF - b * h * k * sF
    return X0, Y0, dX0, dY0, d2X0, d2Y0
end

q_nu_minus_M(c_e, s_e, E2) = s_e + 2 * atan(s_e, 1 + sqrt(1 - E2) - c_e)

function q_at_tasc(h::GQ{0,Float64}, k::GQ{0,Float64})
    _, _, c_e, s_e = solve_F(zero(h), h, k)
    return q_nu_minus_M(c_e, s_e, h * h + k * k)
end

precession_delta(λ, q, q_star, κ) = κ * (λ + q - q_star)

function rotate_XY(X0, Y0, δ)
    sδ, cδ = sincos(δ)
    return X0 * cδ - Y0 * sδ, Y0 * cδ + X0 * sδ
end

sini_from_cosi(c) = sqrt((1 - c) * (1 + c))

function mean_longitude(Δt::GQ{1,Float64}, PB::GQ{1,Float64}, p::GQ{0,Float64})
    phase = Δt / PB
    λ = 2π * (phase - 0.5 * p * phase * phase)
    λdot = 2π * (1 / PB - p * Δt / (PB * PB))
    return λ, λdot
end

function pulsar_mass(n, x_star, mc, c)
    f = (n * n) * (x_star * x_star * x_star) / M_SUN
    s = sini_from_cosi(c)
    mp = (mc * s)^(3 / 2) / sqrt(f) - mc
    return mp, s
end

function g_gamma_gr(n, mp, mc)
    nTsun = n * M_SUN
    return M_SUN / nTsun^(1 / 3) * mc * (mp + 2 * mc) / (mp + mc)^(4 / 3)
end

kappa_gr(x_star, mc, s, E2) = 3 * (M_SUN / x_star) * mc * s / (1 - E2)

function pbdot_gw(n, mp, mc, E2)
    E4 = E2 * E2
    fe = (1 + (73 / 24) * E2 + (37 / 96) * E4) / (1 - E2)^(7 / 2)
    return -192 * π / 5 * (n * M_SUN)^(5 / 3) * fe * mp * mc / (mp + mc)^(1 / 3)
end

function sky_motion(α, δ, μα, μδ, dt)
    sα, cα = sincos(α)
    sδ, cδ = sincos(δ)
    n_p = (cα * cδ, sα * cδ, sδ)
    ndot_p = (-sα * μα - cα * sδ * μδ, cα * μα - sα * sδ * μδ, cδ * μδ)
    r = (n_p[1] + ndot_p[1] * dt, n_p[2] + ndot_p[2] * dt, n_p[3] + ndot_p[3] * dt)
    rnorm = sqrt(dot(r, r))
    n = (r[1] / rnorm, r[2] / rnorm, r[3] / rnorm)
    radial = dot(ndot_p, n)
    ndot = (
        (ndot_p[1] - radial * n[1]) / rnorm,
        (ndot_p[2] - radial * n[2]) / rnorm,
        (ndot_p[3] - radial * n[3]) / rnorm,
    )
    return n, ndot
end

function _native_astrometry(ddr::BinaryDDR, params::NamedTuple)
    return ddr.ecliptic_coordinates ?
           (params.ELONG, params.ELAT, params.PMELONG, params.PMELAT) :
           (params.RAJ, params.DECJ, params.PMRA, params.PMDEC)
end

function _has_ddr_astrometry(ddr::BinaryDDR, params::NamedTuple)
    common = haskey(params, :PX) && haskey(params, :TGEO) && haskey(params, :POSEPOCH)
    native =
        ddr.ecliptic_coordinates ?
        (
            haskey(params, :ELONG) &&
            haskey(params, :ELAT) &&
            haskey(params, :PMELONG) &&
            haskey(params, :PMELAT)
        ) :
        (
            haskey(params, :RAJ) &&
            haskey(params, :DECJ) &&
            haskey(params, :PMRA) &&
            haskey(params, :PMDEC)
        )
    return common && native
end

function _ddr_sky_basis(ddr::BinaryDDR, params::NamedTuple)
    α, δ, μα, μδ = _native_astrometry(ddr, params)
    n0, ndot = sky_motion(α, δ, μα, μδ, params.TGEO - params.POSEPOCH)
    ρ = sqrt(n0[1] * n0[1] + n0[2] * n0[2])
    east = (-n0[2] / ρ, n0[1] / ρ, zero(n0[1]))
    north = (
        n0[2] * east[3] - n0[3] * east[2],
        n0[3] * east[1] - n0[1] * east[3],
        n0[1] * east[2] - n0[2] * east[1],
    )
    return n0, ndot, east, north
end

function _ecliptic_to_icrs(v)
    sϵ, cϵ = sincos(OBL)
    x, y, z = v
    return x, cϵ * y - sϵ * z, sϵ * y + cϵ * z
end

function _icrs_to_ecliptic(v)
    sϵ, cϵ = sincos(OBL)
    x, y, z = v
    return x, cϵ * y + sϵ * z, -sϵ * y + cϵ * z
end

IJ_from_v(v_I, v_J, Ω) = (-v_I * sin(Ω) + v_J * cos(Ω), v_I * cos(Ω) + v_J * sin(Ω))

v_from_mu_parallax(μ_I, μ_J, px, d_I, d_J, Δt_K) =
    (μ_I * Δt_K - px * d_I, μ_J * Δt_K - px * d_J)

function _ddr_geometry(ddr::BinaryDDR, toa::TOA, tcorr::GQ{1,Float64}, params::NamedTuple)
    _, ndot, east, north = _ddr_sky_basis(ddr, params)
    observer =
        ddr.ecliptic_coordinates ? _icrs_to_ecliptic(toa.ephem.ssb_obs_pos) :
        toa.ephem.ssb_obs_pos
    μ_I = dot(east, ndot)
    μ_J = dot(north, ndot)
    d_I = dot(observer, east)
    d_J = dot(observer, north)
    v_I, v_J = v_from_mu_parallax(μ_I, μ_J, params.PX, d_I, d_J, tcorr - params.TGEO)
    return IJ_from_v(v_I, v_J, params.KOM)
end

function _galactic_direction(ddr::BinaryDDR, params::NamedTuple)
    n0, _, _, _ = _ddr_sky_basis(ddr, params)
    n_icrs = ddr.ecliptic_coordinates ? _ecliptic_to_icrs(n0) : n0
    n_gal = (
        dot(ICRS_TO_GAL[1], n_icrs),
        dot(ICRS_TO_GAL[2], n_icrs),
        dot(ICRS_TO_GAL[3], n_icrs),
    )
    l = atan(n_gal[2], n_gal[1])
    b = atan(n_gal[3], sqrt(n_gal[1] * n_gal[1] + n_gal[2] * n_gal[2]))
    return l, b
end

function _galactic_acceleration_los(distance_to_pulsar, l, b)
    sb, cb = sincos(b)
    sl, cl = sincos(l)
    β = (distance_to_pulsar / DDR_R0) * cb - cl
    denominator = sl * sl + β * β
    a_planar =
        -cb * (DDR_THETA0_OVER_C * DDR_THETA0_OVER_C / DDR_R0) * (cl + β / denominator)

    a_z(z) = -4π * DDR_RHO0 * DDR_Z0 * z / sqrt(z * z + DDR_Z0 * DDR_Z0)
    z_psr = DDR_ZSUN + distance_to_pulsar * sb
    a_vertical = (a_z(z_psr) - a_z(DDR_ZSUN)) * sb
    return a_planar + a_vertical
end

function _compose_p(ddr::BinaryDDR, n, mp, mc, E2, params::NamedTuple)
    p_gw = dimensionless(0.0)
    if ddr.pbdot_kinematic
        p_gw = pbdot_gw(n, mp, mc, E2)
        p = p_gw + params.XPBDOT
    else
        p = params.PBDOT
    end

    p_shk = dimensionless(0.0)
    p_gal = dimensionless(0.0)
    if ddr.use_kine
        _, _, μα, μδ = _native_astrometry(ddr, params)
        μ2 = μα * μα + μδ * μδ
        p_shk = params.PB * μ2 / params.PX
        l, b = _galactic_direction(ddr, params)
        p_gal = params.PB * _galactic_acceleration_los(1 / params.PX, l, b)
        p += p_shk + p_gal
    end
    return p, p_shk, p_gal, p_gw
end

function _roemer(x, c, s, I, J, X, Y)
    a = x / s
    Z = sqrt(1 + I * I + J * J)
    return (x * Y + a * c * I * Y + a * J * X) / Z
end

function _shapiro_B_S_squared_norm(c_e, X, Y, c, s, I, J, Ω)
    ρ = 1 - c_e
    sΩ, cΩ = sincos(Ω)
    Rx = (X * cΩ - Y * c * sΩ) / ρ
    Ry = (X * sΩ + Y * c * cΩ) / ρ
    Rz = Y * s / ρ
    Z = sqrt(1 + I * I + J * J)
    v_I = -I * sΩ + J * cΩ
    v_J = I * cΩ + J * sΩ
    nax = v_I / Z
    nay = v_J / Z
    naz = 1 / Z
    dx = nax - Rx
    dy = nay - Ry
    dz = naz - Rz
    return 0.5 * ρ * (dx * dx + dy * dy + dz * dz)
end

function DDRState(ddr::BinaryDDR, toa::TOA, toacorr::TOACorrection, params::NamedTuple)
    tcorr = corrected_toa_value(toa, toacorr, Float64)
    if !haskey(params, :TASC) ||
       !haskey(params, :EPS1) ||
       !haskey(params, :EPS2) ||
       !haskey(params, :COSI) ||
       !haskey(params, :A1) ||
       !haskey(params, :A1DOT) ||
       !haskey(params, :M2)
        return _invalid_ddr_state()
    end

    Δt = tcorr - params.TASC
    h = params.EPS1
    k = params.EPS2
    c = params.COSI
    E2 = h * h + k * k
    if !_ddr_allfinite((h, k, c, Δt, E2)) || abs(value(c)) >= 1 || value(E2) > 0.99^2
        return _invalid_ddr_state()
    end

    s = sini_from_cosi(c)
    x_star = params.A1
    m2 = params.M2
    mc = m2 / M_SUN
    x = x_star + Δt * params.A1DOT
    if !_ddr_allfinite((s, x_star, m2, mc, x, params.A1DOT))
        return _invalid_ddr_state()
    end

    if ddr.use_fbx
        if !haskey(params, :FB) || isempty(params.FB)
            return _invalid_ddr_state()
        end
        f0 = params.FB[1]
        if !_ddr_isfinite(f0) || value(f0) <= 0
            return _invalid_ddr_state()
        end
        n = 2π * f0
    else
        if !haskey(params, :PB) || !_ddr_isfinite(params.PB) || value(params.PB) <= 0
            return _invalid_ddr_state()
        end
        n = 2π / params.PB
    end

    need_mass = ddr.use_pk || ddr.pbdot_kinematic
    if need_mass
        if value(x_star) <= 0 || value(x) <= 0 || value(m2) <= 0
            return _invalid_ddr_state()
        end
        mp, s = pulsar_mass(n, x_star, mc, c)
        if !_ddr_allfinite((mp, s)) || value(mp) <= 0
            return _invalid_ddr_state()
        end
    else
        if value(x_star) < 0 || value(x) < 0 || value(m2) < 0
            return _invalid_ddr_state()
        end
        mp = dimensionless(0.0)
    end

    if ddr.use_pk
        g_gamma = g_gamma_gr(n, mp, mc)
        κ = kappa_gr(x_star, mc, s, E2)
    else
        if !haskey(params, :GGAMMA) ||
           !haskey(params, :OMDOT) ||
           !_ddr_allfinite((params.GGAMMA, params.OMDOT))
            return _invalid_ddr_state()
        end
        g_gamma = params.GGAMMA
        κ = params.OMDOT / n
    end
    if !_ddr_allfinite((g_gamma, κ))
        return _invalid_ddr_state()
    end

    if ddr.use_fbx
        if length(params.FB) == 1
            λ = 2π * Δt * params.FB[1]
            λdot = 2π * params.FB[1]
        else
            λ = 2π * taylor_horner_integral(Δt, params.FB, dimensionless(0.0))
            λdot = 2π * taylor_horner(Δt, params.FB)
        end
    else
        if ddr.pbdot_kinematic
            if !haskey(params, :XPBDOT) || !_ddr_isfinite(params.XPBDOT)
                return _invalid_ddr_state()
            end
        elseif !haskey(params, :PBDOT) || !_ddr_isfinite(params.PBDOT)
            return _invalid_ddr_state()
        end
        if ddr.use_kine
            if !_has_ddr_astrometry(ddr, params) ||
               !_ddr_isfinite(params.PX) ||
               value(params.PX) <= 0
                return _invalid_ddr_state()
            end
            astrometry = _native_astrometry(ddr, params)
            if !_ddr_allfinite((astrometry..., params.TGEO, params.POSEPOCH))
                return _invalid_ddr_state()
            end
        end
        p, _, _, _ = _compose_p(ddr, n, mp, mc, E2, params)
        if !_ddr_isfinite(p)
            return _invalid_ddr_state()
        end
        λ, λdot = mean_longitude(Δt, params.PB, p)
    end
    if !_ddr_allfinite((λ, λdot)) || value(λdot) <= 0
        return _invalid_ddr_state()
    end

    F, _, c_e, s_e = solve_F(λ, h, k)
    if !_ddr_allfinite((F, c_e, s_e))
        return _invalid_ddr_state()
    end
    X0, Y0, dX0, dY0, d2X0, d2Y0 = static_XY(F, h, k)
    q = q_nu_minus_M(c_e, s_e, E2)
    q_star = q_at_tasc(h, k)
    δ = precession_delta(λ, q, q_star, κ)
    X, Y = rotate_XY(X0, Y0, δ)
    dX, dY = rotate_XY(dX0, dY0, δ)
    d2X, d2Y = rotate_XY(d2X0, d2Y0, δ)
    if !_ddr_allfinite((X, Y, dX, dY, d2X, d2Y, q, q_star, δ))
        return _invalid_ddr_state()
    end

    if ddr.use_geo
        if !_has_ddr_astrometry(ddr, params) ||
           !haskey(params, :KOM) ||
           !_ddr_allfinite((params.PX, params.KOM, params.TGEO, params.POSEPOCH)) ||
           value(params.PX) <= 0
            return _invalid_ddr_state()
        end
        astrometry = _native_astrometry(ddr, params)
        if !_ddr_allfinite(astrometry)
            return _invalid_ddr_state()
        end
        I, J = _ddr_geometry(ddr, toa, tcorr, params)
        if !_ddr_allfinite((I, J))
            return _invalid_ddr_state()
        end
        B_S = _shapiro_B_S_squared_norm(c_e, X, Y, c, s, I, J, params.KOM)
    else
        I = dimensionless(0.0)
        J = dimensionless(0.0)
        B_S = _shapiro_B_S_squared_norm(c_e, X, Y, c, s, I, J, dimensionless(0.0))
    end
    if !_ddr_isfinite(B_S) || value(B_S) <= 0
        return _invalid_ddr_state()
    end

    return DDRState(
        x,
        n,
        c,
        s,
        c_e,
        s_e,
        g_gamma,
        X,
        Y,
        dX,
        dY,
        d2X,
        d2Y,
        I,
        J,
        B_S,
        m2,
        true,
    )
end

function rømer_einstein_delay(::BinaryDDR, st::DDRState)::GQ
    Δ_rom = _roemer(st.x, st.c, st.s, st.I, st.J, st.X, st.Y)
    return Δ_rom + st.g_gamma * st.s_e
end

function d_rømer_einstein_delay_d_F(::BinaryDDR, st::DDRState)::GQ
    Δ_romp = _roemer(st.x, st.c, st.s, st.I, st.J, st.dX, st.dY)
    return Δ_romp + st.g_gamma * st.c_e
end

function d2_rømer_einstein_delay_d_F2(::BinaryDDR, st::DDRState)::GQ
    Δ_rompp = _roemer(st.x, st.c, st.s, st.I, st.J, st.d2X, st.d2Y)
    return Δ_rompp - st.g_gamma * st.s_e
end

shapiro_delay(::BinaryDDR, st::DDRState)::GQ = -2 * st.m2 * log(st.B_S)

function correct_toa(ddr::BinaryDDR, toa::TOA, toacorr::TOACorrection, params::NamedTuple)
    st = DDRState(ddr, toa, toacorr, params)
    if !st.valid
        return correct_toa_delay(toacorr; delay = time(NaN), doppler = dimensionless(NaN))
    end

    d = rømer_einstein_delay(ddr, st)
    dp = d_rømer_einstein_delay_d_F(ddr, st)
    dpp = d2_rømer_einstein_delay_d_F2(ddr, st)
    nhat = st.n / (1 - st.c_e)
    d_inv =
        d * (
            1 - nhat * dp + (nhat * dp) * (nhat * dp) + 0.5 * (nhat * nhat) * d * dpp -
            0.5 * st.s_e / (1 - st.c_e) * (nhat * nhat) * d * dp
        )
    delay = d_inv + shapiro_delay(ddr, st)
    doppler = nhat * dp
    return correct_toa_delay(toacorr; delay = delay, doppler = doppler)
end
