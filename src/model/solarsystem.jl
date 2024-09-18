export SolarSystem, correct_toa

const OBL = dimensionless(0.4090926006005829)

const AU = distance(499.00478383615643)

const M_SUN = mass(4.92549094830932e-06)
const M_MERCURY = mass(8.176988758067135e-13)
const M_VENUS = mass(1.205680558494223e-11)
const M_EARTH = mass(1.4975623478139775e-11)
const M_MARS = mass(1.5895305231436198e-12)
const M_JUPITER = mass(4.702819050227708e-09)
const M_SATURN = mass(1.408128810019423e-09)
const M_URANUS = mass(2.1505895513637613e-10)
const M_NEPTUNE = mass(2.5373119991867603e-10)

"""Solar system Rømer delay, Parallax delay, and Shapiro delay.

Corresponds to `AstrometryEquatorial`, `AstrometryEcliptic`, 
and `SolarSystemShapiro` in `PINT`.

Reference:
    [Edwards+ 2006](http://doi.org/10.1111/j.1365-2966.2006.10870.x)
"""
struct SolarSystem <: DelayComponent
    ecliptic_coordinates::Bool
    planet_shapiro::Bool
end

"""Compute the Shapiro delay due to a solar system object."""
function solar_system_shapiro_delay(M::GQ, obs_obj_pos::Tuple, obj_psr_pos::Tuple)
    Lhat = obj_psr_pos
    rvec = obs_obj_pos
    r = sqrt(dot(rvec, rvec))
    Lhat_dot_rvec = dot(Lhat, rvec)
    return -2 * M * log((r - Lhat_dot_rvec) / AU)
end


"""Evaluate proper motion assuming uniform linear motion of the pulsar.

Reference:
    [Gaia Data Release Documentation](https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu3ast/sec_cu3ast_intro/ssec_cu3ast_intro_motion.html)
"""
function evaluate_proper_motion(α0, δ0, pmα, pmδ, dt)
    sinα0, cosα0 = sincos(α0)
    sinδ0, cosδ0 = sincos(δ0)

    x0 = (cosα0 * cosδ0, sinα0 * cosδ0, sinδ0)

    if iszero(pmα) && iszero(pmδ)
        return x0
    end

    xdot_α = (-sinα0 * pmα, cosα0 * pmα, zero(pmα))
    xdot_δ = (-cosα0 * sinδ0 * pmδ, -sinα0 * sinδ0 * pmδ, cosδ0 * pmδ)
    Δx = map((x, y) -> (x + y) * dt, xdot_α, xdot_δ)

    x1 = map(+, x0, Δx)
    x1_mag = sqrt(dot(x1, x1))
    x1hat = map(x -> x / x1_mag, x1)

    return x1hat
end

"""Update the `CorrectedTOA` object with solar system delays and Doppler factor."""
function correct_toa(ss::SolarSystem, ctoa::CorrectedTOA, params::NamedTuple)
    if is_barycentered(ctoa)
        return correct_toa(ctoa)
    end

    long0, lat0 =
        ss.ecliptic_coordinates ? (params.ELONG, params.ELAT) : (params.RAJ, params.DECJ)
    pmlong, pmlat =
        ss.ecliptic_coordinates ? (params.PMELONG, params.PMELAT) :
        (params.PMRA, params.PMDEC)
    px = params.PX
    posepoch = params.POSEPOCH

    dt = corrected_toa_value(ctoa) - posepoch
    Lhat = evaluate_proper_motion(long0, lat0, pmlong, pmlat, dt)

    Rvec = ctoa.toa.ephem.ssb_obs_pos

    if ss.ecliptic_coordinates
        x, y, z = Lhat
        x1 = x
        y1 = cos(OBL) * y - sin(OBL) * z
        z1 = sin(OBL) * y + cos(OBL) * z
        Lhat = (x1, y1, z1)
    end

    Lhat_dot_Rvec = dot(Lhat, Rvec)

    # Rømer delay
    delay = -Lhat_dot_Rvec

    if !iszero(px)
        Rvec_sqr = dot(Rvec, Rvec)
        Rperp_sqr = Rvec_sqr - Lhat_dot_Rvec * Lhat_dot_Rvec

        # Parallax delay
        delay += 0.5 * px * Rperp_sqr
    end

    # Sun Shapiro delay
    delay += solar_system_shapiro_delay(M_SUN, ctoa.toa.ephem.obs_sun_pos, Lhat)

    if ss.planet_shapiro
        masses = (M_JUPITER, M_SATURN, M_VENUS, M_URANUS, M_NEPTUNE)
        positions = (
            ctoa.toa.ephem.obs_jupiter_pos,
            ctoa.toa.ephem.obs_saturn_pos,
            ctoa.toa.ephem.obs_venus_pos,
            ctoa.toa.ephem.obs_uranus_pos,
            ctoa.toa.ephem.obs_neptune_pos,
        )

        for (M, obs_obj_pos) in zip(masses, positions)
            # Planet Shapiro delays
            delay += solar_system_shapiro_delay(M, obs_obj_pos, Lhat)
        end
    end

    vvec = ctoa.toa.ephem.ssb_obs_vel
    Lhat_dot_vvec = dot(Lhat, vvec)
    # Doppler factor due to solar system motion
    # doppler = -d_roemerdelay / d_t
    # This applies to the frequency as freq*(1-doppler)
    doppler = Lhat_dot_vvec

    return correct_toa(ctoa; delay = delay, doppler = doppler, ssb_psr_pos = Lhat)
end

function show(io::IO, ss::SolarSystem)
    coordstr = ss.ecliptic_coordinates ? "Ecliptic" : "Equatorial"
    print(io, "SolarSystem($coordstr, planet_shapiro=$(ss.planet_shapiro))")
end
