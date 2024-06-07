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

struct SolarSystem <: DelayComponent
    ecliptic_coordinates::Bool
    planet_shapiro::Bool
end

read_params_from_dict(::SolarSystem, params::Dict) = (
    POSEPOCH = params[:POSEPOCH][1],
    LONG = params[:LONG][1],
    LAT = params[:LAT][1],
    PMLONG = params[:PMLONG][1],
    PMLAT = params[:PMLAT][1],
    PX = params[:PX][1],
)

"""Compute the Shapiro delay due to a solar system object."""
function solar_system_shapiro_delay(M::GQ, obs_obj_pos::Tuple, obj_psr_pos::Tuple)
    Lhat = obj_psr_pos
    rvec = obs_obj_pos
    r = sqrt(dot(rvec, rvec))
    Lhat_dot_rvec = dot(Lhat, rvec)
    return -2 * M * log((r - Lhat_dot_rvec) / AU)
end

function correct_toa(ss::SolarSystem, ctoa::CorrectedTOA, params::NamedTuple)
    if ctoa.barycentered
        return correct_toa(ctoa)
    end

    long0 = params.LONG
    lat0 = params.LAT

    # TODO: Do this properly.
    if params.PMLONG == dimensionless(0.0) && params.PMLAT == dimensionless(0.0)
        long, lat = long0, lat0
    else
        dt = corrected_toa_value(ctoa) - params.POSEPOCH
        long = long0 + params.PMLONG * dt / cos(params.LAT)
        lat = lat0 + params.PMLAT * dt
    end

    Lhat = (cos(long) * cos(lat), sin(long) * cos(lat), sin(lat))
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

    if params.PX != 0
        Rvec_sqr = dot(Rvec, Rvec)
        Rperp_sqr = Rvec_sqr - Lhat_dot_Rvec^2

        # Parallax delay
        delay += 0.5 * params.PX * Rperp_sqr
    end

    # Sun Shapiro delay
    delay += solar_system_shapiro_delay(M_SUN, toa.ephem.obs_sun_pos, Lhat)

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

    return correct_toa(ctoa; delay = delay, doppler = doppler, barycentered = true)
end
