export SolarSystemEphemeris

"""Struct containing solar system ephemeris vectors corresponding to a TOA.

These are computed using PINT.

References:
    [Hobbs+ 2006](http://doi.org/10.1111/j.1365-2966.2006.10302.x)
    [Luo+ 2021](http://doi.org/10.3847/1538-4357/abe62f)
"""
struct SolarSystemEphemeris
    ssb_obs_pos::NTuple{3,GQ{1,Float64}}
    ssb_obs_vel::NTuple{3,GQ{0,Float64}}
    obs_sun_pos::NTuple{3,GQ{1,Float64}}
    obs_jupiter_pos::NTuple{3,GQ{1,Float64}}
    obs_saturn_pos::NTuple{3,GQ{1,Float64}}
    obs_venus_pos::NTuple{3,GQ{1,Float64}}
    obs_uranus_pos::NTuple{3,GQ{1,Float64}}
    obs_neptune_pos::NTuple{3,GQ{1,Float64}}
    obs_earth_pos::NTuple{3,GQ{1,Float64}}

    function SolarSystemEphemeris(
        ssb_obs_pos,
        ssb_obs_vel,
        obs_sun_pos,
        obs_jupiter_pos,
        obs_saturn_pos,
        obs_venus_pos,
        obs_uranus_pos,
        obs_neptune_pos,
        obs_earth_pos,
    )
        @assert dot(ssb_obs_vel, ssb_obs_vel) < 1 "Magnitude of ssb_obs_vel should be less than 1."

        return new(
            ssb_obs_pos,
            ssb_obs_vel,
            obs_sun_pos,
            obs_jupiter_pos,
            obs_saturn_pos,
            obs_venus_pos,
            obs_uranus_pos,
            obs_neptune_pos,
            obs_earth_pos,
        )
    end
end
