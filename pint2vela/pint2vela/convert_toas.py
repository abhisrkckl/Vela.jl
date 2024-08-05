import numpy as np

from pint.toa import TOAs

from .vela import jl, vl, to_jldd

day_to_s = 86400


def pint_toa_to_vela(toas: TOAs, idx: int, tzr: bool = False):
    assert toas.planets
    assert toas.get_pulse_numbers() is not None

    tdb_ld = toas.table["tdbld"].value[idx] * day_to_s
    # tdb_ld1, tdb_ld2 = np.modf(tdb_ld)
    # tdb = vl.time(vl.Double64(tdb_ld2, tdb_ld1))
    tdb = vl.time(to_jldd(tdb_ld))

    phase = (
        toas.table["pulse_number"].value[idx]
        - toas.table["delta_pulse_number"].value[idx]
    )
    phase1, phase2 = np.modf(phase)
    phase = vl.dimensionless(vl.Double64(phase2, phase1))

    err = vl.time(toas.get_errors()[idx].si.value)
    freq = vl.frequency(toas.get_freqs()[idx].si.value)

    ssb_obs_pos = jl.map(
        vl.distance,
        jl.Tuple(toas.table["ssb_obs_pos"].quantity[idx].to_value("lightsecond")),
    )
    ssb_obs_vel = jl.map(
        vl.speed,
        jl.Tuple(toas.table["ssb_obs_vel"].quantity[idx].to_value("lightsecond/s")),
    )
    obs_sun_pos = jl.map(
        vl.distance,
        jl.Tuple(toas.table["obs_sun_pos"].quantity[idx].to_value("lightsecond")),
    )

    obs_jupiter_pos = jl.map(
        vl.distance,
        jl.Tuple(toas.table["obs_jupiter_pos"].quantity[idx].to_value("lightsecond")),
    )
    obs_saturn_pos = jl.map(
        vl.distance,
        jl.Tuple(toas.table["obs_saturn_pos"].quantity[idx].to_value("lightsecond")),
    )
    obs_venus_pos = jl.map(
        vl.distance,
        jl.Tuple(toas.table["obs_venus_pos"].quantity[idx].to_value("lightsecond")),
    )
    obs_uranus_pos = jl.map(
        vl.distance,
        jl.Tuple(toas.table["obs_uranus_pos"].quantity[idx].to_value("lightsecond")),
    )
    obs_neptune_pos = jl.map(
        vl.distance,
        jl.Tuple(toas.table["obs_neptune_pos"].quantity[idx].to_value("lightsecond")),
    )
    obs_earth_pos = jl.map(
        vl.distance,
        jl.Tuple(toas.table["obs_earth_pos"].quantity[idx].to_value("lightsecond")),
    )

    ephem = vl.SolarSystemEphemeris(
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

    return vl.TOA(tdb, err, freq, phase, False, tzr, ephem, idx + 1)


def pint_toas_to_vela(toas: TOAs):
    return jl.Vector[vl.TOA]([pint_toa_to_vela(toas, idx) for idx in range(len(toas))])
