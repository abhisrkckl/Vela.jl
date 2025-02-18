import numpy as np
from pint import DMconst
from pint.toa import TOAs

from .vela import jl, to_jldd, vl

day_to_s = 86400


def pint_toa_to_vela(toas: TOAs, idx: int, epoch_mjd: float):
    """Construct a `Vela.TOA` object from a `PINT` `TOAs` object and an index."""

    assert (
        toas.planets
    ), "Planerary ephemeris not found in `TOAs` object. Use `planets=True` while reading in the TOAs."
    assert (
        toas.get_pulse_numbers() is not None
    ), "Pulse numbers not found in `TOAs` object. Call `toas.compute_pulse_numbers()`."

    tdb_ld = (toas.table["tdbld"].value[idx] - epoch_mjd) * day_to_s
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

    ephem = vl.SolarSystemEphemeris(
        ssb_obs_pos,
        ssb_obs_vel,
        obs_sun_pos,
        obs_jupiter_pos,
        obs_saturn_pos,
        obs_venus_pos,
        obs_uranus_pos,
        obs_neptune_pos,
    )

    return vl.TOA(tdb, err, freq, phase, ephem, idx + 1)


def pint_wbtoa_to_vela(toas: TOAs, idx: int, epoch_mjd: float):
    """Construct a `Vela.WidebandTOA`s from a `PINT` `TOAs` object containing
    wideband data and an index."""

    assert (
        toas.is_wideband()
    ), "Expected a wideband `TOAs` object here. Make sure that all TOAs have the `-ppdm` and `-ppdme` flags."

    vtoa = pint_toa_to_vela(toas, idx, epoch_mjd)

    dm = vl.GQ[-1]((DMconst * toas.get_dms()[idx]).si.value)
    dmerr = vl.GQ[-1]((DMconst * toas.get_dm_errors()[idx]).si.value)
    dminfo = vl.DMInfo(dm, dmerr)

    return vl.WidebandTOA(vtoa, dminfo)


def pint_toas_to_vela(toas: TOAs, epoch_mjd: float):
    """Construct a Julia `Vector` of `Vela.TOA` or `Vela.WidebandTOA` objects from a
    `PINT` TOAs object."""

    p2v_toas = pint_wbtoa_to_vela if toas.is_wideband() else pint_toa_to_vela
    TOAType = vl.WidebandTOA if toas.is_wideband() else vl.TOA
    return jl.Vector[TOAType](
        [p2v_toas(toas, idx, epoch_mjd) for idx in range(len(toas))]
    )
