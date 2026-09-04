import astropy.units as u
import numpy as np
from pint import DMconst
from pint.toa import TOAs

from .vela import jl, to_jldd, vl

day_to_s = 86400


def pint_nbtoas_to_vela(toas: TOAs, epoch_mjd: float):
    """Construct a `Vela.TOA` object from a `PINT` `TOAs` object and an index."""

    assert (
        toas.planets
    ), "Planerary ephemeris not found in `TOAs` object. Use `planets=True` while reading in the TOAs."
    assert (
        toas.get_pulse_numbers() is not None
    ), "Pulse numbers not found in `TOAs` object. Call `toas.compute_pulse_numbers()`."

    if toas.tzr:
        assert len(toas) == 1, "TZR TOA must have length 1."

    tdb_lds = (toas.table["tdbld"].value - epoch_mjd) * day_to_s

    pulse_numbers = toas.table["pulse_number"].value
    delta_pulse_numbers = toas.table["delta_pulse_number"].value

    phases_ = pulse_numbers - delta_pulse_numbers

    phase1s, phase2s = np.modf(phases_)

    errs_ = toas.get_errors().to_value(u.s)
    freqs_ = toas.get_freqs().to_value(u.Hz)

    ssb_obs_poss_ = toas.table["ssb_obs_pos"].quantity.to_value(u.lightsecond)
    ssb_obs_vels_ = toas.table["ssb_obs_vel"].quantity.to_value(u.lightsecond / u.s)
    obs_sun_poss_ = toas.table["obs_sun_pos"].quantity.to_value(u.lightsecond)
    obs_jupiter_poss_ = toas.table["obs_jupiter_pos"].quantity.to_value(u.lightsecond)
    obs_saturn_poss_ = toas.table["obs_saturn_pos"].quantity.to_value(u.lightsecond)
    obs_venus_poss_ = toas.table["obs_venus_pos"].quantity.to_value(u.lightsecond)
    obs_uranus_poss_ = toas.table["obs_uranus_pos"].quantity.to_value(u.lightsecond)
    obs_neptune_poss_ = toas.table["obs_neptune_pos"].quantity.to_value(u.lightsecond)

    vtoas = []
    for ii in range(len(toas)):
        tdb = vl.time(to_jldd(tdb_lds[ii]))
        phase = vl.dimensionless(vl.Double64(phase2s[ii], phase1s[ii]))
        err = vl.time(errs_[ii])
        freq = vl.frequency(freqs_[ii])

        ssb_obs_pos = jl.map(
            vl.distance,
            jl.Tuple(ssb_obs_poss_[ii, :]),
        )
        ssb_obs_vel = jl.map(
            vl.speed,
            jl.Tuple(ssb_obs_vels_[ii, :]),
        )
        obs_sun_pos = jl.map(
            vl.distance,
            jl.Tuple(obs_sun_poss_[ii, :]),
        )
        obs_jupiter_pos = jl.map(
            vl.distance,
            jl.Tuple(obs_jupiter_poss_[ii, :]),
        )
        obs_saturn_pos = jl.map(
            vl.distance,
            jl.Tuple(obs_saturn_poss_[ii, :]),
        )
        obs_venus_pos = jl.map(
            vl.distance,
            jl.Tuple(obs_venus_poss_[ii, :]),
        )
        obs_uranus_pos = jl.map(
            vl.distance,
            jl.Tuple(obs_uranus_poss_[ii, :]),
        )
        obs_neptune_pos = jl.map(
            vl.distance,
            jl.Tuple(obs_neptune_poss_[ii, :]),
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

        idx = ii + 1 if not toas.tzr else 0
        vtoas.append(vl.TOA(tdb, err, freq, phase, ephem, idx))

    return jl.Vector[vl.TOA](vtoas)


def pint_wbtoas_to_vela(toas: TOAs, epoch_mjd: float):
    """Construct a `Vela.WidebandTOA`s from a `PINT` `TOAs` object containing
    wideband data and an index."""

    assert (
        toas.is_wideband()
    ), "Expected a wideband `TOAs` object here. Make sure that all TOAs have the `-ppdm` and `-ppdme` flags."

    vtoas = pint_nbtoas_to_vela(toas, epoch_mjd)

    dms = (DMconst * toas.get_dms()).si.value
    dmerrs = (DMconst * toas.get_dm_errors()).si.value
    dminfos = [
        vl.DMInfo(vl.GQ[-1](dm), vl.GQ[-1](dmerr)) for dm, dmerr in zip(dms, dmerrs)
    ]

    return jl.Vector[vl.WidebandTOA](
        [vl.WidebandTOA(vtoa, dminfo) for vtoa, dminfo in zip(vtoas, dminfos)]
    )


def pint_toas_to_vela(toas: TOAs, epoch_mjd: float):
    """Construct a Julia `Vector` of `Vela.TOA` or `Vela.WidebandTOA` objects from a
    `PINT` TOAs object."""

    return (
        pint_nbtoas_to_vela(toas, epoch_mjd)
        if not toas.is_wideband()
        else pint_wbtoas_to_vela(toas, epoch_mjd)
    )
