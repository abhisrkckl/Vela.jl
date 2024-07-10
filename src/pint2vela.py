#!/usr/bin/env python

import json
import sys
from typing import Any, Dict, List, Tuple

import h5py as h5
import numpy as np
from astropy import units as u
from astropy.table import Table
from pint.logging import setup as setup_log
from pint.models import PhaseOffset, TimingModel, get_model_and_toas
from pint.models.parameter import (
    AngleParameter,
    MJDParameter,
    floatParameter,
)
from pint.toa import TOAs

from juliacall import Main as jl
jl.seval("import Vela")
vl = jl.Vela

day_to_s = 86400

def pint_toa_to_vela(toas: TOAs, idx: int, tzr: bool = False):
    assert toas.planets
    assert toas.get_pulse_numbers() is not None

    tdb = toas.table["tdbld"].value[idx] * day_to_s
    tdb = vl.time(jl.parse(vl.Float128, str(tdb)))

    phase = toas.table["pulse_number"].value[idx] - toas.table["delta_pulse_number"].value[idx]
    phase = vl.dimensionless(jl.parse(vl.Float128, str(phase)))

    err = vl.time(toas.get_errors()[idx].si.value)
    freq = vl.frequency(toas.get_freqs()[idx].si.value)

    ssb_obs_pos = jl.map(vl.distance, jl.Tuple(toas.table["ssb_obs_pos"].quantity[idx].to_value("lightsecond")))
    ssb_obs_vel = jl.map(vl.speed, jl.Tuple(toas.table["ssb_obs_vel"].quantity[idx].to_value("lightsecond/s")))
    obs_sun_pos = jl.map(vl.distance, jl.Tuple(toas.table["obs_sun_pos"].quantity[idx].to_value("lightsecond")))

    obs_jupiter_pos = jl.map(vl.distance, jl.Tuple(toas.table["obs_jupiter_pos"].quantity[idx].to_value("lightsecond")))
    obs_saturn_pos = jl.map(vl.distance, jl.Tuple(toas.table["obs_saturn_pos"].quantity[idx].to_value("lightsecond")))
    obs_venus_pos = jl.map(vl.distance, jl.Tuple(toas.table["obs_venus_pos"].quantity[idx].to_value("lightsecond")))
    obs_uranus_pos = jl.map(vl.distance, jl.Tuple(toas.table["obs_uranus_pos"].quantity[idx].to_value("lightsecond")))
    obs_neptune_pos = jl.map(vl.distance, jl.Tuple(toas.table["obs_neptune_pos"].quantity[idx].to_value("lightsecond")))
    obs_earth_pos = jl.map(vl.distance, jl.Tuple(toas.table["obs_earth_pos"].quantity[idx].to_value("lightsecond")))

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

    return vl.TOA(tdb, err, freq, phase, False, tzr, ephem)


def pint_toas_to_vela(toas: TOAs):
    return jl.Vector[vl.TOA](
        [pint_toa_to_vela(toas, idx) for idx in range(len(toas))]
    )


def fix_params(model: TimingModel) -> None:
    assert model.PEPOCH.value is not None

    for param in model.params:
        if (
            param.endswith("EPOCH")
            and isinstance(model[param], MJDParameter)
            and model[param].value is None
        ):
            model[param].quantity = model.PEPOCH.quantity

    if "PhaseOffset" not in model.components:
        model.add_component(PhaseOffset())

    model.PHOFF.frozen = False


def _get_multiparam_elements(
    model: TimingModel, multi_param_names: List[str]
) -> List[Dict[str, Any]]:
    elements = []
    for pxname in multi_param_names:
        pxparam = model[pxname]

        if pxparam.value is None:
            break

        scale_factor = pxparam.tcb2tdb_scale_factor
        value = (
            pxparam.value * day_to_s
            if isinstance(pxparam, MJDParameter)
            else (pxparam.quantity * scale_factor).si.value
        )
        dim = pxparam.effective_dimensionality

        original_units = str(pxparam.units)
        unit_conversion_factor = (pxparam.units * scale_factor / u.s**dim).to_value(
            u.dimensionless_unscaled, equivalencies=u.dimensionless_angles()
        )

        elements.append(
            {
                "name": pxparam.name,
                "default_value": float(value),
                "dimension": dim,
                "frozen": pxparam.frozen,
                "original_units": original_units,
                "unit_conversion_factor": float(unit_conversion_factor),
            }
        )

    return elements


def params_from_model(
    model: TimingModel,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    ignore_params = [
        "START",
        "FINISH",
        "RM",
        "CHI2",
        "CHI2R",
        "TRES",
        "DMRES",
        "TZRMJD",
        "TZRFRQ",
        "SWM",
    ]

    # Parameters that are defined as single parameters in PINT
    # but are really multi-parameters
    pseudo_single_params = ["DM", "CM"]

    assert all(psp not in ignore_params for psp in pseudo_single_params)
    assert all(
        psp not in model or not hasattr(model[psp], "prefix")
        for psp in pseudo_single_params
    )

    # Process single parameters
    single_params = []
    for param_name in model.params:
        param = model[param_name]

        # Skip things that are not single parameters
        if (
            param_name in ignore_params
            or param_name in pseudo_single_params
            or not isinstance(
                param,
                (
                    floatParameter,
                    MJDParameter,
                    AngleParameter,
                ),
            )
            or hasattr(param, "prefix")
            or param.quantity is None
        ):
            continue

        scale_factor = param.tcb2tdb_scale_factor
        value = (
            param.value * day_to_s
            if isinstance(param, MJDParameter)
            else (param.quantity * scale_factor).si.value
        )
        dim = param.effective_dimensionality

        original_units = str(param.units)
        unit_conversion_factor = (param.units * scale_factor / u.s**dim).to_value(
            u.dimensionless_unscaled, equivalencies=u.dimensionless_angles()
        )

        single_params.append(
            {
                "name": param_name,
                "default_value": float(value),
                "dimension": dim,
                "frozen": param.frozen,
                "original_units": original_units,
                "unit_conversion_factor": float(unit_conversion_factor),
            }
        )

    # Process pseudo_single_params
    multi_params = []
    processed_multi_params = []
    for param_name in pseudo_single_params:
        if param_name in model and model[param_name].quantity is not None:
            prefix_param_names = [param_name] + list(
                model.get_prefix_mapping(param_name).values()
            )

            elements = _get_multiparam_elements(model, prefix_param_names)

            multi_params.append({"name": param_name, "elements": elements})
            processed_multi_params.append(param_name)

    # Process ordinary multi parameters
    for param_name in model.params:
        param = model[param_name]

        # Skip things that are not ordinary multi parameters
        if (
            param_name in ignore_params
            or not hasattr(param, "prefix")
            or param.prefix in pseudo_single_params
            or param.prefix in processed_multi_params
            or param.quantity is None
        ):
            continue

        prefix_param_names = list(model.get_prefix_mapping(param.prefix).values())

        elements = _get_multiparam_elements(model, prefix_param_names)

        multi_params.append({"name": param.prefix, "elements": elements})
        processed_multi_params.append(param.prefix)

    return single_params, multi_params


def components_from_model(model: TimingModel) -> list:
    component_names = list(model.components.keys())

    components = []
    for component_name in component_names:
        if component_name in ["AbsPhase", "SolarSystemShapiro"]:
            continue
        elif component_name == "Spindown":
            components.append({"name": "Spindown"})
        elif component_name.startswith("Astrometry"):
            components.append(
                {
                    "name": "SolarSystem",
                    "ecliptic_coordinates": (component_name == "AstrometryEcliptic"),
                    "planet_shapiro": model.PLANET_SHAPIRO.value,
                }
            )
        elif component_name == "SolarWindDispersion":
            components.append(
                {"name": "SolarWindDispersion", "model": int(model.SWM.value)}
            )
        elif component_name == "DispersionDM":
            components.append({"name": "DispersionTaylor"})
        elif component_name == "TroposphereDelay" and model.CORRECT_TROPOSPHERE.value:
            components.append({"name": "Troposphere"})
        else:
            components.append({"name": component_name})
    return components


def infodict_from_model(model: TimingModel) -> dict:
    return {
        "pulsar_name": model.PSR.value if model.PSR.value is not None else "",
        "ephem": model.EPHEM.value,
        "clock": model.CLOCK.value,
        "units": model.UNITS.value,
    }


if __name__ == "__main__":
    par, tim, prefix = sys.argv[1:]

    setup_log(level="WARNING")

    model: TimingModel
    toas: TOAs

    model, toas = get_model_and_toas(par, tim, planets=True, add_tzr_to_model=True)
    toas.compute_pulse_numbers(model)
    fix_params(model)

    filename = f"{prefix}.hdf5"

    toas_table = toas_to_table(toas)

    components_str = json.dumps(components_from_model(model))
    params_str = json.dumps(params_from_model(model))
    info_str = json.dumps(infodict_from_model(model))

    tzr_toa: TOAs = model.get_TZR_toa(toas)
    tzr_toa.compute_pulse_numbers(model)
    tzr_table = toas_to_table(tzr_toa)

    toas_table.write(filename, path="TOAs", format="hdf5", overwrite=True)

    with h5.File(filename, "a") as f:
        f.create_dataset("Components", data=components_str)
        f.create_dataset("Parameters", data=params_str)
        f.create_dataset("Info", data=info_str)

    tzr_table.write(filename, path="TZRTOA", format="hdf5", append=True)
