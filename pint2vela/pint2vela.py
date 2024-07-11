from typing import List

from astropy import units as u
import numpy as np
from pint.logging import setup as setup_log
from pint.models import PhaseOffset, TimingModel, get_model_and_toas
from pint.models.parameter import (
    AngleParameter,
    MJDParameter,
    floatParameter,
    Parameter,
)
from pint.toa import TOAs

from juliacall import Main as jl

jl.seval("import Vela")
vl = jl.Vela

day_to_s = 86400


def pint_toa_to_vela(toas: TOAs, idx: int, tzr: bool = False):
    assert toas.planets
    assert toas.get_pulse_numbers() is not None

    tdb_ld = toas.table["tdbld"].value[idx]
    tdb_str = np.format_float_positional(tdb_ld, unique=True)
    tdb = vl.time(jl.parse(vl.Double64, tdb_str) * day_to_s)

    phase = (
        toas.table["pulse_number"].value[idx]
        - toas.table["delta_pulse_number"].value[idx]
    )
    phase = vl.dimensionless(jl.parse(vl.Double64, str(phase)))

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

    return vl.TOA(tdb, err, freq, phase, False, tzr, ephem)


def pint_toas_to_vela(toas: TOAs):
    return jl.Vector[vl.TOA]([pint_toa_to_vela(toas, idx) for idx in range(len(toas))])


def pint_parameter_to_vela(param: Parameter):
    scale_factor = param.tcb2tdb_scale_factor
    default_value = (
        param.value * day_to_s
        if isinstance(param, MJDParameter)
        else (param.quantity * scale_factor).si.value
    )
    dim = param.effective_dimensionality
    original_units = str(param.units)
    unit_conversion_factor = (param.units * scale_factor / u.s**dim).to_value(
        u.dimensionless_unscaled, equivalencies=u.dimensionless_angles()
    )

    return vl.Parameter(
        jl.Symbol(param.name),
        vl.GQ(default_value, dim),
        param.frozen,
        original_units,
        unit_conversion_factor,
    )


def _get_multiparam_elements(model: TimingModel, multi_param_names: List[str]):
    elements = []
    for pxname in multi_param_names:
        pxparam = model[pxname]

        if pxparam.value is None:
            break

        elements.append(pint_parameter_to_vela(pxparam))

    return jl.Vector[vl.Parameter](elements)


def pint_parameters_to_vela(model: TimingModel):
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

        single_params.append(pint_parameter_to_vela(param))
    single_params = jl.Vector[vl.Parameter](single_params)

    # Process pseudo_single_params
    multi_params = []
    processed_multi_params = []
    for param_name in pseudo_single_params:
        if param_name in model and model[param_name].quantity is not None:
            prefix_param_names = [param_name] + list(
                model.get_prefix_mapping(param_name).values()
            )

            elements = _get_multiparam_elements(model, prefix_param_names)

            multi_params.append(vl.MultiParameter(jl.Symbol(param_name), elements))
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

        multi_params.append(vl.MultiParameter(jl.Symbol(param.prefix), elements))
        processed_multi_params.append(param.prefix)
    multi_params = jl.Vector[vl.MultiParameter](multi_params)

    return single_params, multi_params


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


def pint_components_to_vela(model: TimingModel):
    component_names = list(model.components.keys())

    components = []
    for component_name in component_names:
        if component_name in ["AbsPhase", "SolarSystemShapiro"]:
            continue
        elif component_name == "Spindown":
            components.append(vl.Spindown())
        elif component_name == "PhaseOffset":
            components.append(vl.PhaseOffset())
        elif component_name.startswith("Astrometry"):
            ecliptic_coordinates = component_name == "AstrometryEcliptic"
            components.append(
                vl.SolarSystem(ecliptic_coordinates, model.PLANET_SHAPIRO.value)
            )
        elif component_name == "SolarWindDispersion":
            components.append(vl.SolarWindDispersion(int(model.SWM.value)))
        elif component_name == "DispersionDM":
            components.append(vl.DispersionTaylor())
        elif component_name == "TroposphereDelay" and model.CORRECT_TROPOSPHERE.value:
            components.append(vl.Troposphere())
        # else:
        #     components.append({"name": component_name})
    components = jl.Tuple(components)

    return components


def pint_model_to_vela(model: TimingModel, toas: TOAs):
    toas.compute_pulse_numbers(model)
    fix_params(model)

    pulsar_name = model.PSR.value if model.PSR.value is not None else ""

    components = pint_components_to_vela(model)

    single_params, multi_params = pint_parameters_to_vela(model)
    param_handler = vl.ParamHandler(single_params, multi_params)

    tzr_toa = model.get_TZR_toa(toas)
    tzr_toa.compute_pulse_numbers(model)
    tzr_toa = pint_toa_to_vela(tzr_toa, 0, tzr=True)

    return vl.TimingModel(
        pulsar_name,
        model.EPHEM.value,
        model.CLOCK.value,
        model.UNITS.value,
        components,
        param_handler,
        tzr_toa,
    )


def read_model_and_toas(parfile: str, timfile: str):
    setup_log(level="WARNING")
    mp, tp = get_model_and_toas(
        parfile,
        timfile,
        planets=True,
        allow_tcb=True,
        allow_T2=True,
        add_tzr_to_model=True,
    )

    model = pint_model_to_vela(mp, tp)
    toas = pint_toas_to_vela(tp)

    return model, toas
