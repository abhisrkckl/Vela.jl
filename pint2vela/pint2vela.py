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
    compute_effective_dimensionality,
)
from pint.toa import TOAs

from .vela import jl, vl
from .convert_toas import pint_toa_to_vela, pint_toas_to_vela, day_to_s

def pint_parameter_to_vela(param: Parameter):
    scale_factor = (
        param.tcb2tdb_scale_factor if param.tcb2tdb_scale_factor is not None else 1
    )
    default_value = (
        param.value * day_to_s
        if isinstance(param, MJDParameter)
        else (param.quantity * scale_factor).si.value
    )

    dim = (
        param.effective_dimensionality
        if param.tcb2tdb_scale_factor is not None
        else compute_effective_dimensionality(param.quantity, 1)
    )

    original_units = str(param.units)
    unit_conversion_factor = (param.units * scale_factor / u.s**dim).to_value(
        u.dimensionless_unscaled, equivalencies=u.dimensionless_angles()
    )

    return (
        vl.Parameter(
            jl.Symbol(param.name),
            vl.GQ(default_value, dim),
            param.frozen,
            original_units,
            unit_conversion_factor,
        )
        if param.name != "F0"
        else vl.Parameter(
            jl.Symbol(param.name),
            vl.GQ(0.0, dim),
            param.frozen,
            original_units,
            unit_conversion_factor,
        )
    )


def _get_multiparam_elements(model: TimingModel, multi_param_names: List[str]):
    elements = []
    for pxname in multi_param_names:
        pxparam = model[pxname]

        if pxparam.value is None:
            break

        elements.append(pint_parameter_to_vela(pxparam))

    return jl.Vector[vl.Parameter](elements)


pseudo_single_params = ["DM", "CM"]


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
    single_params.append(
        vl.Parameter(
            jl.Symbol("F_"),
            vl.GQ(model.F0.quantity.si.value, -1),
            True,
            str(model.F0.units),
            1.0,
        )
    )
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


def pint_components_to_vela(model: TimingModel, toas: TOAs):
    component_names = list(model.components.keys())

    components = []

    # The order below is important.
    # It goes from the observatory to the pulsar.
    # The general order is DelayComponents -- PhaseComponents -- NoiseComponents.

    if "TroposphereDelay" in component_names and model.CORRECT_TROPOSPHERE.value:
        components.append(vl.Troposphere())

    if "AstrometryEcliptic" in component_names:
        components.append(vl.SolarSystem(True, model.PLANET_SHAPIRO.value))
    elif "AstrometryEquatorial" in component_names:
        components.append(vl.SolarSystem(False, model.PLANET_SHAPIRO.value))

    if "SolarWindDispersion" in component_names and model.NE_SW.value != 0:
        components.append(vl.SolarWindDispersion(int(model.SWM.value)))

    if "DispersionDM" in component_names:
        components.append(vl.DispersionTaylor())

    if "Spindown" in component_names:
        components.append(vl.Spindown())

    if "PhaseOffset" in component_names:
        components.append(vl.PhaseOffset())

    if "PhaseJump" in component_names:
        masks = []
        for jump in model.get_jump_param_objects():
            mask = np.repeat(False, len(toas))
            mask[jump.select_toa_mask(toas)] = True
            assert any(mask), f"{jump.name} has no TOAs!"
            masks.append(mask)
        masks = jl.BitMatrix(np.array(masks))
        components.append(vl.PhaseJump(masks))

    if "ScaleToaError" in component_names:
        efac_masks = [model[efac].select_toa_mask(toas) for efac in model.EFACs]
        equad_masks = [model[equad].select_toa_mask(toas) for equad in model.EQUADs]

        efac_index_mask = []
        equad_index_mask = []
        for ii in range(len(toas)):
            efac_idxs = np.argwhere([(ii in mask) for mask in efac_masks]).flatten()
            equad_idxs = np.argwhere([(ii in mask) for mask in equad_masks]).flatten()

            assert len(efac_idxs) in [0, 1], "EFAC groups must not overlap."
            assert len(equad_idxs) in [0, 1], "EQUAD groups must not overlap."

            efac_index_mask.append(efac_idxs.item() + 1 if len(efac_idxs) == 1 else 0)
            equad_index_mask.append(
                equad_idxs.item() + 1 if len(equad_idxs) == 1 else 0
            )

        efac_index_mask = jl.Vector[jl.UInt](efac_index_mask)
        equad_index_mask = jl.Vector[jl.UInt](equad_index_mask)

        assert len(set(efac_index_mask)) in [
            len(model.EFACs),
            len(model.EFACs) + 1,
        ], "EFAC without any TOAs found!"
        assert len(set(equad_index_mask)) in [
            len(model.EQUADs),
            len(model.EQUADs) + 1,
        ], "EQUAD without any TOAs found!"

        components.append(vl.MeasurementNoise(efac_index_mask, equad_index_mask))

    components = jl.Tuple(components)

    return components


def get_default_prior(
    model: TimingModel,
    param_name: str,
    cheat_prior_scale: float,
    custom_prior_dists: dict,
):
    assert param_name in model.free_params

    param = model[param_name]

    scale_factor = (
        param.tcb2tdb_scale_factor if param.tcb2tdb_scale_factor is not None else 1
    )

    assert param.uncertainty is not None and param.uncertainty > 0

    val = (
        param.value * day_to_s
        if isinstance(param, MJDParameter)
        else (param.quantity * scale_factor).si.value
    )
    err = (param.uncertainty * scale_factor).si.value

    pmin = val - cheat_prior_scale * err
    pmax = val + cheat_prior_scale * err

    pdist = custom_prior_dists.get(param_name, jl.Uniform(pmin, pmax))

    if (
        isinstance(param, (floatParameter, MJDParameter, AngleParameter))
        and param_name not in pseudo_single_params
    ):
        pname = jl.Symbol(param_name)
        return vl.SimplePrior[pname](pdist)
    elif hasattr(param, "prefix") and param.prefix not in pseudo_single_params:
        pname = jl.Symbol(param.prefix)
        index = (
            list(model.get_prefix_mapping(param.prefix).values()).index(param_name) + 1
        )
        return vl.SimplePriorMulti[pname, index](pdist)
    elif param_name in pseudo_single_params:
        pname = jl.Symbol(param_name)
        return vl.SimplePriorMulti[pname, 1](pdist)
    elif param.prefix in pseudo_single_params:
        pname = jl.Symbol(param.prefix)
        index = (
            list(model.get_prefix_mapping(param.prefix).values()).index(param_name) + 2
        )
        return vl.SimplePriorMulti[pname, index](pdist)


def get_default_priors(
    model: TimingModel,
    free_params: List[str],
    cheat_prior_scale: float,
    custom_prior_dists: dict,
):
    return tuple(
        get_default_prior(model, param_name, cheat_prior_scale, custom_prior_dists)
        for param_name in free_params
    )


def pint_model_to_vela(
    model: TimingModel,
    toas: TOAs,
    cheat_prior_scale: float = 10.0,
    custom_prior_dists: dict = {},
):
    toas.compute_pulse_numbers(model)
    fix_params(model)

    pulsar_name = model.PSR.value if model.PSR.value is not None else ""

    components = pint_components_to_vela(model, toas)

    single_params, multi_params = pint_parameters_to_vela(model)
    param_handler = vl.ParamHandler(single_params, multi_params)

    free_params = vl.get_free_param_names(param_handler)

    priors = get_default_priors(
        model, free_params, cheat_prior_scale, custom_prior_dists
    )

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
        priors,
    )


def read_model_and_toas(parfile: str, timfile: str):
    """Read a pair of par & tim files and create a `Vela.TimingModel` object and a
    Julia `Vector` of `TOA`s."""
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
