from astropy import units as u
from typing import List

from pint.models import TimingModel
from pint.models.parameter import (
    AngleParameter,
    MJDParameter,
    floatParameter,
    funcParameter,
    Parameter,
    compute_effective_dimensionality,
)
from pint import DMconst

from .convert_toas import day_to_s
from .vela import vl, jl, to_jldd


def get_scale_factor(param):
    if param.tcb2tdb_scale_factor is not None:
        return param.tcb2tdb_scale_factor
    elif isinstance(param, MJDParameter):
        return 1
    elif param.name == "CM" or (
        hasattr(param, "prefix") and param.prefix in ["CM", "CMWXSIN_", "CMWXCOS_"]
    ):
        return DMconst
    elif hasattr(param, "prefix") and param.prefix in ["EFAC", "EQUAD", "ECORR"]:
        return 1
    elif param.name == "TNCHROMIDX":
        return 1
    else:
        raise ValueError(f"Unable to estimate scale factor for {param.name}.")


def pint_parameter_to_vela(param: Parameter):
    scale_factor = get_scale_factor(param)
    default_value = (
        param.value * day_to_s
        if isinstance(param, MJDParameter)
        else (param.quantity * scale_factor).si.value
    )

    dim = (
        compute_effective_dimensionality(param.quantity, scale_factor)
        if not isinstance(param, MJDParameter)
        else 1
    )

    original_units = str(param.units)
    unit_conversion_factor = (param.units * scale_factor / u.s**dim).to_value(
        u.dimensionless_unscaled, equivalencies=u.dimensionless_angles()
    )

    if param.name != "F0":
        return vl.Parameter(
            jl.Symbol(param.name),
            vl.GQ(default_value, dim),
            param.frozen,
            original_units,
            unit_conversion_factor,
        )
    else:
        return vl.Parameter(
            jl.Symbol(param.name),
            vl.GQ(to_jldd(default_value).lo, dim),
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
            or isinstance(param, funcParameter)
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
            vl.GQ(to_jldd(model.F0.quantity.si.value).hi, -1),
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
