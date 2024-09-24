import re
from typing import List

from astropy import units as u
from astropy.time import Time
from pint import DMconst
from pint.models import TimingModel
from pint.models.parameter import (
    AngleParameter,
    MJDParameter,
    Parameter,
    compute_effective_dimensionality,
    floatParameter,
    funcParameter,
)

from .convert_toas import day_to_s
from .vela import jl, to_jldd, vl


fdjump_rx = re.compile("^FD(\\d+)JUMP(\\d+)")


def get_scale_factor(param: Parameter):
    """Get the scale factor that converts a given parameter value from the `PINT`
    units to the `Vela` units.

    In `Vela`'s internal representation, all quantities have dimensions of [T^n]
    with units of second^n for some n∈ℕ. This is because all measurable quantities
    in pulsar timing are times, frequencies, or phases, and time-derivatives thereof.
    All of these can be written with units second^n for some n∈ℕ.

    Most of these scale factors are the same as the TCB <-> TDB conversion factors discussed
    in https://nanograv-pint.readthedocs.io/en/latest/tcb2tdb-factors.html.
    """

    if param.tcb2tdb_scale_factor is not None:
        return param.tcb2tdb_scale_factor
    elif isinstance(param.quantity, Time):
        return 1
    elif param.name in ["CM"] or (
        hasattr(param, "prefix")
        and param.prefix in ["CM", "CMWXSIN_", "CMWXCOS_", "DMJUMP", "DMEQUAD"]
    ):
        return DMconst
    elif (
        param.name in ["TNCHROMIDX"]
        or (
            hasattr(param, "prefix")
            and param.prefix
            in [
                "EFAC",
                "EQUAD",
                "ECORR",
                "DMEFAC",
                "FD",
            ]
        )
        or fdjump_rx.match(param.name)
    ):
        return 1
    else:
        raise ValueError(f"Unable to estimate scale factor for {param.name}.")


def pint_parameter_to_vela(param: Parameter, epoch_mjd: float):
    """Construct a `Vela.Parameter` object from a `PINT` `Parameter` object."""

    scale_factor = get_scale_factor(param)
    default_value = (
        (param.value - epoch_mjd) * day_to_s
        if isinstance(param.quantity, Time)
        else (param.quantity * scale_factor).si.value
    )

    dim = (
        compute_effective_dimensionality(param.quantity, scale_factor)
        if not isinstance(param.quantity, Time)
        else 1
    )

    original_units = str(param.units)
    unit_conversion_factor = (param.units * scale_factor / u.s**dim).to_value(
        u.dimensionless_unscaled, equivalencies=u.dimensionless_angles()
    )

    if param.name != "F0":
        return vl.Parameter(
            jl.Symbol(param.name),
            vl.GQ[dim](default_value),
            param.frozen,
            original_units,
            unit_conversion_factor,
        )
    else:
        return vl.Parameter(
            jl.Symbol(param.name),
            vl.GQ[dim](to_jldd(default_value).lo),
            param.frozen,
            original_units,
            unit_conversion_factor,
        )


def _get_multiparam_elements(model: TimingModel, multi_param_names: List[str]):
    """Construct a Julia `Vector` of `Vela.Parameter` objects from a list of
    `PINT` parameter names. This is usually done for `prefixParameters` and
    `maskParameters`."""

    elements = []
    for pxname in multi_param_names:
        pxparam = model[pxname]

        if pxparam.value is None:
            break

        elements.append(pint_parameter_to_vela(pxparam, float(model.PEPOCH.value)))

    return jl.Vector[vl.Parameter](elements)


# Parameters that are treated as single `Parameter`s in `PINT` due to compatibility reasons,
# but logically should have been `prefixParameter`s.
pseudo_single_params = ["DM", "CM"]


def pint_parameters_to_vela(model: TimingModel):
    """Convert the parameters in a `PINT` `TimingModel` to their `Vela` representation.

    Single parameters are represented as `Vela.Parameter` objects.
    `prefixParameter`s and `maskParameter`s are represented as `Vela.MultiParameter` objects.

    Returns Julia `Vector`s of `Vela.Parameter`s and `Vela.MultiParameter`s.
    """

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
        "H4",  # Handled separately; converted to STIG.
        "DMX",  # The actual DMX parameters are "DMX_".
        "PEPOCH",  # Included separately. We subtract PEPOCH from all TOAs and MJDParameters.
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

        single_params.append(pint_parameter_to_vela(param, float(model.PEPOCH.value)))

    single_params.append(
        vl.Parameter(
            jl.Symbol("F_"),
            vl.GQ[-1](to_jldd(model.F0.quantity.si.value).hi),
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

    # Process ordinary multi parameters (does not include `FDJUMP`s)
    for param_name in model.params:
        param = model[param_name]

        # Skip things that are not ordinary multi parameters
        if (
            param_name in ignore_params
            or not hasattr(param, "prefix")
            or param.prefix in pseudo_single_params
            or param.prefix in processed_multi_params
            or param.quantity is None
            or fdjump_rx.match(param.name)
        ):
            continue

        prefix_param_names = list(model.get_prefix_mapping(param.prefix).values())

        elements = _get_multiparam_elements(model, prefix_param_names)

        multi_params.append(vl.MultiParameter(jl.Symbol(param.prefix), elements))
        processed_multi_params.append(param.prefix)

    # Special treatment for `FDJUMP`s
    if "FDJump" in model.components:
        fdjump_names = [
            fdj
            for fdj in model.components["FDJump"].fdjumps
            if model[fdj].quantity is not None
        ]
        elements = _get_multiparam_elements(model, fdjump_names)
        multi_params.append(vl.MultiParameter(jl.Symbol("FDJUMP"), elements))

    multi_params = jl.Vector[vl.MultiParameter](multi_params)

    return single_params, multi_params
