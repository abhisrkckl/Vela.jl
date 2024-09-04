from typing import List

from pint.models import TimingModel
from pint.models.parameter import floatParameter, MJDParameter, AngleParameter

from .vela import jl, vl
from .convert_toas import day_to_s
from .convert_parameters import pseudo_single_params

DEFAULT_PRIOR_DISTS = {
    "PHOFF": jl.Uniform(-0.5, 0.5),
    "EFAC": jl.Uniform(0.1, 5.0),
    "EQUAD": jl.Uniform(0.0, 1e-4),
    "SINI": jl.Uniform(0.0, 1.0),
    "STIGMA": jl.Uniform(0.0, 1.0),
}


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

    if param_name in custom_prior_dists:
        pdist = custom_prior_dists[param_name]
    elif hasattr(param, "prefix") and param.prefix in custom_prior_dists:
        pdist = custom_prior_dists[param.prefix]
    elif param_name in DEFAULT_PRIOR_DISTS:
        pdist = DEFAULT_PRIOR_DISTS[param_name]
    elif hasattr(param, "prefix") and param.prefix in DEFAULT_PRIOR_DISTS:
        pdist = DEFAULT_PRIOR_DISTS[param.prefix]
    else:
        val = (
            (
                param.value * day_to_s
                if isinstance(param, MJDParameter)
                else (param.quantity * scale_factor).si.value
            )
            if param_name != "F0"
            else 0.0
        )

        assert (
            param.uncertainty is not None and param.uncertainty > 0
        ), f"Uncertainty not given for {param_name}."
        err = (param.uncertainty * scale_factor).si.value

        pmin = val - cheat_prior_scale * err
        pmax = val + cheat_prior_scale * err

        pdist = jl.Uniform(pmin, pmax)

    if (
        isinstance(param, (floatParameter, MJDParameter, AngleParameter))
        and param_name not in pseudo_single_params
        and not hasattr(param, "prefix")
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
