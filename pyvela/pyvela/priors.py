import json
from typing import IO, List

import numpy as np
from astropy.time import Time
from pint.models import TimingModel
from pint.models.parameter import AngleParameter, MJDParameter, floatParameter

from .parameters import (
    fdjump_rx,
    get_scale_factor,
    get_unit_conversion_factor,
    pseudo_single_params,
)
from .toas import day_to_s
from .vela import jl, vl

DEFAULT_PRIOR_DISTS = {
    "PHOFF": jl.Uniform(-0.5, 0.5),
    "EFAC": jl.LogNormal(0.0, 0.25),
    "KOM": jl.Uniform(0.0, 2 * jl.pi),
    "KIN": vl.KINPriorDistribution(),
    "SINI": vl.SINIPriorDistribution(),
    "STIGMA": vl.STIGMAPriorDistribution(),
    "SHAPMAX": vl.SHAPMAXPriorDistribution(),
    "DMEFAC": jl.LogNormal(0.0, 0.25),
    "PLREDCOS_": jl.Normal(),
    "PLREDSIN_": jl.Normal(),
    "TNREDAMP": jl.Uniform(-16.0, -9.0),
    "TNREDGAM": jl.Uniform(0.0, 7.0),
    "PLDMCOS_": jl.Normal(),
    "PLDMSIN_": jl.Normal(),
    "TNDMAMP": jl.Uniform(-16.0, -9.0),
    "TNDMGAM": jl.Uniform(0.0, 7.0),
    "PLCHROMCOS_": jl.Normal(),
    "PLCHROMSIN_": jl.Normal(),
    "TNCHROMAMP": jl.Uniform(-16.0, -9.0),
    "TNCHROMGAM": jl.Uniform(0.0, 7.0),
}


def get_default_prior(
    model: TimingModel,
    param_name: str,
    epoch_mjd: float,
    cheat_prior_scale: float,
    custom_prior_dists: dict,
):
    """Returns a `Vela.Prior` object corresponding to a free model parameter.
    If the parameter is included in `custom_prior_dists`, the custom prior is
    used. Otherwise, if the parameter is included in `DEFAULT_PRIOR_DISTS`, the default
    prior is used. Otherwise, a 'cheat' prior is used, which corresponds to a uniform
    distribution centered at the default value in the `model`, and a width that is
    `2 * cheat_prior_scale` times the uncertainty quoted in the model."""

    assert param_name in model.free_params

    param = model[param_name]

    scale_factor = get_scale_factor(param)

    if param_name in custom_prior_dists:
        pdist = custom_prior_dists[param_name]
    # elif hasattr(param, "prefix") and param.prefix in custom_prior_dists:
    #     pdist = custom_prior_dists[param.prefix]
    elif param_name in DEFAULT_PRIOR_DISTS:
        pdist = DEFAULT_PRIOR_DISTS[param_name]
    elif hasattr(param, "prefix") and param.prefix in DEFAULT_PRIOR_DISTS:
        pdist = DEFAULT_PRIOR_DISTS[param.prefix]
    else:
        val = (
            (
                (param.value - epoch_mjd) * day_to_s
                if isinstance(param.quantity, Time)
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
    elif hasattr(param, "prefix") and fdjump_rx.match(param.name):
        pname = jl.Symbol("FDJUMP")
        fdjump_names = [
            fdj
            for fdj in model.components["FDJump"].fdjumps
            if model[fdj].quantity is not None
        ]
        index = fdjump_names.index(param.name) + 1
        return vl.SimplePriorMulti[pname, index](pdist)
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
    epoch_mjd: float,
    cheat_prior_scale: float,
    custom_prior_dists: dict,
) -> tuple:
    """Returns a tuple of `Vela.Prior` objects corresponding to each free parameter."""

    params_disallow_custom_priors = {"PLREDCOS_", "PLREDSIN_", "PLDMCOS_", "PLDMSIN_"}
    assert (
        len(params_disallow_custom_priors.intersection(custom_prior_dists.keys())) == 0
    ), f"Custom priors cannot be set for the following parameters: {params_disallow_custom_priors.intersection(custom_prior_dists.keys())}"

    return tuple(
        get_default_prior(
            model, param_name, epoch_mjd, cheat_prior_scale, custom_prior_dists
        )
        for param_name in free_params
    )


def process_custom_priors(custom_priors_raw: dict, model: TimingModel) -> dict:
    output_dict = dict()
    for par in model.free_params:
        prior_info: dict
        if par in custom_priors_raw:
            prior_info = custom_priors_raw[par]
        elif hasattr(model[par], "prefix") and model[par].prefix in custom_priors_raw:
            prior_info = custom_priors_raw[model[par].prefix]
        else:
            continue

        unit_conversion_factor = get_unit_conversion_factor(model[par])

        distr_type = getattr(jl.Distributions, prior_info["distribution"])
        args = np.array(prior_info["args"])

        scaled_args = vl.scale_prior_args(distr_type, args, unit_conversion_factor)

        if "upper" in prior_info and "lower" in prior_info:
            output_dict[par] = jl.Distributions.truncated(
                distr_type(*scaled_args),
                lower=prior_info["lower"],
                upper=prior_info["upper"],
            )
        elif "upper" in prior_info:
            output_dict[par] = jl.Distributions.truncated(
                distr_type(*scaled_args), upper=prior_info["upper"]
            )
        elif "lower" in prior_info:
            output_dict[par] = jl.Distributions.truncated(
                distr_type(*scaled_args), lower=prior_info["lower"]
            )
        else:
            output_dict[par] = distr_type(*scaled_args)

    return output_dict
