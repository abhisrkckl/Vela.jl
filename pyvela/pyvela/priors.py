from typing import List

import astropy.constants as const
import astropy.units as u
import numpy as np
from astropy.time import Time
from pint import DMconst, dmu
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

dmunit = (DMconst * dmu).to_value(u.Hz)  # DM unit in geometric units
Rmax = (52850 * u.lightyear / const.c).to_value(
    u.s
)  # Radius of the Galaxy in geometric units
Msun = (const.G * u.Msun / const.c**3).to_value(u.s)  # Solar mass in geometric units

# Some of these prior distributions are based on physical considerations.
# Others are based on typical values found in millisecond pulsars. They
# may not be valid for slow pulsars.
DEFAULT_PRIOR_DISTS = {
    "PHOFF": jl.Uniform(
        -0.5, 0.5
    ),  # Assumes phase connection and an appropriate TZRMJD value.
    "EFAC": jl.LogNormal(0.0, 0.25),  # Should be close to 1 most of the time.
    "EQUAD": jl.LogUniform(1e-9, 1e-4),  # Ballpark range based on PTA pulsars
    "ECORR": jl.LogUniform(1e-9, 1e-4),  # Ballpark range based on PTA pulsars
    "KOM": jl.Uniform(0.0, 2 * jl.pi),  # Physical prior
    "KIN": vl.KINPriorDistribution(),  # cos(ι) is uniformly distributed in [-1,1]
    "SINI": vl.SINIPriorDistribution(),  # cos(ι) is uniformly distributed in [0,1]
    "STIGMA": vl.STIGMAPriorDistribution(),  # cos(ι) is uniformly distributed in [0,1]
    "SHAPMAX": vl.SHAPMAXPriorDistribution(),  # cos(ι) is uniformly distributed in [0,1]
    "DMEFAC": jl.LogNormal(0.0, 0.25),  # Should be close to 1 most of the time.
    "DMEQUAD": jl.LogUniform(
        1e-8 * dmunit, 1e-2 * dmunit
    ),  # Ballpark range based on PTA pulsars
    "PLREDCOS_": jl.Normal(),  # Prior is part of the parameter definition
    "PLREDSIN_": jl.Normal(),  # Prior is part of the parameter definition
    "TNREDAMP": jl.Uniform(-18.0, -9.0),  # Ballpark range based on PTA pulsars
    "TNREDGAM": jl.Uniform(0.0, 7.0),  # Ballpark range based on PTA pulsars
    "PLDMCOS_": jl.Normal(),  # Prior is part of the parameter definition
    "PLDMSIN_": jl.Normal(),  # Prior is part of the parameter definition
    "TNDMAMP": jl.Uniform(-18.0, -9.0),  # Ballpark range based on PTA pulsars
    "TNDMGAM": jl.Uniform(0.0, 7.0),  # Ballpark range based on PTA pulsars
    "PLCHROMCOS_": jl.Normal(),  # Prior is part of the parameter definition
    "PLCHROMSIN_": jl.Normal(),  # Prior is part of the parameter definition
    "TNCHROMAMP": jl.Uniform(-18.0, -9.0),  # Ballpark range based on PTA pulsars
    "TNCHROMGAM": jl.Uniform(0.0, 7.0),  # Ballpark range based on PTA pulsars
    "DMX_": jl.Normal(0, 1e-2 * dmunit),  # Ballpark range based on PTA pulsars
    "PX": vl.PXPriorDistribution(
        2 * Rmax
    ),  # Uniformly distributed within a sphere of radius 2*Rmax.
    "M2": jl.LogUniform(1e-9 * Msun, 100 * Msun),  # Between a moon and a massive star.
    "H3": jl.LogUniform(
        1e-9 * Msun, 100 * Msun
    ),  # Same range as M2 because H3 = M2 * STIGMA**3 and STIGMA ∈ [0,1].
    "RAJ": jl.Uniform(0, 2 * jl.pi),  # duh!
    "DECJ": jl.LatitudePriorDistribution(),  # sin(DECJ) is uniformly distributed in [-1,1].
    "ELONG": jl.Uniform(0, 2 * jl.pi),  # duh!
    "ELAT": jl.LatitudePriorDistribution(),  # sin(ELAT) is uniformly distributed in [-1,1].
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

    assert (
        param_name in model.free_params
    ), "Refusing to construct prior for a non-free parameter."

    param = model[param_name]

    scale_factor = get_scale_factor(param)

    if param_name in custom_prior_dists:
        pdist = custom_prior_dists[param_name]
        source = vl.USER_DEFINED_PRIOR
    elif hasattr(param, "prefix") and param.prefix in custom_prior_dists:
        pdist = custom_prior_dists[param.prefix]
        source = vl.USER_DEFINED_PRIOR
    elif param_name in DEFAULT_PRIOR_DISTS:
        pdist = DEFAULT_PRIOR_DISTS[param_name]
        source = vl.DEFAULT_PRIOR
    elif hasattr(param, "prefix") and param.prefix in DEFAULT_PRIOR_DISTS:
        pdist = DEFAULT_PRIOR_DISTS[param.prefix]
        source = vl.DEFAULT_PRIOR
    elif param_name in ["TASC", "T0"]:
        # Both are `MJDParameter`s
        # The binary model is periodic in TASC/T0 with an approximate period of PB.
        # This prior picks one period.
        val = (param.value - epoch_mjd) * day_to_s
        PB = float(model.pb()[0].to_value(u.s))
        pmin = val - PB / 2
        pmax = val + PB / 2
        pdist = jl.Uniform(pmin, pmax)
        source = vl.DEFAULT_PRIOR
    elif hasattr(param, "prefix") and param.prefix in ["JUMP"]:
        # JUMP phase will be in [-0.5, 0.5] if the TOAs are phase connected.
        P = float((1 / model["F0"].quantity).to_value(u.s))
        pmin = -P / 2
        pmax = P / 2
        pdist = jl.Uniform(pmin, pmax)
        source = vl.DEFAULT_PRIOR
    else:
        # "Cheat" priors.
        # Must be careful while using these.
        val = (
            (
                (param.value - epoch_mjd) * day_to_s
                if isinstance(param.quantity, Time)
                else (param.quantity * scale_factor).si.value
            )
            if param_name != "F0"
            else 0.0
        )

        assert param.uncertainty is not None and param.uncertainty > 0, (
            f"Unable to construct prior for {param_name}. This can be resolved by "
            f"(a) defining a prior in the prior file or "
            f"(b) providing the frequentist uncertainty in the par file so that a 'cheat' prior can be used."
        )
        err = (param.uncertainty * scale_factor).si.value

        pmin = val - cheat_prior_scale * err
        pmax = val + cheat_prior_scale * err

        pdist = jl.Uniform(pmin, pmax)
        source = vl.CHEAT_PRIOR

    if (
        isinstance(param, (floatParameter, MJDParameter, AngleParameter))
        and param_name not in pseudo_single_params
        and not hasattr(param, "prefix")
    ):
        pname = jl.Symbol(param_name)
        return vl.SimplePrior[pname](pdist, source)
    elif hasattr(param, "prefix") and fdjump_rx.match(param.name):
        pname = jl.Symbol("FDJUMP")
        fdjump_names = [
            fdj
            for fdj in model.components["FDJump"].fdjumps
            if model[fdj].quantity is not None
        ]
        index = fdjump_names.index(param.name) + 1
        return vl.SimplePriorMulti[pname, index](pdist, source)
    elif hasattr(param, "prefix") and param.prefix not in pseudo_single_params:
        pname = jl.Symbol(param.prefix)
        index = (
            list(model.get_prefix_mapping(param.prefix).values()).index(param_name) + 1
        )
        return vl.SimplePriorMulti[pname, index](pdist, source)
    elif param_name in pseudo_single_params:
        pname = jl.Symbol(param_name)
        return vl.SimplePriorMulti[pname, 1](pdist, source)
    elif param.prefix in pseudo_single_params:
        pname = jl.Symbol(param.prefix)
        index = (
            list(model.get_prefix_mapping(param.prefix).values()).index(param_name) + 2
        )
        return vl.SimplePriorMulti[pname, index](pdist, source)


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
    ), (
        f"Custom priors cannot be set for the following parameters: {params_disallow_custom_priors.intersection(custom_prior_dists.keys())}. "
        f"The priors for these parameters are defined by the model."
    )

    return tuple(
        get_default_prior(
            model, param_name, epoch_mjd, cheat_prior_scale, custom_prior_dists
        )
        for param_name in free_params
    )


def process_custom_priors(custom_priors_raw: dict, model: TimingModel) -> dict:
    """Generate a Dict[str, jl.Distribution] from the input prior dictionary, usually
    read from a JSON file. This function unit conversions and the interpretation of
    `upper`/`lower` attributes."""
    output_dict = {}
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
                lower=(prior_info["lower"] * unit_conversion_factor),
                upper=(prior_info["upper"] * unit_conversion_factor),
            )
        elif "upper" in prior_info:
            output_dict[par] = jl.Distributions.truncated(
                distr_type(*scaled_args),
                upper=(prior_info["upper"] * unit_conversion_factor),
            )
        elif "lower" in prior_info:
            output_dict[par] = jl.Distributions.truncated(
                distr_type(*scaled_args),
                lower=(prior_info["lower"] * unit_conversion_factor),
            )
        else:
            output_dict[par] = distr_type(*scaled_args)

    return output_dict
