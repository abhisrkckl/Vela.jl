from typing import List, Optional, Tuple
from pint.logging import setup as setup_log
from pint.models import PhaseOffset, TimingModel, get_model_and_toas
from pint.models.parameter import MJDParameter
from pint.toa import TOAs
from pint.binaryconvert import convert_binary

from .vela import vl
from .convert_components import pint_components_to_vela
from .convert_toas import pint_toa_to_vela, pint_toas_to_vela
from .convert_parameters import pint_parameters_to_vela
from .priors import get_default_priors
from .ecorr import ecorr_sort


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

    model["PHOFF"].frozen = False
    model["PHOFF"].uncertainty_value = 0.1

    if (
        "H4" in model
        and model["H4"].quantity is not None
        and model["STIGMA"].quantity is None
    ):
        model["STIGMA"].quantity = model["H4"].quantity / model["H3"].quantity
        model["STIGMA"].frozen = model["H4"].frozen
        model["H4"].frozen = True

    zeroable_params = [
        "M2",
        "SINI",
        "PBDOT",
        "XPBDOT",
        "A1DOT",
        "EPS1DOT",
        "EPS2DOT",
        "H3",
        "STIGMA",
    ]
    for p in zeroable_params:
        if p in model and model[p].quantity is None:
            model[p].value = 0


def get_kernel(
    model: TimingModel,
    ecorr_toa_ranges: List[Tuple[int, int]],
    ecorr_indices: List[int],
):
    if not model.has_correlated_errors:
        return vl.WhiteNoiseKernel()
    elif not model.has_time_correlated_errors:
        ecorr_groups = vl.Vector(
            [
                vl.EcorrGroup(start, stop, index)
                for (start, stop), index in zip(ecorr_toa_ranges, ecorr_indices)
            ]
        )
        return vl.EcorrKernel(ecorr_groups)
    raise NotImplementedError("Correlated noise kernels are not yet implemented.")


def pint_model_to_vela(
    model: TimingModel,
    toas: TOAs,
    cheat_prior_scale: float,
    custom_prior_dists: dict,
    ecorr_toa_ranges: Optional[List[Tuple[int, int]]] = None,
    ecorr_indices: Optional[List[Tuple[int]]] = None,
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

    kernel = get_kernel(model, ecorr_toa_ranges, ecorr_indices)

    return vl.TimingModel(
        pulsar_name,
        model.EPHEM.value,
        model.CLOCK.value,
        model.UNITS.value,
        components,
        kernel,
        param_handler,
        tzr_toa,
        priors,
    )


def read_model_and_toas(
    parfile: str, timfile: str, cheat_prior_scale=20, custom_prior_dicts={}
):
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

    if "BinaryBT" in mp.components:
        mp = convert_binary(mp, "DD")

    if "EcorrNoise" in mp.components:
        tp, ecorr_toa_ranges, ecorr_indices = ecorr_sort(mp, tp)
    else:
        ecorr_toa_ranges, ecorr_indices = None, None

    model = pint_model_to_vela(
        mp,
        tp,
        cheat_prior_scale,
        custom_prior_dicts,
        ecorr_toa_ranges=ecorr_toa_ranges,
        ecorr_indices=ecorr_indices,
    )
    toas = pint_toas_to_vela(tp)

    return model, toas


def par_tim_to_jlso(parfile: str, timfile: str, jlsofile: str):
    mv, tv = read_model_and_toas(parfile, timfile)
    vl.save_pulsar_data(jlsofile, mv, tv)
