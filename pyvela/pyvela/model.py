from typing import List, Optional, Tuple

import numpy as np
from pint.models import PhaseOffset, TimingModel
from pint.models.parameter import MJDParameter, maskParameter
from pint.toa import TOAs

from .dmx import get_dmx_mask
from .gp_noise import PLChromNoiseGP, PLDMNoiseGP, PLRedNoiseGP
from .parameters import pint_parameters_to_vela
from .priors import get_default_priors
from .toas import day_to_s, pint_toa_to_vela
from .vela import jl, vl


def read_mask(toas: TOAs, params: List[maskParameter]) -> np.ndarray:
    """Read a TOA mask from a `maskParameter` in a `Vela`-friendly
    representation."""

    masks = []
    for param in params:
        mask = np.repeat(False, len(toas))
        mask[param.select_toa_mask(toas)] = True
        assert any(
            mask
        ), f"Mask parameter {param.name} has no TOAs! Please modify the par file to avoid such parameters."
        masks.append(mask)
    return np.array(masks)


def is_exclusive_mask(mask: np.ndarray) -> bool:
    """Check if the mask is exclusive. An exclusive mask is where one TOA
    belongs to only one group.

    For example, `EFAC`s and `EQUAD`s are generally exclusive, whereas `JUMP`s
    are sometimes not.
    """
    return all(map(lambda x: x in [0, 1], mask.sum(axis=0)))


def get_exclusive_mask(mask: np.ndarray) -> np.ndarray:
    """Convert a mask to its exclusive representation. Throws an error
    if the input is not exclusive."""
    result = []
    for m in mask.T:
        wh = np.where(m)[0]
        if len(wh) == 0:
            result.append(0)
        elif len(wh) == 1:
            result.append(wh.item() + 1)
        else:
            raise ValueError("The mask is not exclusive!")
    assert (
        len(result) == mask.shape[1]
    ), "Shape of the constructed (exclusive) index mask is inconsistent with its bit mask representation. This is a bug."
    return np.array(result)


def pint_components_to_vela(model: TimingModel, toas: TOAs):
    """Construct a tuple containing the `Vela.Component` objects from a
    pair of `PINT` `TimingModel` and `TOAs` objects.

    Most `Vela.Component` types have a one-to-one correspondence with their
    `PINT` counterparts. Exceptions include `SolarSystemShapiro`, which is
    implemented in `Vela.SolarSystem` along with astrometric delays.

    Note that only the TOA-uncorrelated `PINT` `Components` have their
    `Vela.Component` counterparts. The TOA-correlated components such as
    `EcorrNoise` are represented as `Vela.Kernel`s.

    Unlike `PINT` `Component`s, `Vela.Component`s are sometimes tightly coupled
    to a set of TOAs for performance reasons. Examples include `Vela.MeasurementNoise`
    and `Vela.PhaseJump`, where the TOA selection masks are precomputed and stored
    within the `Vela.Component` object.
    """

    component_names = list(model.components.keys())

    components = []

    # The order below is important.
    # It goes from the observatory to the pulsar.
    # The general order is DelayComponents -- PhaseComponents -- NoiseComponents.

    # if "TroposphereDelay" in component_names and model.CORRECT_TROPOSPHERE.value:
    #     components.append(vl.Troposphere())

    if "AstrometryEcliptic" in component_names:
        components.append(vl.SolarSystem(True, model["PLANET_SHAPIRO"].value))
    elif "AstrometryEquatorial" in component_names:
        components.append(vl.SolarSystem(False, model["PLANET_SHAPIRO"].value))

    if "SolarWindDispersion" in component_names and not (
        model["NE_SW"].value == 0 and model["NE_SW"].frozen
    ):
        components.append(vl.SolarWindDispersion())

    if "DispersionDM" in component_names:
        components.append(vl.DispersionTaylor())

    if "DispersionDMX" in component_names:
        dmx_mask = get_dmx_mask(model, toas)
        components.append(vl.DispersionPiecewise(dmx_mask))
    elif "DMWaveX" in component_names:
        components.append(vl.DMWaveX())
    elif "PLDMNoiseGP" in component_names:
        components.append(vl.PowerlawDispersionNoiseGP(int(model["TNDMC"].value)))

    if "FDJumpDM" in component_names:
        fdjumpdms = list(
            map(lambda pname: model[pname], model.components["FDJumpDM"].fdjump_dms)
        )
        masks0 = read_mask(toas, fdjumpdms)
        dmoff = (
            vl.ExclusiveDispersionOffset(jl.Vector[jl.UInt](get_exclusive_mask(masks0)))
            if is_exclusive_mask(masks0)
            else vl.DispersionOffset(jl.BitMatrix(masks0))
        )
        components.append(dmoff)

    if "DispersionJump" in component_names:
        dmjumps = list(
            map(lambda pname: model[pname], model.components["DispersionJump"].dm_jumps)
        )
        masks0 = read_mask(toas, dmjumps)
        dmjump = (
            vl.ExclusiveDispersionJump(jl.Vector[jl.UInt](get_exclusive_mask(masks0)))
            if is_exclusive_mask(masks0)
            else vl.DispersionJump(jl.BitMatrix(masks0))
        )
        components.append(dmjump)

    if "ChromaticCM" in component_names:
        components.append(vl.ChromaticTaylor())

    if "ChromaticCMX" in component_names:
        cmx_mask = get_dmx_mask(model, toas, param_prefix="CMX_")
        components.append(vl.ChromaticPiecewise(cmx_mask))
    elif "CMWaveX" in component_names:
        components.append(vl.CMWaveX())
    elif "PLChromNoiseGP" in component_names:
        components.append(vl.PowerlawChromaticNoiseGP(int(model.TNCHROMC.value)))

    if model.BINARY.value is not None:
        assert (model["PB"].quantity is not None) != (
            model["FB0"].quantity is not None
        ), "Expecting one and only one of PB and FB0. Please check the par file."
        use_fbx = model["FB0"].quantity is not None
        if "BinaryELL1" in component_names:
            components.append(vl.BinaryELL1(use_fbx))
        elif "BinaryELL1H" in component_names:
            components.append(vl.BinaryELL1H(use_fbx))
        elif "BinaryELL1k" in component_names:
            components.append(vl.BinaryELL1k(use_fbx))
        elif "BinaryDD" in component_names:
            components.append(vl.BinaryDD(use_fbx))
        elif "BinaryDDH" in component_names:
            components.append(vl.BinaryDDH(use_fbx))
        elif "BinaryDDS" in component_names:
            components.append(vl.BinaryDDS(use_fbx))
        elif "BinaryDDK" in component_names:
            assert (
                "AstrometryEcliptic" in component_names
                or "AstrometryEquatorial" in component_names
            ), "`AstrometryEcliptic` or `AstrometryEquatorial` must be present in the model when `BinaryDDK` is used. Please check the par file."
            ecliptic_coords = "AstrometryEcliptic" in component_names
            components.append(vl.BinaryDDK(use_fbx, ecliptic_coords))
        else:
            raise NotImplementedError(
                f"BINARY {model.BINARY.value} not (yet?) implemented."
            )

    if "FD" in component_names:
        components.append(vl.FrequencyDependent())

    if "FDJump" in component_names:
        assert model["FDJUMPLOG"].value, "Only 'FDJUMPLOG Y' is supported currently."

        fdjumps = [
            model[fdj]
            for fdj in model.components["FDJump"].fdjumps
            if model[fdj].quantity is not None
        ]
        mask0 = read_mask(toas, fdjumps)
        exponents = [model.components["FDJump"].get_fd_index(fd.name) for fd in fdjumps]

        components.append(
            vl.FrequencyDependentJump(
                vl.BitMatrix(mask0),
                vl.Vector[vl.UInt](exponents),
            )
        )

    if "WaveX" in component_names:
        components.append(vl.WaveX())
    elif "PLRedNoiseGP" in component_names:
        components.append(vl.PowerlawRedNoiseGP(int(model["TNREDC"].value)))

    if "Spindown" in component_names:
        components.append(vl.Spindown())

    if "Glitch" in component_names:
        components.append(vl.Glitch())

    if "PhaseOffset" in component_names:
        components.append(vl.PhaseOffset())

    if "PhaseJump" in component_names:
        masks0 = read_mask(toas, model.get_jump_param_objects())
        phase_jump = (
            vl.ExclusivePhaseJump(jl.Vector[jl.UInt](get_exclusive_mask(masks0)))
            if is_exclusive_mask(masks0)
            else vl.PhaseJump(jl.BitMatrix(masks0))
        )
        components.append(phase_jump)

    if "ScaleToaError" in component_names:
        efac_mask0 = read_mask(
            toas, [model[ef] for ef in model.EFACs if model.EFACs[ef][0] is not None]
        )
        equad_mask0 = read_mask(
            toas, [model[eq] for eq in model.EQUADs if model.EQUADs[eq][0] is not None]
        )

        assert len(efac_mask0) == 0 or is_exclusive_mask(
            efac_mask0
        ), "Non-exclusive EFAC masks are not supported. Check the par file for overlapping EFACs."
        assert len(equad_mask0) == 0 or is_exclusive_mask(
            equad_mask0
        ), "Non-exclusive EQUAD masks are not supported. Check the par file for overlapping EQUADs."

        efac_mask = (
            jl.Vector[jl.UInt](get_exclusive_mask(efac_mask0))
            if len(efac_mask0) > 0
            else jl.Vector[jl.UInt](np.zeros(len(toas)))
        )
        equad_mask = (
            jl.Vector[jl.UInt](get_exclusive_mask(equad_mask0))
            if len(equad_mask0) > 0
            else jl.Vector[jl.UInt](np.zeros(len(toas)))
        )

        components.append(vl.MeasurementNoise(efac_mask, equad_mask))

    if "ScaleDmError" in component_names:
        dmefac_mask0 = read_mask(
            toas,
            [
                model[dmef]
                for dmef in model.DMEFACs
                if model.DMEFACs[dmef][0] is not None
            ],
        )
        dmequad_mask0 = read_mask(
            toas,
            [
                model[dmeq]
                for dmeq in model.DMEQUADs
                if model.DMEQUADs[dmeq][0] is not None
            ],
        )

        assert len(dmefac_mask0) == 0 or is_exclusive_mask(
            dmefac_mask0
        ), "Non-exclusive DMEFAC masks are not supported. Check the par file for overlapping DMEFACs."
        assert len(dmequad_mask0) == 0 or is_exclusive_mask(
            dmequad_mask0
        ), "Non-exclusive DMEQUAD masks are not supported. Check the par file for overlapping DMEQUADs."

        dmefac_mask = (
            jl.Vector[jl.UInt](get_exclusive_mask(dmefac_mask0))
            if len(dmefac_mask0) > 0
            else jl.Vector[jl.UInt](np.zeros(len(toas)))
        )
        dmequad_mask = (
            jl.Vector[jl.UInt](get_exclusive_mask(dmequad_mask0))
            if len(dmequad_mask0) > 0
            else jl.Vector[jl.UInt](np.zeros(len(toas)))
        )

        components.append(vl.DispersionMeasurementNoise(dmefac_mask, dmequad_mask))

    return jl.Tuple(components)


def fix_params(model: TimingModel) -> None:
    """Fix the parameters of a `PINT` `TimingModel` to make it `Vela`-friendly.

    It does the following.
        1. Ensures that all the `EPOCH`s are set.
        2. Ensures that `PHOFF` is included and is free.
        3. Converts `H4` to `STIGMA`
        4. Sets the unset parameter values to 0 where possible.
    """

    assert model["PEPOCH"].value is not None, "PEPOCH is not given in the par file."

    for param in model.params:
        if (
            param.endswith("EPOCH")
            and isinstance(model[param], MJDParameter)
            and model[param].value is None
        ):
            model[param].quantity = model["PEPOCH"].quantity

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
    toas: TOAs,
    ecorr_toa_ranges: List[Tuple[int, int]],
    ecorr_indices: List[int],
):
    """Construct a `Vela.Kernel` object. It may be a white noise kernel or
    an ECORR kernel. Time-correlated noise kernels are not yet suppoted."""
    if not model.has_correlated_errors:
        return vl.WhiteNoiseKernel()
    elif not model.has_time_correlated_errors:
        ecorr_mask0 = read_mask(
            toas, [model[ef] for ef in model.EFACs if model.EFACs[ef][0] is not None]
        )

        assert len(ecorr_mask0) == 0 or is_exclusive_mask(
            ecorr_mask0
        ), "Non-exclusive ECORRs are not supported. Check the par file for overlapping ECORRs."

        ecorr_groups = vl.Vector(
            [
                vl.EcorrGroup(start, stop, index)
                for (start, stop), index in zip(ecorr_toa_ranges, ecorr_indices)
            ]
        )
        return vl.EcorrKernel(ecorr_groups)

    raise NotImplementedError("Time-correlated noise kernels are not yet implemented.")


def fix_red_noise_components(model: TimingModel, toas: TOAs):
    f1 = 1 / toas.get_Tspan()
    epoch = model["PEPOCH"].quantity

    if "PLRedNoise" in model.components:
        plred_gp = PLRedNoiseGP(model.components["PLRedNoise"], f1, epoch)
        model.remove_component("PLRedNoise")
        model.add_component(plred_gp)

    if "PLDMNoise" in model.components:
        pldm_gp = PLDMNoiseGP(model.components["PLDMNoise"], f1, epoch)
        model.remove_component("PLDMNoise")
        model.add_component(pldm_gp)

    if "PLChromNoise" in model.components:
        pldm_chrom = PLChromNoiseGP(model.components["PLChromNoise"], f1, epoch)
        model.remove_component("PLChromNoise")
        model.add_component(pldm_chrom)


def pint_model_to_vela(
    model: TimingModel,
    toas: TOAs,
    cheat_prior_scale: float,
    custom_prior_dists: dict,
    ecorr_toa_ranges: Optional[List[Tuple[int, int]]] = None,
    ecorr_indices: Optional[List[Tuple[int]]] = None,
):
    """Construct a `Vela.TimingModel` from a `PINT` `TimingModel`."""

    epoch_mjd = float(model["PEPOCH"].value)

    toas.compute_pulse_numbers(model)

    fix_params(model)

    fix_red_noise_components(model, toas)

    pulsar_name = model["PSR"].value if model["PSR"].value is not None else ""

    components = pint_components_to_vela(model, toas)

    single_params, multi_params = pint_parameters_to_vela(model)
    param_handler = vl.ParamHandler(single_params, multi_params)

    free_params = vl.get_free_param_names(param_handler)

    priors = get_default_priors(
        model, free_params, epoch_mjd, cheat_prior_scale, custom_prior_dists
    )

    tzr_toa = model.get_TZR_toa(toas)
    tzr_toa.compute_pulse_numbers(model)
    tzr_toa = pint_toa_to_vela(tzr_toa, -1, epoch_mjd)

    kernel = get_kernel(model, toas, ecorr_toa_ranges, ecorr_indices)

    return vl.TimingModel(
        pulsar_name,
        model["EPHEM"].value,
        model["CLOCK"].value,
        model["UNITS"].value,
        vl.time(epoch_mjd * day_to_s),
        components,
        kernel,
        param_handler,
        tzr_toa,
        priors,
    )
