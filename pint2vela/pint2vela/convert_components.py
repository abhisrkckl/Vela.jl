from typing import List

import numpy as np
from pint.models import TimingModel
from pint.models.parameter import maskParameter
from pint.toa import TOAs

from .vela import jl, vl


def read_mask(toas: TOAs, params: List[maskParameter]):
    """Read a TOA mask from a `maskParameter` in a `Vela`-friendly
    representation."""

    masks = []
    for param in params:
        mask = np.repeat(False, len(toas))
        mask[param.select_toa_mask(toas)] = True
        assert any(mask), f"{param.name} has no TOAs!"
        masks.append(mask)
    return np.array(masks)


def is_exclusive_mask(mask: np.ndarray):
    """Check if the mask is exclusive. An exclusive mask is where one TOA
    belongs to only one group.

    For example, `EFAC`s and `EQUAD`s are generally exclusive, whereas `JUMP`s
    are sometimes not.
    """
    return all(map(lambda x: x in [0, 1], mask.sum(axis=0)))


def get_exclusive_mask(mask: np.ndarray):
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
    assert len(result) == mask.shape[1]
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
        components.append(vl.SolarSystem(True, model.PLANET_SHAPIRO.value))
    elif "AstrometryEquatorial" in component_names:
        components.append(vl.SolarSystem(False, model.PLANET_SHAPIRO.value))

    if "SolarWindDispersion" in component_names and not (
        model.NE_SW.value == 0 and model.NE_SW.frozen
    ):
        components.append(vl.SolarWindDispersion())

    if "DispersionDM" in component_names:
        components.append(vl.DispersionTaylor())

    if "DMWaveX" in component_names:
        components.append(vl.DMWaveX())

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

    if "CMWaveX" in component_names:
        components.append(vl.CMWaveX())

    if model.BINARY.value is not None:
        assert (model["PB"].quantity is not None) != (model["FB0"].quantity is not None)
        use_fbx = model["FB0"].quantity is not None
        if "BinaryELL1" in component_names:
            components.append(vl.BinaryELL1(use_fbx))
        elif "BinaryELL1H" in component_names:
            components.append(vl.BinaryELL1H(use_fbx))
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
            )
            ecliptic_coords = "AstrometryEcliptic" in component_names
            components.append(vl.BinaryDDK(use_fbx, ecliptic_coords))
        else:
            raise NotImplementedError(
                f"BINARY {model.BINARY.value} not (yet?) implemented."
            )

    if "FD" in component_names:
        components.append(vl.FrequencyDependent())

    if "WaveX" in component_names:
        components.append(vl.WaveX())

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

        assert len(efac_mask0) == 0 or is_exclusive_mask(efac_mask0)
        assert len(equad_mask0) == 0 or is_exclusive_mask(equad_mask0)

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

        assert len(dmefac_mask0) == 0 or is_exclusive_mask(dmefac_mask0)
        assert len(dmequad_mask0) == 0 or is_exclusive_mask(dmequad_mask0)

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

    components = jl.Tuple(components)

    return components
