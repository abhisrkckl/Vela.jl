from typing import List
import numpy as np

from pint.models import TimingModel
from pint.models.parameter import Parameter
from pint.toa import TOAs

from .vela import vl, jl


def read_mask(toas: TOAs, params: List[Parameter]):
    masks = []
    for param in params:
        mask = np.repeat(False, len(toas))
        mask[param.select_toa_mask(toas)] = True
        assert any(mask), f"{param.name} has no TOAs!"
        masks.append(mask)
    return np.array(masks)


def is_exclusive_mask(mask: np.ndarray):
    return all(map(lambda x: x in [0, 1], mask.sum(axis=0)))


def get_exclusive_mask(mask: np.ndarray):
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

    if "FD" in component_names:
        components.append(vl.FrequencyDependent())

    if "WaveX" in component_names:
        components.append(vl.WaveX())

    if "Spindown" in component_names:
        components.append(vl.Spindown())

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

    components = jl.Tuple(components)

    return components
