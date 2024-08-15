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

    if "TroposphereDelay" in component_names and model.CORRECT_TROPOSPHERE.value:
        components.append(vl.Troposphere())

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

    if "ChromaticCM" in component_names:
        components.append(vl.ChromaticTaylor())

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
        efac_mask0 = read_mask(toas, [model[ef] for ef in model.EFACs])
        equad_mask0 = read_mask(toas, [model[eq] for eq in model.EQUADs])

        assert is_exclusive_mask(efac_mask0)
        assert is_exclusive_mask(equad_mask0)

        efac_mask = jl.Vector[jl.UInt](get_exclusive_mask(efac_mask0))
        equad_mask = jl.Vector[jl.UInt](get_exclusive_mask(equad_mask0))

        components.append(vl.MeasurementNoise(efac_mask, equad_mask))

    components = jl.Tuple(components)

    return components
