import numpy as np

from pint.models import TimingModel
from pint.toa import TOAs

from .vela import vl, jl


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

    if "SolarWindDispersion" in component_names and model.NE_SW.value != 0:
        components.append(vl.SolarWindDispersion(int(model.SWM.value)))

    if "DispersionDM" in component_names:
        components.append(vl.DispersionTaylor())

    if "Spindown" in component_names:
        components.append(vl.Spindown())

    if "PhaseOffset" in component_names:
        components.append(vl.PhaseOffset())

    if "PhaseJump" in component_names:
        masks = []
        for jump in model.get_jump_param_objects():
            mask = np.repeat(False, len(toas))
            mask[jump.select_toa_mask(toas)] = True
            assert any(mask), f"{jump.name} has no TOAs!"
            masks.append(mask)
        masks = jl.BitMatrix(np.array(masks))
        components.append(vl.PhaseJump(masks))

    if "ScaleToaError" in component_names:
        efac_masks = [model[efac].select_toa_mask(toas) for efac in model.EFACs]
        equad_masks = [model[equad].select_toa_mask(toas) for equad in model.EQUADs]

        efac_index_mask = []
        equad_index_mask = []
        for ii in range(len(toas)):
            efac_idxs = np.argwhere([(ii in mask) for mask in efac_masks]).flatten()
            equad_idxs = np.argwhere([(ii in mask) for mask in equad_masks]).flatten()

            assert len(efac_idxs) in [0, 1], "EFAC groups must not overlap."
            assert len(equad_idxs) in [0, 1], "EQUAD groups must not overlap."

            efac_index_mask.append(efac_idxs.item() + 1 if len(efac_idxs) == 1 else 0)
            equad_index_mask.append(
                equad_idxs.item() + 1 if len(equad_idxs) == 1 else 0
            )

        efac_index_mask = jl.Vector[jl.UInt](efac_index_mask)
        equad_index_mask = jl.Vector[jl.UInt](equad_index_mask)

        assert len(set(efac_index_mask)) in [
            len(model.EFACs),
            len(model.EFACs) + 1,
        ], "EFAC without any TOAs found!"
        assert len(set(equad_index_mask)) in [
            len(model.EQUADs),
            len(model.EQUADs) + 1,
        ], "EQUAD without any TOAs found!"

        components.append(vl.MeasurementNoise(efac_index_mask, equad_index_mask))

    components = jl.Tuple(components)

    return components
