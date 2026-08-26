import numpy as np
from pint.models import TimingModel
from pint.toa import TOAs
from pint.utils import get_prefix_timeranges


def get_dmx_mask(
    model: TimingModel, toas: TOAs, param_prefix: str = "DMX_"
) -> np.ndarray:
    """Get a Vela-compatible DMX/CMX mask given a timing model.
    The output is an ndarray containing the DMX/CMX index for each TOA.
    Throws an exception if overlapping DMX ranges are found."""

    indices, r1s, r2s = get_prefix_timeranges(model, param_prefix)
    r1s = r1s.mjd
    r2s = r2s.mjd

    assert all(r1s[1:] > r2s[:-1]), f"Overlapping {param_prefix} ranges found!"

    output_mask = np.repeat(0, len(toas))
    mjds = toas.get_mjds().value
    for ii, r1, r2 in zip(indices, r1s, r2s):
        bool_mask = np.logical_and(mjds >= r1, mjds < r2)
        assert np.any(bool_mask), f"No TOAs found for {param_prefix}{ii:4d}."
        output_mask[bool_mask] = ii
    assert np.all(
        output_mask > 0
    ), f"Some TOAs don't belong to any {param_prefix} ranges."

    return output_mask
