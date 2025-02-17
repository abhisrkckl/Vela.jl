from typing import List

import numpy as np
from astropy.time import Time
from pint.models import TimingModel
from pint.toa import TOAs
from pint.utils import find_prefix_bytime


def get_dmx_mask(
    model: TimingModel, toas: TOAs, param_prefix: str = "DMX_"
) -> np.ndarray:
    mask: List[int] = []
    for ii in range(len(toas)):
        t: Time = toas.table["mjd"][ii]
        idxs = find_prefix_bytime(model, param_prefix, t)
        if np.isscalar(idxs):
            mask.append(int(idxs))
        elif len(idxs) == 0:
            mask.append(0)
        else:
            raise ValueError(
                f"Multiple {param_prefix} ranges found for TOA number {ii} ({str(t)})."
            )
    assert len(mask) == len(toas)
    return np.array(mask, dtype=int)
