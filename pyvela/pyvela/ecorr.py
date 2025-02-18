from typing import List, Tuple

import astropy.units as u
import numpy as np
from pint.models import TimingModel
from pint.models.noise_model import get_ecorr_nweights
from pint.toa import TOAs


def ecorr_weight_indices(model: TimingModel, toas: TOAs) -> np.ndarray:
    """Returns an array containing the ECORR parameter indices corresponding to each ECORR group.
    The indices start from 1 to match Julia indexing. 0 represents no ECORR group affiliation.
    """
    ecorrs = model.components["EcorrNoise"].get_ecorrs()

    ts = (toas.table["tdbld"].quantity * u.day).to_value(u.s)
    nweights = [get_ecorr_nweights(ts[ec.select_toa_mask(toas)]) for ec in ecorrs]
    nc = sum(nweights)

    weight_indices = np.zeros(nc, dtype=int)
    nctot = 0
    for ec, nn in zip(ecorrs, nweights):
        weight_indices[nctot : nn + nctot] = ec.index
        nctot += nn

    # weight_indices = np.insert(weight_indices, 0, 0)

    return weight_indices


def ecorr_sort(
    model: TimingModel, toas: TOAs
) -> Tuple[TOAs, List[Tuple[int]], np.ndarray]:
    """Returns a `TOAs` object that is sorted according to ECORR groups.
    The TOAs without any ECORR group affiliation appear at the beginning.
    Also returns a list containing the start and end indices of each ECORR
    group and an array containing the ECORR parameter indices for each group."""

    assert (
        "EcorrNoise" in model.components
    ), "`ecorr_sort()` was called without any ECORRs in the model."

    ecorr_masks = model.components["EcorrNoise"].get_noise_basis(toas).T.astype(bool)

    toa_indices = np.arange(len(toas))

    ecorr_sort_mask = []
    toa_ranges = []

    weight_indices = ecorr_weight_indices(model, toas)

    # TOAs which do not belong to any ECORR group
    no_ecorr_mask = np.logical_not(np.any(ecorr_masks, axis=0))
    if np.any(no_ecorr_mask):
        ecorr_sort_mask.extend(toa_indices[no_ecorr_mask])
        toa_ranges.append((1, int(sum(no_ecorr_mask))))
        weight_indices = np.insert(weight_indices, 0, 0)

    for ecmask in ecorr_masks:
        ecorr_sort_mask.extend(toa_indices[ecmask])
        toa_ranges.append(
            (len(ecorr_sort_mask) - int(sum(ecmask)) + 1, len(ecorr_sort_mask))
        )

    assert len(ecorr_sort_mask) == len(
        toas
    ), "Shape of the ECORR mask is inconsistent with the number of TOAs. This is a bug."
    assert len(toa_ranges) == len(
        weight_indices
    ), "Number of ECORR groups is inconsistent with the weight indices. This is a bug."

    toas_sorted = toas[ecorr_sort_mask]
    toas_sorted.table["index"] = np.arange(len(toas))

    return toas_sorted, toa_ranges, weight_indices
