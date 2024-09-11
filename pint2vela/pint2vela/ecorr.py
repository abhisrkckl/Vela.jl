from pint.models import TimingModel
from pint.models.noise_model import get_ecorr_nweights
from pint.toa import TOAs

import numpy as np
import astropy.units as u


def ecorr_weight_indices(model: TimingModel, toas: TOAs):
    ecorrs = model.components["EcorrNoise"].get_ecorrs()

    ts = (toas.table["tdbld"].quantity * u.day).to(u.s).value
    nweights = [get_ecorr_nweights(ts[ec.select_toa_mask(toas)]) for ec in ecorrs]
    nc = sum(nweights)

    weight_indices = np.zeros(nc, dtype=int)
    nctot = 0
    for ec, nn in zip(ecorrs, nweights):
        weight_indices[nctot : nn + nctot] = ec.index
        nctot += nn

    # weight_indices = np.insert(weight_indices, 0, 0)

    return weight_indices


def ecorr_sort(model: TimingModel, toas: TOAs):
    assert "EcorrNoise" in model.components

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

    assert len(ecorr_sort_mask) == len(toas)
    assert len(toa_ranges) == len(weight_indices)

    toas_sorted = toas[ecorr_sort_mask]
    toas_sorted.table["index"] = np.arange(len(toas))

    return toas_sorted, toa_ranges, weight_indices
