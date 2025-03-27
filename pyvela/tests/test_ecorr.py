import os
from io import StringIO

import numpy as np
import astropy.units as u
import pytest
from pint.models import get_model, get_model_and_toas
from pint.simulation import make_fake_toas_uniform

from pyvela.ecorr import ecorr_sort
from pyvela.model import get_kernel
from pyvela.spnta import SPNTA

datadir = f"{os.path.dirname(os.path.realpath(__file__))}/datafiles"


@pytest.fixture(scope="module")
def data_B1855p09():
    parfile = f"{datadir}/B1855+09_NANOGrav_9yv1.par"
    timfile = f"{datadir}/B1855+09_NANOGrav_9yv1.tim"
    mp, tp = get_model_and_toas(parfile, timfile, usepickle=True)
    spnta = SPNTA(parfile, timfile)

    return mp, tp, spnta


def test_ecorr_sort(data_B1855p09):
    m, t, _ = data_B1855p09

    t_sorted, ranges, indices = ecorr_sort(m, t)

    assert len(t_sorted) == len(t)
    assert np.all(np.diff(t_sorted.table["index"]) == 1)

    assert len(ranges) == len(indices)

    for i in range(1, len(ranges)):
        assert ranges[i][0] == ranges[i - 1][1] + 1

    assert (
        len(np.unique(m.components["EcorrNoise"].get_noise_weights(t)))
        == len(np.unique(indices)) - 1
    )


def test_read_data(data_B1855p09):
    spnta: SPNTA
    mp, tp, spnta = data_B1855p09

    t_sorted, ranges, indices = ecorr_sort(mp, tp)

    assert len(t_sorted) == len(spnta.toas)
    assert len(spnta.model.kernel.ecorr_groups) == len(indices)
    assert spnta.model.kernel.ecorr_groups[0].index == 0
