from pint2vela.ecorr import ecorr_sort
from pint2vela import read_model_and_toas
from pint.models import get_model_and_toas

import numpy as np
import pytest


@pytest.fixture(scope="module")
def data_B1855p09():
    parfile = "datafiles/B1855+09_NANOGrav_9yv1.par"
    timfile = "datafiles/B1855+09_NANOGrav_9yv1.tim"
    mp, tp = get_model_and_toas(parfile, timfile)
    mv, tv = read_model_and_toas(parfile, timfile)

    return mp, tp, mv, tv


def test_ecorr_sort(data_B1855p09):
    m, t, _, _ = data_B1855p09

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
    mp, tp, mv, tv = data_B1855p09

    t_sorted, ranges, indices = ecorr_sort(mp, tp)

    assert len(t_sorted) == len(tv)
    assert len(mv.kernel.ecorr_groups) == len(indices)
    assert mv.kernel.ecorr_groups[0].index == 0
