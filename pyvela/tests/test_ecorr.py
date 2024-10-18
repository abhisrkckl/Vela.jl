from io import StringIO
import os

import numpy as np
import pytest
from pint.models import get_model, get_model_and_toas

from pyvela.ecorr import ecorr_sort
from pyvela.model import get_kernel
from pyvela.spnta import SPNTA

datadir = os.path.dirname(os.path.realpath(__file__)) + "/datafiles"


@pytest.fixture(scope="module")
def data_B1855p09():
    parfile = f"{datadir}/B1855+09_NANOGrav_9yv1.par"
    timfile = f"{datadir}/B1855+09_NANOGrav_9yv1.tim"
    mp, tp = get_model_and_toas(parfile, timfile)
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


def test_plrednoise_error():
    par = """
        PSR                                  SIM2
        EPHEM                               DE440
        CLOCK                        TT(BIPM2019)
        UNITS                                 TDB
        START              53000.9999999566961806
        FINISH             56985.0000000464502315
        DILATEFREQ                              N
        DMDATA                                  N
        NTOA                                 2000
        CHI2                   1999.9595338710205
        CHI2R                  1.0039957499352512
        TRES                1.5197086404438400383
        RAJ                      5:00:00.00000104 1 0.00000309788366042071
        DECJ                    14:59:59.99982230 1 0.00026526190479642293
        PMRA                                  0.0
        PMDEC                                 0.0
        PX                                    0.0
        F0                  100.00000000000000391 1 6.1147023828669208493e-14
        F1              -9.999985179158219319e-16 1 1.3675433277483732467e-21
        F2                                    0.0
        PEPOCH             55000.0000000000000000
        EFAC            tel gbt        1.2663639225043157 1 0.021424340920438052
        ECORR           tel gbt        0.8409712702910475 1 0.049670237117127584
        PLANET_SHAPIRO                          N
        DM                                   15.0
        DM1                                   0.0
        TZRMJD             55000.0000000000000000
        TZRSITE                               gbt
        TZRFRQ                             1400.0
        PHOFF              4.2709307823244897e-07 1 9.1199685066235e-06
        TNREDAMP                              -14
        TNREDGAM                              3.5
        TNREDC                                 15
    """
    m = get_model(StringIO(par))

    with pytest.raises(NotImplementedError):
        get_kernel(m, None, [], [])
