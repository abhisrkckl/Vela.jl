import numpy as np
import pytest

from pint2vela import read_model_and_toas
from pint2vela.vela import vl


@pytest.fixture
def data_NGC6440E():
    mv, tv = read_model_and_toas("datafiles/NGC6440E.par", "datafiles/NGC6440E.tim")
    params = vl.read_param_values_to_vector(
        mv.param_handler, mv.param_handler._default_params_tuple
    )
    return mv, tv, params


def test_data(data_NGC6440E):
    mv, tv, params = data_NGC6440E

    assert len(tv) == 62

    assert len(mv.components) == 4

    pnames = vl.get_free_param_names(mv.param_handler)
    assert set(pnames) == {
        "RAJ",
        "DECJ",
        "DM",
        "PHOFF",
        "F0",
        "F1",
    }

    assert len(params) == 6

    prnames = [str(vl.param_name(pr)) for pr in mv.priors]
    assert len(mv.priors) == len(pnames)
    assert all([pn.startswith(prn) for pn, prn in zip(pnames, prnames)])


def test_chi2(data_NGC6440E):
    mv, tv, params = data_NGC6440E
    calc_chi2 = vl.get_chi2_func(mv, tv)
    assert calc_chi2(params) / len(tv) < 1.1


def test_likelihood(data_NGC6440E):
    mv, tv, params = data_NGC6440E
    calc_lnlike = vl.get_lnlike_func(mv, tv)
    assert np.isfinite(calc_lnlike(params))


def test_prior(data_NGC6440E):
    mv, _, params = data_NGC6440E
    calc_lnprior = vl.get_lnprior_func(mv)
    assert np.isfinite(calc_lnprior(mv.param_handler._default_params_tuple))
    assert calc_lnprior(params) == calc_lnprior(mv.param_handler._default_params_tuple)
