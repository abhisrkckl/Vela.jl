import numpy as np
import pytest

from pint.models import get_model_and_toas
from pint2vela import read_model_and_toas
from pint2vela.vela import vl

datasets = [
    "NGC6440E",
    "pure_rotator",
    "sim1",
    "sim_jump",
    "sim_jump_ex",
    "sim_sw",
    "sim_fdjumpdm",
    "sim_cm",
    "sim3",
    "sim4",
]


@pytest.fixture(params=datasets, scope="module")
def model_and_toas(request):
    dataset = request.param
    mv, tv = read_model_and_toas(f"datafiles/{dataset}.par", f"datafiles/{dataset}.tim")
    params = vl.read_param_values_to_vector(
        mv.param_handler, mv.param_handler._default_params_tuple
    )

    m, t = get_model_and_toas(f"datafiles/{dataset}.par", f"datafiles/{dataset}.tim")

    return mv, tv, params, m, t


def test_data(model_and_toas):
    mv, tv, params, m, t = model_and_toas

    assert len(tv) == len(t)

    assert len(mv.components) <= len(m.components)

    pnames = vl.get_free_param_names(mv.param_handler)
    assert set(pnames) == set(m.free_params)

    assert len(params) == len(pnames)

    prnames = [str(vl.param_name(pr)) for pr in mv.priors]
    assert len(mv.priors) == len(pnames)
    assert all([pn.startswith(prn) for pn, prn in zip(pnames, prnames)])


def test_chi2(model_and_toas):
    mv, tv, params, _, _ = model_and_toas
    calc_chi2 = vl.get_chi2_func(mv, tv)
    assert calc_chi2(params) / len(tv) < 1.1


def test_likelihood(model_and_toas):
    mv, tv, params, _, _ = model_and_toas
    calc_lnlike = vl.get_lnlike_func(mv, tv)
    assert np.isfinite(calc_lnlike(params))


def test_prior(model_and_toas):
    mv, _, params, _, _ = model_and_toas
    calc_lnprior = vl.get_lnprior_func(mv)
    assert np.isfinite(calc_lnprior(mv.param_handler._default_params_tuple))
    assert calc_lnprior(params) == calc_lnprior(mv.param_handler._default_params_tuple)
