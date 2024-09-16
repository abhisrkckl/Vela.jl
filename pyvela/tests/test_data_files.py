import numpy as np
import pytest
from juliacall import Main as jl
from pint.models import get_model_and_toas

from pyvela import read_model_and_toas
from pyvela.vela import vl
from pyvela.convert_parameters import fdjump_rx

jl.seval("using BenchmarkTools")
jl.seval("get_alloc(func, args...) = @ballocated(($func)(($args)...))")

datasets = [
    "NGC6440E",
    "pure_rotator",
    "sim1",
    "sim_dmwn",
    "sim_jump",
    "sim_jump_ex",
    "sim_sw",
    "sim_sw.wb",
    "sim_fdjumpdm",
    "sim_dmjump",
    "sim_cm",
    "sim2",
    "sim3",
    "sim4",
    "sim_fd",
    "sim_fdjump",
    "sim6",
    "sim_dd",
    "sim_ddk",
    "sim_glitch",
    "J0613-0200.InPTA.NB",
    "J1857+0943.InPTA.NB",
    "J0613-0200.sim",
    "J1856-3754.sim",
    "J1802-2124.sim",
    "J0955-6150.sim",
    "J1208-5936.sim",
    "J2302+4442.sim",
    "J1227-6208.sim",
]


@pytest.fixture(params=datasets, scope="module")
def model_and_toas(request):
    dataset = request.param
    mv, tv = read_model_and_toas(
        f"datafiles/{dataset}.par",
        f"datafiles/{dataset}.tim",
        custom_prior_dicts={
            "PHOFF": jl.Uniform(-0.5, 0.5),
            "EFAC": jl.Uniform(0.5, 2.5),
        },
    )
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
    if "H4" not in m.free_params:
        assert set(pnames) == set(m.free_params).union({"PHOFF"})
    else:
        assert set(pnames) == set(m.free_params).union({"PHOFF", "STIGMA"}).difference(
            {"H4"}
        )

    assert len(params) == len(pnames)

    prnames = [str(vl.param_name(pr)) for pr in mv.priors]
    assert len(mv.priors) == len(pnames)
    assert all(
        [pn.startswith(prn) or fdjump_rx.match(pn) for pn, prn in zip(pnames, prnames)]
    )


def test_chi2(model_and_toas):
    mv, tv, params, m, t = model_and_toas
    calc_chi2 = vl.get_chi2_func(mv, tv)
    assert ("PHOFF" not in m) or (
        calc_chi2(params) / len(tv) / (1 + int(t.is_wideband())) < 1.2
    )


def test_likelihood(model_and_toas):
    mv, tv, params, _, _ = model_and_toas
    calc_lnlike = vl.get_lnlike_func(mv, tv)
    assert np.isfinite(calc_lnlike(params))


def test_prior(model_and_toas):
    mv, _, params, _, _ = model_and_toas
    calc_lnprior = vl.get_lnprior_func(mv)
    assert np.isfinite(calc_lnprior(mv.param_handler._default_params_tuple))
    assert calc_lnprior(params) == calc_lnprior(mv.param_handler._default_params_tuple)


def test_alloc(model_and_toas):
    mv, tv, _, _, _ = model_and_toas
    params = mv.param_handler._default_params_tuple
    calc_lnlike = vl.get_lnlike_serial_func(mv, tv)
    assert jl.get_alloc(calc_lnlike, params) == 0