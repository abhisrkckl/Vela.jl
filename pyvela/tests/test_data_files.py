import os
import numpy as np
import pytest
from juliacall import Main as jl
from pint.models import get_model_and_toas

from pyvela import SPNTA
from pyvela.model import fix_red_noise_components
from pyvela.vela import vl
from pyvela.parameters import fdjump_rx

jl.seval("using BenchmarkTools")
jl.seval("get_alloc(func, args...) = @ballocated(($func)(($args)...))")

datadir = os.path.dirname(os.path.realpath(__file__)) + "/datafiles"

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
    "sim3.gp",
    "sim4",
    "sim4.gp",
    "sim_fd",
    "sim_fdjump",
    "sim6",
    "sim6.gp",
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
def spnta(request):
    dataset = request.param

    custom_prior_dict={
        "PHOFF": jl.Uniform(-0.5, 0.5),
        "EFAC": jl.Uniform(0.5, 2.5),
    }

    spnta = SPNTA(
        f"{datadir}/{dataset}.par",
        f"{datadir}/{dataset}.tim",
        custom_prior_dict=custom_prior_dict
    )
    params = spnta.maxlike_params

    m, t = get_model_and_toas(f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim")

    if (
        len(
            {"PLRedNoise", "PLDMNoise", "PLChromNoise"}.intersection(
                m.components.keys()
            )
        )
        > 0
    ):
        fix_red_noise_components(m, t)

    return spnta, params, m, t


def test_data(spnta):
    spnta: SPNTA
    spnta, params, m, t = spnta

    assert len(spnta.toas) == len(t)

    assert len(spnta.model.components) <= len(m.components)

    pnames = spnta.param_names
    if "H4" not in m.free_params:
        assert set(pnames) == set(m.free_params).union({"PHOFF"})
    else:
        assert set(pnames) == set(m.free_params).union({"PHOFF", "STIGMA"}).difference(
            {"H4"}
        )

    assert len(params) == len(pnames)

    prnames = [str(vl.param_name(pr)) for pr in spnta.model.priors]
    assert len(spnta.model.priors) == len(pnames)
    assert all(
        [pn.startswith(prn) or fdjump_rx.match(pn) for pn, prn in zip(pnames, prnames)]
    )


def test_chi2(spnta):
    spnta, params, m, t = spnta
    calc_chi2 = vl.get_chi2_func(spnta.model, spnta.toas)

    if (
        len(
            {"PLRedNoiseGP", "PLDMNoiseGP", "PLChromNoiseGP"}.intersection(
                m.components.keys()
            )
        )
        == 0
    ):
        assert ("PHOFF" not in m) or (
            calc_chi2(params) / len(spnta.toas) / (1 + int(t.is_wideband())) < 1.2
        )


def test_likelihood(spnta):
    spnta, params, _, _ = spnta
    calc_lnlike = vl.get_lnlike_func(spnta.model, spnta.toas)
    assert np.isfinite(calc_lnlike(params))


def test_prior(spnta):
    spnta, params, _, _ = spnta
    calc_lnprior = vl.get_lnprior_func(spnta.model)
    assert np.isfinite(calc_lnprior(spnta.model.param_handler._default_params_tuple))
    assert calc_lnprior(params) == calc_lnprior(spnta.model.param_handler._default_params_tuple)


def test_alloc(spnta):
    spnta, _, _, _ = spnta
    params = spnta.model.param_handler._default_params_tuple
    calc_lnlike = vl.get_lnlike_serial_func(spnta.model, spnta.toas)
    assert jl.get_alloc(calc_lnlike, params) == 0


def test_posterior(spnta):
    spnta, _, _, _ = spnta
    parv = np.array(
        vl.read_param_values_to_vector(
            spnta.model.param_handler, spnta.model.param_handler._default_params_tuple
        )
    )
    maxlike_params_v = np.array([parv, parv, parv, parv])

    lnpost = vl.get_lnpost_func(spnta.model, spnta.toas, True)

    lnpvals = lnpost(maxlike_params_v)

    assert np.allclose(lnpvals, lnpvals[0])
