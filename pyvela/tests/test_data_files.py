import os
import time

import numpy as np
import pytest
from juliacall import Main as jl
from pint.models import get_model_and_toas

from pyvela.model import fix_red_noise_components
from pyvela.parameters import fdjump_rx
from pyvela.spnta import SPNTA
from pyvela.vela import vl

jl.seval("using BenchmarkTools")
jl.seval("get_alloc(func, args...) = @ballocated(($func)(($args)...))")

datadir = os.path.dirname(os.path.realpath(__file__)) + "/datafiles"

datasets = [
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
    "sim_ddk",
    "sim_glitch",
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

    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"

    m, t = get_model_and_toas(f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim")

    custom_priors = f"{datadir}/custom_priors.json" if "EFAC1" in m.free_params else {}

    spnta = SPNTA(
        parfile,
        timfile,
        custom_priors=custom_priors,
    )

    if (
        len(
            {"PLRedNoise", "PLDMNoise", "PLChromNoise"}.intersection(
                m.components.keys()
            )
        )
        > 0
    ):
        fix_red_noise_components(m, t)

    return spnta, m, t


def test_data(model_and_toas):
    spnta: SPNTA
    spnta, m, t = model_and_toas

    assert len(spnta.toas) == len(t)

    assert len(spnta.model.components) <= len(m.components)

    pnames = vl.get_free_param_names(spnta.model.param_handler)
    if "H4" not in m.free_params:
        assert set(pnames) == set(m.free_params).union({"PHOFF"})
    else:
        assert set(pnames) == set(m.free_params).union({"PHOFF", "STIGMA"}).difference(
            {"H4"}
        )

    assert len(spnta.maxlike_params) == len(pnames)

    prnames = [str(vl.param_name(pr)) for pr in spnta.model.priors]
    assert len(spnta.model.priors) == len(pnames)
    assert all(
        [pn.startswith(prn) or fdjump_rx.match(pn) for pn, prn in zip(pnames, prnames)]
    )

    assert all(np.isfinite(spnta.rescale_samples(spnta.maxlike_params)))

    assert spnta.is_wideband() == t.is_wideband()

    assert all(np.isfinite(spnta.get_mjds())) and len(spnta.get_mjds()) == len(t)


def test_chi2(model_and_toas):
    spnta: SPNTA
    spnta, m, t = model_and_toas
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
            calc_chi2(spnta.maxlike_params)
            / len(spnta.toas)
            / (1 + int(t.is_wideband()))
            < 1.2
        )


def test_likelihood(model_and_toas):
    spnta: SPNTA
    spnta, _, _ = model_and_toas
    calc_lnlike = vl.get_lnlike_func(spnta.model, spnta.toas)
    assert np.isfinite(calc_lnlike(spnta.maxlike_params))


def test_prior(model_and_toas):
    spnta: SPNTA
    spnta, _, _ = model_and_toas
    calc_lnprior = vl.get_lnprior_func(spnta.model)
    assert np.isfinite(calc_lnprior(spnta.model.param_handler._default_params_tuple))
    assert calc_lnprior(spnta.maxlike_params) == calc_lnprior(
        spnta.model.param_handler._default_params_tuple
    )


def test_alloc(model_and_toas):
    spnta: SPNTA
    spnta, _, _ = model_and_toas
    params = spnta.model.param_handler._default_params_tuple
    calc_lnlike = vl.get_lnlike_serial_func(spnta.model, spnta.toas)
    assert jl.get_alloc(calc_lnlike, params) == 0


def test_posterior(model_and_toas):
    spnta: SPNTA
    spnta, _, _ = model_and_toas
    parv = spnta.maxlike_params
    maxlike_params_v = np.array([parv, parv, parv, parv])

    lnpost = vl.get_lnpost_func(spnta.model, spnta.toas, True)

    lnpvals = lnpost(maxlike_params_v)

    assert np.allclose(lnpvals, lnpvals[0])


def test_readwrite_jlso(model_and_toas):
    spnta: SPNTA
    spnta, _, _ = model_and_toas

    jlsoname = f"__{spnta.model.pulsar_name}.jlso"

    spnta.save_jlso(jlsoname)
    assert os.path.isfile(jlsoname)

    time.sleep(0.1)

    spnta2 = SPNTA.load_jlso(jlsoname)
    assert len(spnta2.toas) == len(spnta.toas)
    assert set(spnta2.param_names) == set(spnta.param_names)

    os.unlink(jlsoname)
