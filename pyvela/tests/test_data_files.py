import os
import time
from typing import Tuple

import numpy as np
import pytest
from pint.models import get_model_and_toas, TimingModel
from pint.toa import TOAs

from pyvela.model import fix_params, fix_red_noise_components
from pyvela.parameters import fdjump_rx
from pyvela.spnta import SPNTA, convert_model_and_toas
from pyvela.vela import vl

# jl.seval("using BenchmarkTools")
# jl.seval("get_alloc(func, args...) = @ballocated(($func)(($args)...))")

datadir = f"{os.path.dirname(os.path.realpath(__file__))}/datafiles"

datasets = [
    "sim_1",
    "sim_sw.wb",
    "sim3.gp",
    "sim3",
    "sim_fdjump",
    "sim_ddk",
    "sim_glitch",
    "sim_dmx",
    "sim_cmx",
    "sim_ell1k",
    "sim_dmgp_wb",
    "J0613-0200.sim",
    "J1208-5936.sim",
    "J2302+4442.sim",
    "J1227-6208.sim",
]


@pytest.fixture(params=datasets[:3], scope="module")
def model_and_toas(request):
    dataset = request.param

    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"

    m, t = get_model_and_toas(
        f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim", usepickle=True
    )

    custom_priors = f"{datadir}/custom_priors.json"

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
        fix_params(m, t)
        fix_red_noise_components(m, t)

    return spnta, m, t


@pytest.mark.parametrize("dataset", datasets[3:])
def test_read_data(dataset):
    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"
    m, t = get_model_and_toas(parfile, timfile, planets=True)
    model, toas = convert_model_and_toas(
        m, t, m.get_params_of_component_type("NoiseComponent"), False
    )
    assert len(toas) == len(t)
    assert len(model.components) <= len(m.components)


def test_data(model_and_toas: Tuple[SPNTA, TimingModel, TOAs]):
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

    assert len(spnta.default_params) == len(pnames)

    prnames = [str(vl.param_name(pr)) for pr in spnta.model.priors]
    assert len(spnta.model.priors) == len(pnames)
    assert all(
        pn.startswith(prn) or fdjump_rx.match(pn) for pn, prn in zip(pnames, prnames)
    )

    assert all(np.isfinite(spnta.rescale_samples(spnta.default_params)))

    assert spnta.wideband == t.is_wideband()

    assert all(np.isfinite(spnta.mjds)) and len(spnta.mjds) == len(t)

    assert all(np.isfinite(spnta.time_residuals(spnta.default_params)))

    assert all(np.isfinite(spnta.scaled_toa_unceritainties(spnta.default_params)))

    if spnta.wideband:
        assert all(np.isfinite(spnta.dm_residuals(spnta.default_params)))
        assert all(np.isfinite(spnta.scaled_dm_unceritainties(spnta.default_params)))

    assert all(np.isfinite(spnta.model_dm(spnta.default_params)))

    assert (
        all(len(label) > 0 for label in spnta.param_labels)
        and len(spnta.param_labels) == spnta.ndim
    )

    assert len(spnta.param_units) == spnta.ndim


def test_chi2(model_and_toas: Tuple[SPNTA, TimingModel, TOAs]):
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
        assert (
            calc_chi2(spnta.default_params)
            / len(spnta.toas)
            / (1 + int(t.is_wideband()))
            < 1.5
        )


def test_likelihood(model_and_toas):
    spnta: SPNTA
    spnta, _, _ = model_and_toas
    calc_lnlike = vl.get_lnlike_func(spnta.model, spnta.toas)
    assert np.isfinite(calc_lnlike(spnta.default_params))


def test_prior(model_and_toas):
    spnta: SPNTA
    spnta, _, _ = model_and_toas
    calc_lnprior = vl.get_lnprior_func(spnta.model)
    assert np.isfinite(calc_lnprior(spnta.model.param_handler._default_params_tuple))
    assert calc_lnprior(spnta.default_params) == calc_lnprior(
        spnta.model.param_handler._default_params_tuple
    )


# def test_alloc(model_and_toas):
#     spnta: SPNTA
#     spnta, _, _ = model_and_toas
#     params = spnta.model.param_handler._default_params_tuple
#     calc_lnlike = vl.get_lnlike_serial_func(spnta.model, spnta.toas)
#     assert jl.get_alloc(calc_lnlike, params) == 0


def test_posterior(model_and_toas):
    spnta: SPNTA
    spnta, _, _ = model_and_toas
    parv = spnta.default_params
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

    spnta2 = SPNTA.load_jlso(jlsoname, spnta.model_pint.name)
    assert len(spnta2.toas) == len(spnta.toas)
    assert set(spnta2.param_names) == set(spnta.param_names)

    os.unlink(jlsoname)


def test_gp_model_conversion():
    dataset = "sim3.gp"
    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"
    m, t = get_model_and_toas(parfile, timfile, planets=True)

    assert all(c in m.components for c in ["PLRedNoise", "PLDMNoise", "PLChromNoise"])

    fix_params(m, t)

    fix_red_noise_components(m, t)

    assert all(
        c not in m.components for c in ["PLRedNoise", "PLDMNoise", "PLChromNoise"]
    )
    assert all(
        c in m.components for c in ["PLRedNoiseGP", "PLDMNoiseGP", "PLChromNoiseGP"]
    )


@pytest.mark.parametrize("dataset", ["sim3.gp", "sim7.gp"])
def test_gp_model_marg(dataset):
    parfile, timfile = f"{datadir}/{dataset}.par", f"{datadir}/{dataset}.tim"

    spnta1 = SPNTA(parfile, timfile, marginalize_gp_noise=True)
    assert len(spnta1.param_names) == len(spnta1.model_pint.free_params)

    assert np.isfinite(spnta1.lnpost(spnta1.default_params))
    assert np.isclose(
        spnta1.lnpost(spnta1.default_params),
        spnta1.lnpost_vectorized(np.array([spnta1.default_params]))[0],
    )
