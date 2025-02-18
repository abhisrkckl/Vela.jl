import os

import pytest
from pint.models import get_model_and_toas
from pint.models.parameter import floatParameter

from pyvela.parameters import get_scale_factor
from pyvela.spnta import SPNTA

datadir = f"{os.path.dirname(os.path.realpath(__file__))}/datafiles"


def test_unknown_scale_factor():
    bla = floatParameter(
        name="BLA",
        value=0.0,
        units="",
        frozen=False,
        tcb2tdb_scale_factor=None,
        convert_tcb2tdb=False,
    )
    with pytest.raises(ValueError):
        get_scale_factor(bla)


def test_open_prior_file():
    parfile, timfile = f"{datadir}/NGC6440E.par", f"{datadir}/NGC6440E.tim"
    with open(f"{datadir}/custom_priors.json", "r") as prior_file:
        spnta1 = SPNTA(parfile, timfile, custom_priors=prior_file)

    m, t = get_model_and_toas(parfile, timfile, planets=True)
    spnta2 = SPNTA.from_pint(m, t, custom_priors=f"{datadir}/custom_priors.json")

    assert all(pr1 == pr2 for pr1, pr2 in zip(spnta1.model.priors, spnta2.model.priors))


def test_dmx_exception():
    parfile, timfile = f"{datadir}/sim_dmx.par", f"{datadir}/sim_dmx.tim"
    m, t = get_model_and_toas(parfile, timfile, planets=True)
    m["DMXR2_0001"].value = m["DMXR2_0002"].value
    with pytest.raises(ValueError):
        SPNTA.from_pint(m, t)
