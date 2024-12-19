import os
from pint.models.parameter import floatParameter
from pint.models import get_model, get_model_and_toas
from pint.simulation import make_fake_toas_uniform

from pyvela.parameters import get_scale_factor
from pyvela.model import pint_components_to_vela

import pytest
from io import StringIO

from pyvela.spnta import SPNTA

datadir = os.path.dirname(os.path.realpath(__file__)) + "/datafiles"


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


def test_ell1k_exception():
    par = """
        # Based on J0636+5128. Only for testing. Relativistic 
        # effects may not be physically consistent.
        PSRJ            PSRTEST
        ELONG           96.363147518                  1.400e-08
        ELAT            28.24309900                   3.000e-08
        DM              11.108159
        PEPOCH          57277.00
        F0              348.5592316999902             9.000e-13
        F1              -4.1895E-16                   8.000e-20
        POSEPOCH        57277.00
        DMEPOCH         56307.00
        BINARY          ELL1k
        A1              0.0898636            1         6.000e-08
        TASC            57277.01594431       1         1.100e-07
        EPS1            1.4E-6               1         1.000e-05
        EPS2            1.7E-5               1         9.000e-06
        OMDOT           100                  1
        LNEDOT          0                    1
        CLK             TT(BIPM2017)
        EPHEM           DE436
        RM              -7                            1.000e+00
        PX              1.4                           3.000e-01
        FB0             0.00017391195942     1         4.000e-14
        FB1             -7.7E-20             1         3.000e-21
        PMELONG         3.33                          4.000e-02
        PMELAT          -1.44                         1.100e-01
        UNITS           TDB
    """
    model = get_model(StringIO(par))
    fake_toas = make_fake_toas_uniform(
        startMJD=50000, endMJD=51000, ntoas=100, model=model, add_noise=True
    )

    with pytest.raises(NotImplementedError):
        pint_components_to_vela(model, fake_toas)


def test_open_prior_file():
    parfile, timfile = f"{datadir}/NGC6440E.par", f"{datadir}/NGC6440E.tim"
    with open(f"{datadir}/custom_priors.json", "r") as prior_file:
        spnta1 = SPNTA(parfile, timfile, custom_priors=prior_file)

    m, t = get_model_and_toas(parfile, timfile, planets=True)
    spnta2 = SPNTA.from_pint(m, t, custom_priors=f"{datadir}/custom_priors.json")

    assert all(
        [pr1 == pr2 for pr1, pr2 in zip(spnta1.model.priors, spnta2.model.priors)]
    )


def test_dmx_exception():
    parfile, timfile = f"{datadir}/sim_dmx.par", f"{datadir}/sim_dmx.tim"
    m, t = get_model_and_toas(parfile, timfile, planets=True)
    m.DMXR2_0001.value = m.DMXR2_0002.value
    with pytest.raises(ValueError):
        SPNTA.from_pint(m, t)
