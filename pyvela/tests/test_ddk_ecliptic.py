"""ICRS ↔ ecliptic DDK residual parity.

PINT's `as_ECL()` rotates coordinates and `KOM`. After the ecliptic DDK
annual-parallax frame fix, Vela residuals must agree across that conversion
at the picosecond-to-nanosecond level. A mixed-frame annual term disagrees
by ~2 μs on `sim_ddk`.
"""

from copy import deepcopy
import os

import numpy as np
from pint.models import get_model_and_toas
from pint.residuals import Residuals

from pyvela.spnta import SPNTA

datadir = f"{os.path.dirname(os.path.realpath(__file__))}/datafiles"


def _freeze(model):
    model = deepcopy(model)
    for p in list(model.free_params):
        model[p].frozen = True
    return model


def _max_abs_mean_sub(a, b):
    d = np.asarray(a) - np.asarray(b)
    return float(np.max(np.abs(d - d.mean())))


def _pint_resids(model, toas):
    return Residuals(toas, model, subtract_mean=False).time_resids.to_value("s")


def _vela_resids(model, toas):
    spnta = SPNTA.from_pint(_freeze(model), toas)
    return spnta.time_residuals(spnta.default_params)


def test_ddk_icrs_ecliptic_parity():
    par = f"{datadir}/sim_ddk.par"
    tim = f"{datadir}/sim_ddk.tim"
    m, t = get_model_and_toas(par, tim, planets=True)
    m_ecl = m.as_ECL()

    pint_icrs = _pint_resids(m, t)
    pint_ecl = _pint_resids(m_ecl, t)
    assert _max_abs_mean_sub(pint_icrs, pint_ecl) < 1e-12

    vela_icrs = _vela_resids(m, t)
    vela_ecl = _vela_resids(m_ecl, t)
    assert _max_abs_mean_sub(vela_icrs, vela_ecl) < 1e-9
    assert _max_abs_mean_sub(vela_ecl, pint_ecl) < 1e-9
    assert _max_abs_mean_sub(vela_icrs, pint_icrs) < 1e-9
