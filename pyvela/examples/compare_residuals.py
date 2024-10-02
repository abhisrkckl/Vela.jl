#!/usr/bin/env python

import sys

from matplotlib import pyplot as plt
from pint.logging import setup as setup_log
from pint.models import get_model_and_toas
from pint.residuals import Residuals

from pyvela import Vela as vl
from pyvela import read_model_and_toas

setup_log(level="WARNING")

parfile, timfile = sys.argv[1], sys.argv[2]

m, t = get_model_and_toas(parfile, timfile)
t.compute_pulse_numbers(m)

res = Residuals(t, m)

mv, tv = read_model_and_toas(parfile, timfile)

if not t.is_wideband():
    rv = list(
        map(vl.value, vl.form_residuals(mv, tv, mv.param_handler._default_params_tuple))
    )
else:
    rv = [
        vl.value(wr[0])
        for wr in vl.form_residuals(mv, tv, mv.param_handler._default_params_tuple)
    ]

plt.errorbar(
    t.get_mjds(),
    res.calc_time_resids().si.value,
    res.get_data_error().si.value,
    ls="",
    marker="+",
    color="red",
    label="PINT",
)
plt.errorbar(
    t.get_mjds(),
    rv,
    res.get_data_error().si.value,
    ls="",
    marker="+",
    color="blue",
    label="Vela",
)
plt.legend()
plt.show()
