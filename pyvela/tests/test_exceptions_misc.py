from pint.models.parameter import floatParameter
from pyvela.parameters import get_scale_factor

import pytest


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
