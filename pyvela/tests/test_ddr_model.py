from io import StringIO

import pytest
from pint.models import get_model
from pint.simulation import make_fake_toas_uniform

from pyvela.model import (
    _refuse_nondefault_ddr_galaxy,
    _refuse_nonzero_ddr_placeholders,
    center_model_epochs,
    fix_params,
    pint_components_to_vela,
)
from pyvela.spnta import convert_model_and_toas
from pyvela.vela import jl, vl


def _model_and_toas(binary_lines: str):
    model = get_model(
        StringIO(
            f"""
PSR TEST
EPHEM DE440
CLOCK TT(BIPM2023)
UNITS TDB
RAJ 12:00:00
DECJ 30:00:00
PMRA 0
PMDEC 0
PX 1
F0 100
PEPOCH 58000
PLANET_SHAPIRO N
{binary_lines}
"""
        )
    )
    toas = make_fake_toas_uniform(58000, 58001, 2, model, obs="@")
    toas.compute_posvels(ephem=model["EPHEM"].value, planets=True)
    fix_params(model, toas)
    return model, toas


def test_ell1_fbx_chart_detection_after_pint_setup():
    model, toas = _model_and_toas(
        """
BINARY ELL1
FB0 1.1574074074074073e-5
A1 5
TASC 58000
EPS1 0.02
EPS2 -0.03
"""
    )

    assert model["FB0"].quantity is not None
    assert model["PB"].quantity is not None
    components = pint_components_to_vela(model, toas)
    assert any(jl.isa(component, vl.BinaryELL1) for component in components)


def test_ddr_component_dispatch():
    model, toas = _model_and_toas(
        """
BINARY DDR
PB 1
A1 5
TASC 58000
EPS1 0.02
EPS2 -0.03
M2 0.8
COSI 0.5
DDRPK Y
DDRPBDOT absorb_gw
DDRGEO N
DDRKINE N
"""
    )

    tasc = model["TASC"].value
    center_model_epochs(model, toas)
    assert model["TASC"].value == tasc

    components = pint_components_to_vela(model, toas)
    assert any(jl.isa(component, vl.BinaryDDR) for component in components)
    vela_model, vela_toas = convert_model_and_toas(model, toas, [], [], {})
    assert len(vela_toas) == len(toas)
    assert any(jl.isa(component, vl.BinaryDDR) for component in vela_model.components)


def test_ddr_rejects_unsupported_parameter_values():
    model, _ = _model_and_toas(
        """
BINARY DDR
PB 1
A1 5
TASC 58000
EPS1 0.02
EPS2 -0.03
M2 0.8
COSI 0.5
DDRPK Y
DDRPBDOT absorb_gw
DDRGEO N
DDRKINE N
"""
    )

    model["DDRR0"].value = 9.0
    with pytest.raises(ValueError, match="DDRR0"):
        _refuse_nondefault_ddr_galaxy(model)

    model["DDRR0"].value = 8.178
    model["EPS1DOT"].value = 1e-20
    with pytest.raises(ValueError, match="EPS1DOT"):
        _refuse_nonzero_ddr_placeholders(model)
