import os

from pyvela.pint2vela import par_tim_to_jlso
from pyvela.vela import vl


def test_convert_jlso():
    par_tim_to_jlso(
        "datafiles/NGC6440E.par", "datafiles/NGC6440E.tim", "__NGC6440E.jlso"
    )
    assert os.path.isfile("__NGC6440E.jlso")

    mv, tv = vl.load_pulsar_data("__NGC6440E.jlso")
    assert len(tv) == 62
    assert set(vl.get_free_param_names(mv.param_handler)) == {
        "RAJ",
        "DECJ",
        "DM",
        "PHOFF",
        "F0",
        "F1",
    }

    os.unlink("__NGC6440E.jlso")
