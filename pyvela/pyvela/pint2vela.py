from pint.binaryconvert import convert_binary
from pint.logging import setup as setup_log
from pint.models import get_model_and_toas

from .convert_model import pint_model_to_vela
from .convert_toas import pint_toas_to_vela
from .ecorr import ecorr_sort
from .vela import vl


def read_model_and_toas(
    parfile: str, timfile: str, cheat_prior_scale=20, custom_prior_dicts={}
):
    """Read a pair of par & tim files and create a `Vela.TimingModel` object and a
    Julia `Vector` of `TOA`s."""

    setup_log(level="WARNING")
    mp, tp = get_model_and_toas(
        parfile,
        timfile,
        planets=True,
        allow_tcb=True,
        allow_T2=True,
        add_tzr_to_model=True,
    )

    if "BinaryBT" in mp.components:
        mp = convert_binary(mp, "DD")

    if "EcorrNoise" in mp.components:
        assert not tp.is_wideband(), "ECORR is not supported for wideband data."
        tp, ecorr_toa_ranges, ecorr_indices = ecorr_sort(mp, tp)
    else:
        ecorr_toa_ranges, ecorr_indices = None, None

    model = pint_model_to_vela(
        mp,
        tp,
        cheat_prior_scale,
        custom_prior_dicts,
        ecorr_toa_ranges=ecorr_toa_ranges,
        ecorr_indices=ecorr_indices,
    )
    toas = pint_toas_to_vela(tp, float(mp.PEPOCH.value))

    return model, toas


def par_tim_to_jlso(parfile: str, timfile: str, jlsofile: str):
    """Save a pair of `par` & `tim` files as a `JLSO` file."""
    mv, tv = read_model_and_toas(parfile, timfile)
    vl.save_pulsar_data(jlsofile, mv, tv)
