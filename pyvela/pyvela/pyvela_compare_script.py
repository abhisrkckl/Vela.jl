from argparse import ArgumentParser

from astropy import units as u
from matplotlib import pyplot as plt
from pint import DMconst, dmu
from pint.logging import setup as setup_log
from pint.models import get_model_and_toas
from pint.residuals import Residuals, WidebandTOAResiduals
from pyvela import SPNTA

setup_log(level="WARNING")


def parse_args(argv):
    parser = ArgumentParser(
        prog="pyvela-compare",
        description="Compare residuals obtained using Vela.jl and PINT.",
    )
    parser.add_argument("par_file")
    parser.add_argument("tim_file")
    parser.add_argument("-P", "--prior_file")

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    m, t = get_model_and_toas(args.par_file, args.tim_file, planets=True)
    t.compute_pulse_numbers(m)

    spnta = SPNTA.from_pint(
        m,
        t,
        custom_priors=(args.prior_file if args.prior_file is not None else {}),
    )

    if not t.is_wideband():
        res = Residuals(t, m)

        rv = spnta.time_residuals(spnta.maxlike_params)
        rp = res.calc_time_resids(subtract_mean=False).si.value
        terr = res.get_data_error().si.value

        plt.errorbar(
            t.get_mjds(),
            rp,
            terr,
            ls="",
            marker="+",
            color="red",
            label="PINT",
        )
        plt.errorbar(
            t.get_mjds(),
            rv,
            terr,
            ls="",
            marker="+",
            color="blue",
            label="Vela",
        )
        plt.ylabel("Time residuals (s)")
        plt.xlabel("MJD")
        plt.legend()
        plt.show()
    else:
        res = WidebandTOAResiduals(t, m)

        rv = spnta.time_residuals(spnta.maxlike_params)
        rp = res.toa.calc_time_resids(subtract_mean=False).si.value
        terr = res.toa.get_data_error().si.value

        dv = (spnta.dm_residuals(spnta.maxlike_params) * u.Hz / DMconst).to_value(dmu)
        dp = res.dm.calc_resids().to_value(dmu)
        derr = res.dm.get_data_error().to_value(dmu)

        plt.subplot(211)
        plt.errorbar(
            t.get_mjds(),
            rp,
            terr,
            ls="",
            marker="+",
            color="red",
            label="PINT",
        )
        plt.errorbar(
            t.get_mjds(),
            rv,
            terr,
            ls="",
            marker="+",
            color="blue",
            label="Vela",
        )
        plt.ylabel("Time residuals (s)")
        plt.legend()

        plt.subplot(212)
        plt.errorbar(
            t.get_mjds(),
            dp,
            derr,
            ls="",
            marker="+",
            color="red",
            label="PINT",
        )
        plt.errorbar(
            t.get_mjds(),
            dv,
            derr,
            ls="",
            marker="+",
            color="blue",
            label="Vela",
        )
        plt.ylabel("DM residuals (dmu)")
        plt.xlabel("MJD")
        plt.show()
