from astropy import units as u
from astropy.time import Time
from pint.models.noise_model import PLChromNoise, PLDMNoise, PLRedNoise
from pint.models.parameter import MJDParameter, floatParameter, prefixParameter
from pint.models.timing_model import DelayComponent


class PLRedNoiseGP(DelayComponent):
    """A dummy PINT component that is used to translate `PLRedNoise` to a form
    `Vela.jl` can handle."""

    def __init__(self, plrednoise: PLRedNoise, f1: u.Quantity, epoch: Time):
        super().__init__()

        self.add_param(plrednoise.TNREDAMP)
        self.add_param(plrednoise.TNREDGAM)
        self.add_param(plrednoise.TNREDC)

        self.TNREDAMP.tcb2tdb_scale_factor = 1
        self.TNREDGAM.tcb2tdb_scale_factor = 1

        self.add_param(
            MJDParameter(
                name="PLREDEPOCH",
                description="Epoch of the Powerlaw Fourier GP representation of the achromatic red noise",
                value=epoch.mjd,
                time_scale="tdb",
                tcb2tdb_scale_factor=u.Quantity(1),
                frozen=True,
            )
        )

        self.add_param(
            floatParameter(
                name="PLREDFREQ",
                description="Fundamental frequency of the Powerlaw Fourier GP representation of the achromatic red noise",
                units="1/year",
                value=f1.to_value("1/year"),
                tcb2tdb_scale_factor=u.Quantity(1),
                frozen=True,
            )
        )

        for ii in range(1, int(plrednoise.TNREDC.value) + 1):
            self.add_param(
                prefixParameter(
                    parameter_type="float",
                    name=f"PLREDSIN_{ii:04d}",
                    description="Sine amplitude of the Powerlaw Fourier GP representation of the achromatic red noise",
                    units=u.dimensionless_unscaled,
                    value=0.0,
                    frozen=False,
                    tcb2tdb_scale_factor=u.Quantity(1),
                )
            )

            self.add_param(
                prefixParameter(
                    parameter_type="float",
                    name=f"PLREDCOS_{ii:04d}",
                    description="Cosine amplitude of the Powerlaw Fourier GP representation of the achromatic red noise",
                    units=u.dimensionless_unscaled,
                    value=0.0,
                    frozen=False,
                    tcb2tdb_scale_factor=u.Quantity(1),
                )
            )

        self.delay_funcs_component += [self.dummy_delay]

    def dummy_delay(self, toas, delays):
        return 0 * u.s


class PLDMNoiseGP(DelayComponent):
    """A dummy PINT component that is used to translate `PLDMNoise` to a form
    `Vela.jl` can handle."""

    def __init__(self, pldmnoise: PLDMNoise, f1: u.Quantity, epoch: Time):
        super().__init__()

        self.add_param(pldmnoise.TNDMAMP)
        self.add_param(pldmnoise.TNDMGAM)
        self.add_param(pldmnoise.TNDMC)

        self.TNDMAMP.tcb2tdb_scale_factor = 1
        self.TNDMGAM.tcb2tdb_scale_factor = 1

        self.add_param(
            MJDParameter(
                name="PLDMEPOCH",
                description="Epoch of the Powerlaw Fourier GP representation of the DM noise",
                value=epoch.mjd,
                time_scale="tdb",
                tcb2tdb_scale_factor=u.Quantity(1),
                frozen=True,
            )
        )

        self.add_param(
            floatParameter(
                name="PLDMFREQ",
                description="Fundamental frequency of the Powerlaw Fourier GP representation of the DM noise",
                units="1/year",
                value=f1.to_value("1/year"),
                tcb2tdb_scale_factor=u.Quantity(1),
                frozen=True,
            )
        )

        for ii in range(1, int(pldmnoise.TNDMC.value) + 1):
            self.add_param(
                prefixParameter(
                    parameter_type="float",
                    name=f"PLDMSIN_{ii:04d}",
                    description="Sine amplitude of the Powerlaw Fourier GP representation of the DM noise",
                    units=u.dimensionless_unscaled,
                    value=0.0,
                    frozen=False,
                    tcb2tdb_scale_factor=u.Quantity(1),
                )
            )

            self.add_param(
                prefixParameter(
                    parameter_type="float",
                    name=f"PLDMCOS_{ii:04d}",
                    description="Cosine amplitude of the Powerlaw Fourier GP representation of the DM noise",
                    units=u.dimensionless_unscaled,
                    value=0.0,
                    frozen=False,
                    tcb2tdb_scale_factor=u.Quantity(1),
                )
            )

        self.delay_funcs_component += [self.dummy_delay]

    def dummy_delay(self, toas, delays):
        return 0 * u.s


class PLChromNoiseGP(DelayComponent):
    """A dummy PINT component that is used to translate `PLDMNoise` to a form
    `Vela.jl` can handle."""

    def __init__(self, plchromnoise: PLChromNoise, f1: u.Quantity, epoch: Time):
        super().__init__()

        self.add_param(plchromnoise.TNCHROMAMP)
        self.add_param(plchromnoise.TNCHROMGAM)
        self.add_param(plchromnoise.TNCHROMC)

        self.TNCHROMAMP.tcb2tdb_scale_factor = 1
        self.TNCHROMGAM.tcb2tdb_scale_factor = 1

        self.add_param(
            MJDParameter(
                name="PLCHROMEPOCH",
                description="Epoch of the Powerlaw Fourier GP representation of the chromatic noise",
                value=epoch.mjd,
                time_scale="tdb",
                tcb2tdb_scale_factor=u.Quantity(1),
                frozen=True,
            )
        )

        self.add_param(
            floatParameter(
                name="PLCHROMFREQ",
                description="Fundamental frequency of the Powerlaw Fourier GP representation of the chromatic noise",
                units="1/year",
                value=f1.to_value("1/year"),
                tcb2tdb_scale_factor=u.Quantity(1),
                frozen=True,
            )
        )

        for ii in range(1, int(plchromnoise.TNCHROMC.value) + 1):
            self.add_param(
                prefixParameter(
                    parameter_type="float",
                    name=f"PLCHROMSIN_{ii:04d}",
                    description="Sine amplitude of the Powerlaw Fourier GP representation of the chromatic noise",
                    units=u.dimensionless_unscaled,
                    value=0.0,
                    frozen=False,
                    tcb2tdb_scale_factor=u.Quantity(1),
                )
            )

            self.add_param(
                prefixParameter(
                    parameter_type="float",
                    name=f"PLCHROMCOS_{ii:04d}",
                    description="Cosine amplitude of the Powerlaw Fourier GP representation of the chromatic noise",
                    units=u.dimensionless_unscaled,
                    value=0.0,
                    frozen=False,
                    tcb2tdb_scale_factor=u.Quantity(1),
                )
            )

        self.delay_funcs_component += [self.dummy_delay]

    def dummy_delay(self, toas, delays):
        return 0 * u.s
