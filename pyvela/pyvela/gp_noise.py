from astropy import units as u
from astropy.time import Time
from pint.models.noise_model import PLChromNoise, PLDMNoise, PLRedNoise
from pint.models.parameter import MJDParameter, prefixParameter
from pint.models.timing_model import DelayComponent


class PLRedNoiseGP(DelayComponent):
    """A dummy PINT component that is used to translate `PLRedNoise` to a form
    with unmarginalized amplitudes."""

    def __init__(self, plrednoise: PLRedNoise, epoch: Time):
        super().__init__()

        self.add_param(plrednoise.TNREDAMP)
        self.add_param(plrednoise.TNREDGAM)
        self.add_param(plrednoise.TNREDC)
        self.add_param(plrednoise.PLREDFREQ)
        self.add_param(plrednoise.TNREDFLOG)
        self.add_param(plrednoise.TNREDFLOG_FACTOR)

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

        Nharm = int(plrednoise.TNREDC.value) + (
            int(plrednoise.TNREDFLOG.value)
            if plrednoise.TNREDFLOG.value is not None
            else 0
        )
        for ii in range(1, Nharm + 1):
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
    with unmarginalized amplitudes."""

    def __init__(self, pldmnoise: PLDMNoise, epoch: Time):
        super().__init__()

        self.add_param(pldmnoise.TNDMAMP)
        self.add_param(pldmnoise.TNDMGAM)
        self.add_param(pldmnoise.TNDMC)
        self.add_param(pldmnoise.PLDMFREQ)
        self.add_param(pldmnoise.TNDMFLOG)
        self.add_param(pldmnoise.TNDMFLOG_FACTOR)

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

        Nharm = int(pldmnoise.TNDMC.value) + (
            int(pldmnoise.TNDMFLOG.value) if pldmnoise.TNDMFLOG.value is not None else 0
        )
        for ii in range(1, Nharm + 1):
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
    with unmarginalized amplitudes."""

    def __init__(self, plchromnoise: PLChromNoise, epoch: Time):
        super().__init__()

        self.add_param(plchromnoise.TNCHROMAMP)
        self.add_param(plchromnoise.TNCHROMGAM)
        self.add_param(plchromnoise.TNCHROMC)
        self.add_param(plchromnoise.PLCHROMFREQ)
        self.add_param(plchromnoise.TNCHROMFLOG)
        self.add_param(plchromnoise.TNCHROMFLOG_FACTOR)

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

        Nharm = int(plchromnoise.TNCHROMC.value) + (
            int(plchromnoise.TNCHROMFLOG.value)
            if plchromnoise.TNCHROMFLOG.value is not None
            else 0
        )
        for ii in range(1, Nharm + 1):
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
