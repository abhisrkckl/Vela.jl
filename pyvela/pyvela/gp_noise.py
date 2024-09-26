from pint.models.timing_model import DelayComponent
from pint.models.noise_model import PLRedNoise
from pint.models.parameter import floatParameter, prefixParameter, MJDParameter

from astropy import units as u
from astropy.time import Time

class PLRedNoiseGP(DelayComponent):
    """A dummy PINT component that is used to translate `PLRedNoise` to a form
    `Vela.jl` can handle."""

    def __init__(self, plrednoise: PLRedNoise, f1: u.Quantity, epoch: Time):
        self.add_param(plrednoise.TNREDAMP)
        self.add_param(plrednoise.TNREDGAM)

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

        for ii in range(1, self.TNREDC.value + 1):
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
        raise NotImplementedError("This is a dummy component for dealing with `Vela.jl` only and cannot be evaluated.")
