from functools import cached_property

from .spnta import SPNTA
from .vela import Vela as vl, jl


class SPNA:
    def __init__(self, spnta: SPNTA):
        self.spnta = spnta

    @cached_property
    def residuals(self):
        return self.spnta.time_residuals(self.spnta.default_params)

    @cached_property
    def white_noise_components(self):
        return vl.get_white_noise_components(self.spnta.model)
