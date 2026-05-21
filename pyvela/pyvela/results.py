from functools import cached_property
import json
import os

from pint.models import TimingModel, get_model
from pint.toa import TOAs, get_TOAs

import numpy as np

from .spnta import SPNTA


class SPNTAResults:
    def __init__(self, result_dir: str):
        self.result_dir = result_dir

    @cached_property
    def summary(self) -> dict:
        with open(f"{self.result_dir}/summary.json", "r") as f:
            return json.load(f)

    @cached_property
    def model_input(self) -> TimingModel:
        filename = self.summary["input"]["par_file"]
        return get_model(f"{self.result_dir}/{filename}", allow_T2=True, allow_tcb=True)

    @cached_property
    def toas_input(self) -> TOAs:
        filename = self.summary["input"]["tim_file"]
        return get_TOAs(
            f"{self.result_dir}/{filename}", planets=True, model=self.model_input
        )

    @cached_property
    def psrname(self) -> str:
        return self.model_input["PSR"].value

    @cached_property
    def epoch(self) -> float:
        return np.genfromtxt(f"{self.result_dir}/epoch.txt")

    @cached_property
    def model_median(self) -> TimingModel:
        filename = f"{self.psrname}_median.par"
        return get_model(f"{self.result_dir}/{filename}", allow_T2=True, allow_tcb=True)

    @cached_property
    def samples(self) -> np.ndarray:
        return np.load(f"{self.result_dir}/samples.npy")

    @cached_property
    def samples_raw(self) -> np.ndarray:
        return np.load(f"{self.result_dir}/samples_raw.npy")

    @cached_property
    def _residuals(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/residuals.txt")

    @cached_property
    def mjds(self) -> np.ndarray:
        return self._residuals[:, 0]

    @cached_property
    def time_residuals(self) -> np.ndarray:
        return self._residuals[:, 1]

    @cached_property
    def whitened_time_residuals(self) -> np.ndarray:
        return self._residuals[:, 2]

    @cached_property
    def scaled_toa_uncertainties(self) -> np.ndarray:
        return self._residuals[:, 3]

    @cached_property
    def is_wideband(self) -> np.ndarray:
        return self._residuals.shape[1] == 7

    @cached_property
    def dm_residuals(self) -> np.ndarray:
        return self._residuals[:, 4] if self.is_wideband else None

    @cached_property
    def whitened_dm_residuals(self) -> np.ndarray:
        return self._residuals[:, 5] if self.is_wideband else None

    @cached_property
    def scaled_dm_uncertainties(self) -> np.ndarray:
        return self._residuals[:, 6] if self.is_wideband else None

    @cached_property
    def prior_info(self) -> dict:
        with open(f"{self.result_dir}/prior_info.json", "r") as f:
            return json.load(f)

    @cached_property
    def prior_evals(self) -> np.ndarray:
        return np.load(f"{self.result_dir}/prior_evals.npy")

    @cached_property
    def param_stds(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/params_std.txt")

    @cached_property
    def param_medians(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/params_median.txt")

    @cached_property
    def param_units(self) -> np.ndarray:
        return np.genfromtxt(
            f"{self.result_dir}/param_units.txt", dtype=str, delimiter="~"
        )

    @cached_property
    def param_scale_factors(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/param_scale_factors.txt")

    @cached_property
    def param_prefixes(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/param_prefixes.txt", dtype=str)

    @cached_property
    def param_names(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/param_names.txt", dtype=str)

    @cached_property
    def param_maxpost_values(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/param_maxpost_values.txt")

    @cached_property
    def param_default_values(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/param_default_values.txt")

    @cached_property
    def param_true_values(self) -> np.ndarray:
        if os.path.isfile(f"{self.result_dir}/param_true_values.txt"):
            return np.genfromtxt(f"{self.result_dir}/param_true_values.txt")
        else:
            return None

    @cached_property
    def param_autocorr(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/param_autocorr.txt")

    @cached_property
    def marginalized_param_default_values(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/marginalized_param_default_values.txt")

    @cached_property
    def marginalized_param_names(self) -> np.ndarray:
        return np.genfromtxt(
            f"{self.result_dir}/marginalized_param_names.txt", dtype=str
        )

    @cached_property
    def marginalized_param_medians(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/marginalized_params_median.txt")

    @cached_property
    def marginalized_param_maxpost_values(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/marginalized_param_maxpost_values.txt")

    @cached_property
    def marginalized_param_scale_factors(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/marginalized_param_scale_factors.txt")

    @cached_property
    def marginalized_param_stds(self) -> np.ndarray:
        return np.genfromtxt(f"{self.result_dir}/marginalized_params_std.txt")
