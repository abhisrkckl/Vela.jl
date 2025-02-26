from functools import cached_property

import numpy as np
import astropy.units as u
from scipy.linalg import cholesky, solve_triangular, cho_solve

from pint.models.noise_model import NoiseComponent, PLRedNoise, PLDMNoise, PLChromNoise
from pint.residuals import Residuals

from .spnta import SPNTA


class NonCentralParametrization:
    def __init__(self, spnta: SPNTA):
        self.spnta = spnta

    @property
    def model(self):
        return self.spnta.model_pint

    @property
    def toas(self):
        assert (
            self.spnta.toas_pint is not None
        ), "`toas_pint` is None! This usually happens when an SPNTA object is read from a JLSO file. Try reading it from par & tim files directly."
        return self.spnta.toas_pint

    @cached_property
    def efac_masks(self):
        return (
            {ef: self.model[ef].select_toa_mask(self.toas) for ef in self.model.EFACs}
            if hasattr(self.model, "EFACs")
            else {}
        )

    @cached_property
    def equad_masks(self):
        return (
            {eq: self.model[eq].select_toa_mask(self.toas) for eq in self.model.EQUADs}
            if hasattr(self.model, "EQUADs")
            else {}
        )

    @cached_property
    def param_names(self):
        noise_params = self.model.get_params_of_component_type("NoiseComponent")
        return [p for p in self.model.free_params if p in noise_params]

    @cached_property
    def default_param_values(self):
        return np.array([self.model[p].value for p in self.param_names])

    @cached_property
    def toa_errors(self):
        return self.toas.get_errors()

    @cached_property
    def param_index_mapping(self):
        return {p: ii for ii, p in enumerate(self.param_names)}

    def get_Ndiag(self, params):
        sigma = self.toa_errors
        efac = np.ones(len(sigma)) * u.dimensionless_unscaled
        equad = np.zeros(len(sigma)) * u.us

        for ef, efmask in self.efac_masks.items():
            idx = self.param_index_mapping[ef]
            efac[efmask] = params[idx] * u.dimensionless_unscaled

        for eq, eqmask in self.equad_masks.items():
            idx = self.param_index_mapping[eq]
            equad[eqmask] = params[idx] * u.us

        return (sigma**2 + equad**2) * efac**2

    @cached_property
    def _designmatrix_params_units(self):
        return self.model.designmatrix(self.toas)

    @cached_property
    def designmatrix(self):
        return self._designmatrix_params_units[0]

    @cached_property
    def designmatrix_params(self):
        return self._designmatrix_params_units[1]

    @cached_property
    def designmatrix_units(self):
        return self._designmatrix_params_units[2]

    @cached_property
    def noise_designmatrix(self):
        return self.model.noise_model_designmatrix(self.toas)

    @cached_property
    def noise_designmatrix_params(self):
        param_names = []
        for nc in self.model.NoiseComponent_list:
            nc: NoiseComponent
            if isinstance(nc, PLRedNoise):
                for ii in range(1, int(nc.TNREDC.value) + 1):
                    param_names.append(f"PLREDSIN_{ii:04d}")
                    param_names.append(f"PLREDCOS_{ii:04d}")
            elif isinstance(nc, PLDMNoise):
                for ii in range(1, int(nc.TNDMC.value) + 1):
                    param_names.append(f"PLDMSIN_{ii:04d}")
                    param_names.append(f"PLDMCOS_{ii:04d}")
            elif isinstance(nc, PLChromNoise):
                for ii in range(1, int(nc.TNCHROMC.value) + 1):
                    param_names.append(f"PLCHROMSIN_{ii:04d}")
                    param_names.append(f"PLCHROMCOS_{ii:04d}")

        return param_names

    @cached_property
    def noise_designmatrix_units(self):
        return [u.dimensionless_unscaled] * len(self.noise_designmatrix_params)

    def get_full_designmatrix(self, params):
        M_tm = self.designmatrix
        M_nm = self.noise_designmatrix
        sigma = np.sqrt(self.get_noise_weights(params))
        return np.hstack((M_tm, M_nm * sigma)) if M_nm is not None else M_tm

    @cached_property
    def residuals(self):
        return Residuals(self.toas, self.model).calc_time_resids()

    @cached_property
    def weights(self):
        return np.ones(len(self.designmatrix_params)) * 1e40

    def get_noise_weights(self, params):
        weights = []
        for nc in self.model.NoiseComponent_list:
            nc: NoiseComponent
            if isinstance(nc, PLRedNoise):
                log10_A = params[self.param_index_mapping["TNREDAMP"]]
                gamma = params[self.param_index_mapping["TNREDGAM"]]
                nharms = int(self.model.TNREDC.value)
            elif isinstance(nc, PLDMNoise):
                log10_A = params[self.param_index_mapping["TNDMAMP"]]
                gamma = params[self.param_index_mapping["TNDMGAM"]]
                nharms = int(self.model.TNDMC.value)
            elif isinstance(nc, PLChromNoise):
                log10_A = params[self.param_index_mapping["TNCHROMAMP"]]
                gamma = params[self.param_index_mapping["TNCHROMGAM"]]
                nharms = int(self.model.TNCHROMC.value)

            A = 10**log10_A

            f1 = 1 / self.toas.get_Tspan()
            fyr = 1 / u.year
            freqs = np.arange(1, nharms + 1) * f1
            freqs = np.repeat(freqs, 2)

            weights_nc = (
                A * A / (12 * np.pi**2 * fyr**3) * f1 * (fyr / freqs) ** gamma
            ).to_value(u.s**2)

            weights.extend(weights_nc)

        return np.array(weights)

    @cached_property
    def full_weights(self):
        phidiag_tm = self.weights
        phidiag_nm = np.ones(len(self.noise_designmatrix_params))
        return np.append(phidiag_tm, phidiag_nm)

    @cached_property
    def default_timing_model_params(self):
        return np.array(
            [
                float(self.model[p].value) * int(p != "F0")
                for p in self.model.free_params
                if p not in self.param_names
            ]
        )

    @cached_property
    def default_all_params(self):
        return np.append(
            self.default_timing_model_params,
            np.zeros(len(self.noise_designmatrix_params)),
        )

    @cached_property
    def timing_model_param_names(self):
        return [p for p in self.model.free_params if p not in self.param_names]

    @cached_property
    def all_param_names(self):
        return np.array(
            self.timing_model_param_names
            + self.noise_designmatrix_params
            + self.param_names
        )

    def transform_timing_model_params(self, all_params):
        ndim_tm = len(self.timing_model_param_names) + len(
            self.noise_designmatrix_params
        )
        alpha = all_params[:ndim_tm]
        theta = all_params[ndim_tm:]

        Ndiag = self.get_Ndiag(theta).to_value(u.s**2)
        M = self.get_full_designmatrix(theta)

        Ninv_M = M / Ndiag[:, None]
        MT_Ninv_M = M.T @ Ninv_M
        Phidiag = self.full_weights
        Sigmainv = np.diag(1 / Phidiag) + MT_Ninv_M
        Sigmainv_cf = cholesky(Sigmainv, lower=False)  # Sigmainv = L^T L

        y = self.residuals.to_value(u.s)
        MT_Ninv_y = Ninv_M.T @ y
        ahat = cho_solve((Sigmainv_cf, False), MT_Ninv_y)

        da = solve_triangular(Sigmainv_cf, alpha, lower=False)
        a0 = self.default_all_params
        a = a0 - ahat + da

        return np.append(a, theta)

    @cached_property
    def param_reorder_mask(self):
        return np.array(
            [
                np.where(self.all_param_names == p)[0].item()
                for p in self.spnta.param_names
            ]
        )

    def lnpost_ncp(self, params):
        params1 = (
            self.transform_timing_model_params(params)[self.param_reorder_mask]
            * self.spnta.scale_factors
        )
        return self.spnta.lnpost(params1)

    def get_ncp_sample(self):
        alpha = np.random.randn(
            len(self.timing_model_param_names) + len(self.noise_designmatrix_params)
        )

        params = (
            self.spnta.prior_transform(np.random.rand(self.spnta.ndim))
            / self.spnta.scale_factors
        )
        theta = [
            params[list(self.spnta.param_names).index(p)] for p in self.param_names
        ]

        return np.append(alpha, theta)
