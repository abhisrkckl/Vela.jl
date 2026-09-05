from functools import cached_property
from typing import List, Optional
from copy import deepcopy

import numpy as np
import emcee

from pint.models import TimingModel
from pint.toa import TOAs

from .spnta import SPNTA


class SDTSampler:
    def __init__(
        self,
        spnta: SPNTA,
        data_tempering_factor: float = 0.75,
        ntoa_min: int = 32,
        nwalkers_per_param=5,
    ):
        self.spnta = spnta
        self.data_tempering_factor = data_tempering_factor
        self.ntoa_min = ntoa_min
        self.ndim = self.spnta.ndim
        self.nwalkers = self.ndim * nwalkers_per_param

    @cached_property
    def spnta_subsets(self) -> List[SPNTA]:
        tzrtoa = self.spnta.model_pint.get_TZR_toa(self.spnta.toas_pint)
        tzrtoa.compute_pulse_numbers(self.spnta.model_pint)

        model_ = deepcopy(self.spnta.model_pint)

        spntas = []
        spnta1 = self.spnta
        while True:
            print(len(spnta1.toas_pint))
            spntas.append(spnta1)
            toas1 = get_toas_subset(
                self.spnta.model_pint_modified,
                spnta1.toas_pint,
                self.data_tempering_factor,
                self.ntoa_min,
            )

            if toas1 is None:
                break

            for par in ["TNREDC", "TNDMC", "TNCHROMC"]:
                if par in model_:
                    model_[par].value = max(
                        int(round(model_[par].value * self.data_tempering_factor)), 4
                    )

            spnta1 = SPNTA.from_pint(
                model_,
                toas1,
                analytic_marginalized_params=self.spnta.analytic_marginalized_params,
                custom_priors=(
                    self.spnta.custom_priors_dict
                    if hasattr(self.spnta, "custom_priors_dict")
                    else {}
                ),
                tzrtoa=tzrtoa,
            )
        spntas.reverse()

        return spntas

    @cached_property
    def samplers(self):
        return [
            emcee.EnsembleSampler(
                self.nwalkers,
                self.ndim,
                spnta_.lnpost_vectorized,
                vectorize=True,
                moves=[(emcee.moves.StretchMove(), 0.5), (emcee.moves.DEMove(), 0.5)],
            )
            for spnta_ in self.spnta_subsets
        ]

    @cached_property
    def final_sampler(self):
        return self.samplers[-1]

    def run_mcmc(self, x0, nsteps_initial=6000, nsteps_mid=1000, nsteps_final=6000):
        for ii, sampler in enumerate(self.samplers):
            if ii == 0:
                nsteps = nsteps_initial
            elif ii == len(self.spnta_subsets) - 1:
                nsteps = nsteps_final
            else:
                nsteps = nsteps_mid

            sampler.run_mcmc(x0, nsteps, progress=True)

            if ii < len(self.samplers) - 1:
                ch = sampler.get_chain(discard=nsteps // 2, flat=True)
                idx = np.random.choice(len(ch), size=self.nwalkers, replace=False)

                # x0 = sampler.get_chain()[-1,:,:]
                x0 = ch[idx, :]

    def get_chains(self, thin=10):
        return [sampler.get_chain(thin=thin) for sampler in self.samplers]

    def get_full_chain(self, thin=10):
        chains = self.get_chains(thin=thin)
        return np.concatenate(chains, axis=0)

    def get_final_chain(self, flat=True, thin=10, discard=0):
        return self.final_sampler.get_chain(flat=flat, thin=thin, discard=discard)


def get_toas_subset(
    model: TimingModel, toas: TOAs, data_tempering_factor: float, ntoa_min: int
) -> Optional[SPNTA]:
    ntoas = len(toas)

    selected_toa_idxs = []
    if "ScaleToaError" in model.components:
        for efac in model.EFACs:
            mask = model[efac].select_toa_mask(toas)
            nmask = len(mask)
            ntoas_sel = max(int(nmask * data_tempering_factor), min(ntoa_min, nmask))
            mask_selected = np.random.choice(mask, ntoas_sel, replace=False)
            selected_toa_idxs.extend(mask_selected)
    else:
        ntoas_sel = max(int(ntoas * data_tempering_factor), min(ntoa_min, ntoas))
        selected_toa_idxs.extend(np.random.choice(ntoas, ntoas_sel, replace=False))
    selected_toa_idxs.sort()

    return toas[selected_toa_idxs] if len(selected_toa_idxs) < len(toas) else None
