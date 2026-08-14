from functools import cached_property
from typing import List, Optional

import numpy as np
from joblib import Parallel, delayed

from pint.models import TimingModel
from pint.toa import TOAs

from .spnta import SPNTA

class SDTSampler:
    def __init__(self, spnta: SPNTA, data_tempering_factor: float = 0.75, ntoa_min: int = 64):
        self.spnta = spnta
        self.data_tempering_factor = data_tempering_factor
        self.ntoa_min = ntoa_min

    @cached_property
    def spnta_subsets(self) -> List[SPNTA]:
        toass = []
        toas1 = self.spnta.toas_pint
        while toas1 is not None:
            print("ntoas = ", len(toas1))
            toass.append(toas1)
            toas1 = get_toas_subset(self.spnta.model_pint_modified, toas1, self.data_tempering_factor, self.ntoa_min)
        toass.reverse()

        return Parallel(n_jobs=min(8, len(toass)))(
            delayed(SPNTA.from_pint)(
                self.spnta.model_pint_modified,
                toas,
                analytic_marginalized_params=self.spnta.analytic_marginalized_params,
                custom_priors=self.spnta.custom_priors_dict if hasattr(self.spnta, "custom_priors_dict") else {}
            )
            for toas in toass
        )


def get_toas_subset(model: TimingModel, toas: TOAs, data_tempering_factor: float, ntoa_min: int) -> Optional[SPNTA]:
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
        selected_toa_idxs.extend(
            np.random.choice(ntoas, ntoas_sel, replace=False)
        )
    selected_toa_idxs.sort()

    return toas[selected_toa_idxs] if len(selected_toa_idxs) < len(toas) else None