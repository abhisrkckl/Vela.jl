import json
import os
import shutil
from typing import Optional
import numpy as np
from .spnta import SPNTA


class SPNTAResults:
    def __init__(
        self,
        spnta: SPNTA,
        chain: np.ndarray,
        info_dict: dict,
        truth_parfile: Optional[str] = None,
    ):
        self.parfile = spnta.parfile
        self.timfile = spnta.timfile

        self.info_dict = info_dict

        self.custom_priors_dict = spnta.custom_priors_dict

        self.truth_parfile = truth_parfile

        self.param_names = spnta.param_names
        self.param_prefixes = spnta.param_prefixes
        self.param_units = spnta.param_units
        self.scale_factors = spnta.scale_factors

        self.default_values = spnta.default_params

        self.samples_raw = chain
        self.median_sample = np.median(chain, axis=0)

    @property
    def samples(self):
        return self.samples_raw / self.scale_factors

    def save_to_dir(self, outdir: str):
        pass

    def prepare_outdir(self, outdir: str):
        assert not os.path.isdir(outdir), "The output directory already exists!"

        os.mkdir(outdir)
        shutil.copy(self.parfile, outdir)
        shutil.copy(self.timfile, outdir)

        with open(f"{outdir}/custom_priors.json", "w") as prior_file:
            json.dump(self.custom_priors_dict, prior_file, indent=4)

        with open(f"{outdir}/summary.json", "w") as summary_file:
            json.dump(self.info_dict, summary_file, indent=4)

        if self.truth_parfile is not None:
            shutil.copy(self.truth_parfile, outdir)
