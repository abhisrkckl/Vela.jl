# `VelaFitter`

`VelaFitter` is an implementation of the `PINT` `Fitter` interface using `Vela.jl` as
a backend. This is a thin wrapper over the `pyvela.SPNTA` class. It can be used like 
any other `PINT` `Fitter` for both narrowband and wideband TOAs. The fitting is done 
by sampling the posterior distribution using `emcee`. 

Unlike other `Fitter`s, the `toas` object herein is immutable, and changing its contents 
will not make any difference in the fitting.

An example is given below.

```
from pint.models import get_model_and_toa
from pyvela.fitter import VelaFitter

model, toas = get_model_and_toas("NGC6440E.par", "NGC6440E.tim")

ftr = VelaFitter(
    toas, 
    model, 
    custom_priors={
        "EQUAD": {
            "distribution": "LogUniform", 
            "args": [1e-3, 2.0]
        }
    }
)

ftr.fit_toas()

print(ftr.model)
```