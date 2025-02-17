# Installation

`Vela.jl` can be installed directly from GitHub. We recommend installing it within a dedicated 
`conda` environment.

Please note that `Vela.jl` is only tested against Python 3.12 and Julia 1.10 in Ubuntu at present.

If you don't have Julia installed, please install it following the instructions found 
[here](https://julialang.org/downloads/).

Now install the Python dependencies and set the environment variables. The most important one is
`PYTHON_JULIACALL_HANDLE_SIGNALS`. If it is not set properly you'll get segmentation faults.
```
(base) $ # Setup conda environment
(base) $ conda create -n vela python=3.12
(base) $ conda activate vela
(vela) $ conda install -c conda-forge pyjuliacall pint-pulsar>=1.1 black emcee nestle corner tqdm pytest pytest-xdist
(vela) $ conda env config vars set PYTHON_JULIACALL_HANDLE_SIGNALS=yes
(vela) $ conda env config vars set PYTHON_JULIACALL_THREADS=4
(vela) $ conda env config vars set JULIA_NUM_THREADS=4
(vela) $ conda env config vars set JULIA_CONDAPKG_BACKEND="Null"
```

Now install the Julia packages.
```
(vela) $ julia
julia> import Pkg
julia> Pkg.Registry.add(url="https://github.com/abhisrkckl/julia_registry")
julia> Pkg.add(["LocalRegistry", "JuliaFormatter", "BenchmarkTools", "PythonCall", "Distributions", "DoubleFloats", "GeometricUnits"])
julia> Pkg.add(url="https://github.com/abhisrkckl/Vela.jl")
julia> exit()
```

Install the Python interface `pyvela`.
```
(vela) $ pip install git+https://github.com/abhisrkckl/Vela.jl
```

The number of threads available to `Vela.jl` for parallel processing can be controlled 
using the environment variables `JULIA_NUM_THREADS` (for direct use from Julia) or 
`PYTHON_JULIACALL_THREADS` (for use from within Python).

The `pyvela/examples` directory provides several example datasets and scripts.
A basic example (using the Python wrapper) can be run like this:
```
(vela) $ ./run_example_emcee.py NGC6440E.par NGC6440E.tim
```
