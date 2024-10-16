# Installation

`Vela.jl` can be installed directly from GitHub. We recommend installing it within a dedicated 
`conda` environment.

Please note that `Vela.jl` is only tested against Python 3.12 and Julia 1.10 in Ubuntu at present.

The following commands install the Python and Julia packages that are needed for developing `Vela.jl` and running the examples and tests.

```
(base) $ # Setup conda environment
(base) $ conda create -n vela python=3.12
(base) $ conda activate vela
(vela) $ conda install -c conda-forge julia pyjuliacall pint-pulsar black emcee nestle corner tqdm pytest pytest-xdist
(vela) $ conda env config vars set PYTHON_JULIACALL_HANDLE_SIGNALS=yes
(vela) $ conda env config vars set PYTHON_JULIACALL_THREADS=4
(vela) $ conda env config vars set JULIA_NUM_THREADS=4
(vela) $ conda env config vars set JULIA_CONDAPKG_BACKEND="Null"
(vela) $ 
(vela) $ # Install Julia packages
(vela) $ julia
julia> ] add LocalRegistry, JuliaFormatter, BenchmarkTools
julia> ] registry add https://github.com/abhisrkckl/julia_registry
julia> ] add https://github.com/abhisrkckl/Vela.jl
julia> exit()
(vela) $ 
(vela) $ # Install Python wrapper
(vela) $ pip install git+https://github.com/abhisrkckl/Vela.jl
```

The number of threads available to `Vela.jl` for parallel processing can be controlled 
using the environment variables `JULIA_NUM_THREADS` (for direct use from Julia) or `PYTHON_JULIACALL_THREADS` (for use from within Python).

The `pint2vela/examples` directory provides several example datasets and scripts.
A basic example (using the Python wrapper) can be run like this:
```
(vela) $ ./run_example_emcee.py NGC6440E.par NGC6440E.tim
```
