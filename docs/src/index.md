# Vela.jl

A Bayesian pulsar timing and noise analysis package.

The `Vela.jl` project aims to develop a fast and simple-to-use framework for doing Bayesian
pulsar timing & noise analysis. It currently supports both narrowband and wideband TOAs along 
with most commonly used pulsar timing models. It can be used from Julia REPL, scripts, or 
notebooks. It can also be used from Python REPL, scripts, and notebooks with the help of the 
`pint2vela` interface.

It is under active development.

## Getting started

`Vela.jl` can be installed directly from GitHub. We recommend installing it within a dedicated 
`conda` environment.

The following commands install the Python and Julia packages that are needed for developing Vela.jl and running the examples and tests.

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
A basic example can be run like this:
```
(vela) $ ./run_example_emcee.py NGC6440E.par NGC6440E.tim
```
It produces a `corner` plot like this:
![NGC6440E_posterior](NGC6440E_posterior.png)

## Types
```@autodocs
Modules = [Vela]
Order = [:type]
```

## Functions
```@autodocs
Modules = [Vela]
Order = [:function]
```