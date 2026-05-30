# Installation

`Vela.jl` can be installed directly from GitHub. We recommend installing it within a dedicated 
`conda` environment. See the instructions [here](https://www.anaconda.com/docs/getting-started/miniconda/install/linux-install)
to install `miniconda`.

Please note that `Vela.jl` is only tested against Python 3.12 and Julia 1.11 in Ubuntu at present.

If you don't have Julia installed, please install it using `juliaup` following the instructions found 
[here](https://julialang.org/downloads/).

!!! warning
    Some of the dependencies don't work properly if Julia isn't installed using `juliaup`.
    Specifically, avoid installing Julia using `conda`. The following instructions are for
    installing the Python dependencies *only* in the `conda` environment.

The following instructions assume that both `conda` and `julia` are installed.

Now, install the Python dependencies and set the environment variables. The most important one is
`PYTHON_JULIACALL_HANDLE_SIGNALS`. If it is not set properly you'll get segmentation faults.
```
(base) $ # Setup conda environment
(base) $ conda create -n vela python=3.12
(base) $ conda activate vela
(vela) $ conda env config vars set PYTHON_JULIACALL_HANDLE_SIGNALS=yes
(vela) $ conda env config vars set PYTHON_JULIACALL_THREADS=4
(vela) $ conda env config vars set JULIA_NUM_THREADS=4
(vela) $ conda env config vars set JULIA_CONDAPKG_BACKEND="Null"
(vela) $ conda env config vars set PYTHON_JULIACALL_EXE="$(which julia)"
(vela) $ conda env config vars set PYTHON_JULIACALL_PROJECT=$(julia -e 'print(joinpath(DEPOT_PATH[1], "environments", "v$(VERSION.major).$(VERSION.minor)"))')
(vela) $ conda env config vars set PYTHON_JULIAPKG_OFFLINE=true
(vela) $ conda install -c conda-forge pyjuliacall black emcee nestle corner tqdm pytest pytest-xdist
(vela) $ pip install git+https://github.com/nanograv/PINT
```

The number of threads available to `Vela.jl` for parallel processing can be controlled 
using the environment variables `JULIA_NUM_THREADS` (for direct use from Julia) or 
`PYTHON_JULIACALL_THREADS` (for use from within Python). This should be set based on the 
number of CPU cores available in your machine.

Now, install the Julia packages.
```
(vela) $ julia
julia> import Pkg
julia> Pkg.Registry.add("General")
julia> Pkg.Registry.add(url="https://github.com/abhisrkckl/julia_registry")
julia> Pkg.add(["LocalRegistry", "JuliaFormatter", "BenchmarkTools", "PythonCall", "Distributions", "DoubleFloats", "GeometricUnits"])
julia> Pkg.add(url="https://github.com/abhisrkckl/Vela.jl")
julia> exit()
```

Install the Python interface `pyvela`.
```
(vela) $ pip install git+https://github.com/abhisrkckl/Vela.jl
```

The `pyvela/examples` directory provides several example datasets and scripts.
A basic example (using the Python wrapper) can be run like this:
```
(vela) $ ./run_example_emcee.py NGC6440E.par NGC6440E.tim
```

## Testing if installation is successful

Run the following commands to check if the installation is successful.

```
(vela) $ julia -e 'import Pkg; Pkg.test("Vela")'
(vela) $ python -c 'python -c 'from pyvela import SPNTA, __version__; print(__version__)''
```

## Updating `Vela.jl`

To update a `Vela.jl` installation, do the following.
```
(vela) $ julia
julia> import Pkg
julia> Pkg.update(["GeometricUnits", "Vela"])
julia> exit()
```

It is best to reinstall `pyvela`:
```
(vela) $ pip install git+https://github.com/abhisrkckl/Vela.jl
```

Note that both `Vela.jl` and `pyvela` should be updated together. You will get an
exception if a version mismatch is detected.