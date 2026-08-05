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

If you wish to use the `pyvela-poco` script, also run
```
(vela) $ pip install pocomc
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
(vela) $ python -c 'from pyvela import SPNTA, __version__; print(__version__)'
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

## Installing in MacOS

The above instructions should work also in MacOS, except that Mac machines with Apple M-series chips don't support the `long double` type natively. This is not a problem for `Vela.jl` itself, since it uses the `DoubleFloats` package to handle extended precision arithmetic. However, `PINT` relies on the `numpy.longdouble` type, and will not work on Apple M-series machines normally.

Furthermore, because `PyVela` bridges Python and Julia using `juliacall`, any architecture mismatch between a native `arm64` Python environment and Intel-based dependencies will cause dynamic linker crashes.

To resolve this, both Python and Julia must be forced into an isolated Intel emulation environment (`osx-64`) using Rosetta 2 and Conda. 

For a complete, step-by-step walkthrough on how to configure this environment, handle the `numpy` dependencies, and properly initialize the Julia backend, please refer to the comprehensive guide:
📄 **[View the Complete Apple Silicon Installation Guide](./Vela_MacOS.pdf)**

### Troubleshooting: PyVela Version Error in Jupyter
Because `PyVela` is cloned manually rather than installed via a standard `pip` wheel, it lacks standard package metadata. This can cause version-related errors when executing code.

To easily bypass this without modifying the underlying source files, you can manually assign a dummy version string right after importing it in your Jupyter Notebook:

```python
import pyvela
pyvela.__version__ = "1.0.0"


## Using with `apptainer`

An `apptainer` definition file for `Vela.jl` along with its dependencies is available 
[here](https://github.com/abhisrkckl/Vela.jl/blob/main/apptainer/vela.def).