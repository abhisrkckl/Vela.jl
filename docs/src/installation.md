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

Furthermore, because `pyvela` bridges Python and Julia using `juliacall`, any architecture mismatch between a native `arm64` Python environment and Intel-based dependencies (like `x86_64`) will cause the dynamic linker to crash.

To resolve this, both Python and Julia must be forced into an isolated Intel emulation environment using Rosetta 2 and Conda. Below are the complete step-by-step instructions to configure this environment.

### Creating an Intel-Based Conda Environment
To ensure all compiled binaries are compatible, we must override Conda's default behavior and force it to pull Intel-based (`osx-64`) packages.

**Terminal:**
```bash
# Force Conda to use the Intel architecture for this environment
CONDA_SUBDIR=osx-64 conda create -n vela python=3.12 -c conda-forge -y

# Activate the environment
conda activate vela

# Lock the environment to osx-64 permanently
conda config --env --set subdir osx-64

# Verify the environment is correctly emulating the Intel architecture
python -c "import platform; print(platform.machine())"  # Should print "x86_64"
echo "CONDA_SUBDIR: $CONDA_SUBDIR"  # Should print "CONDA_SUBDIR: osx-64"
```

### Installing PINT and an Isolated Intel Julia
By default, `juliapkg` will search the macOS system for a global Julia installation (which is typically `arm64`). To prevent crashes, we must install an isolated, Intel version of Julia directly into the Conda environment alongside our Python tools.

**Terminal:**
```bash
# Install Julia within the Conda environment
conda install -c conda-forge julia -y

# Install PINT 
python -m pip install pint-pulsar
```

To ensure `juliapkg` does not use old, cached global environments, clear its hidden cache directories:

**Terminal:**
```bash
rm -rf ~/miniconda3/envs/vela/julia_env
rm -rf ~/.juliapkg
```

### Downloading and Hot-Patching pyvela
Because `pyvela` is cloned directly from its repository rather than installed via a formal `pip` wheel on PyPI, we must download and install it manually.

First, clone the repository and install it into your active Conda environment:

**Terminal:**
```bash
# Clone the pyvela repository
git clone https://github.com/abhisrkckl/pyvela.git

# Navigate into the newly downloaded directory
cd pyvela

# Install the package into the conda environment
python -m pip install .
```

Because of this manual installation method, it lacks standard package metadata. This causes crashes when `spnta` checks for a version string. We must bypass the metadata check in the `__init__.py` file and pull the version directly from Julia.

**Terminal:**
```bash
# Open the initialization file
nano ~/miniconda3/envs/vela/lib/python3.12/site-packages/pyvela/pyvela/__init__.py
```

Change the version definition line to the following:

**Text Editor (Inside Terminal):**
```python
__version__ = Vela.pkg_version()
```
Save the file and exit (`Ctrl+O`, `Enter`, then `Ctrl+X` if using nano).

### Troubleshooting: Version Mismatch Errors
Because `pyvela` acts as a direct interface to `Vela.jl`, the two packages are designed to be updated together. If their versions fall out of sync, you will receive a version mismatch error when trying to run `import pyvela`.

While the proper fix is to update both packages (see the **Updating Vela.jl** section), you might need a quick workaround if you are actively testing code in a Jupyter Notebook and cannot update immediately. You can temporarily resolve the error by manually assigning a dummy version string directly in your script:

```python
import pyvela
pyvela.__version__ = "1.0.0"
```

### Compiling the Julia Backend
With the Intel architectures aligned, we initialize the Julia backend from within Python to download `Vela.jl` and its required registered/unregistered packages.

Run the following inside a Jupyter Notebook or Python script:

**Python / Jupyter Notebook:**
```python
from juliacall import Main as jl

# Tell Julia to use its native package manager
jl.seval('import Pkg')

# Install standard dependencies
jl.seval('Pkg.add("Distributions")')
jl.seval('Pkg.add("DoubleFloats")')

# Install unregistered backend dependencies
jl.seval('Pkg.add(url="https://github.com/abhisrkckl/GeometricUnits.jl")')
jl.seval('Pkg.add(url="https://github.com/abhisrkckl/Vela.jl")')

# Optional: Update Vela.jl to the latest master branch
jl.seval('Pkg.update("Vela")')
```



## Using with `apptainer`

An `apptainer` definition file for `Vela.jl` along with its dependencies is available 
[here](https://github.com/abhisrkckl/Vela.jl/blob/main/apptainer/vela.def).