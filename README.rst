Vela
----

.. image:: https://github.com/abhisrkckl/Vela.jl/actions/workflows/Vela-CI.yml/badge.svg
   :target: https://github.com/abhisrkckl/Vela.jl/actions
   :alt: CI Status

.. image:: https://codecov.io/gh/abhisrkckl/Vela.jl/graph/badge.svg?token=Y6ES2WTYEV 
   :target: https://codecov.io/gh/abhisrkckl/Vela.jl
   :alt: Coverage

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: LICENCE
   :alt: License: MIT

Installation
------------
The following commands install the Python and Julia packages that are needed for developing
``Vela.jl`` and running the examples and tests. 

::

   (base) $ conda create -n vela python=3.12
   (base) $ conda activate vela
   (vela) $ conda install -c conda-forge julia pyjuliacall pint-pulsar black emcee nestle corner tqdm pytest pytest-xdist
   (vela) $ conda env config vars set PYTHON_JULIACALL_HANDLE_SIGNALS=yes
   (vela) $ conda env config vars set PYTHON_JULIACALL_THREADS=4
   (vela) $ conda env config vars set JULIA_NUM_THREADS=4
   (vela) $ conda env config vars set JULIA_CONDAPKG_BACKEND="Null"
   (vela) $ julia
   julia> ] add LocalRegistry, JuliaFormatter, BenchmarkTools
   julia> ] registry add https://github.com/abhisrkckl/julia_registry
   julia> ] add https://github.com/abhisrkckl/Vela.jl
   julia> exit()
   (vela) $ pip install git+https://github.com/abhisrkckl/Vela.jl

Usage
-----
See ``run_example_nestle.py`` and ``run_example_emcee.py`` in ``pint2vela/examples`` directory
for examples on how to use ``Vela.jl`` with an MCMC sampler and a nested sampler respectively.
