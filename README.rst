Vela
----

.. image:: https://github.com/abhisrkckl/Vela.jl/actions/workflows/CI.yml/badge.svg
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
::

   (base) $ conda create -n vela python=3.11
   (base) $ conda activate vela
   (vela) $ conda install -c conda-forge julia juliacall pint-pulsar  
   (vela) $ julia
   julia> ] add LocalRegistry
   julia> ] registry add https://github.com/abhisrkckl/julia_registry
   julia> ] add https://github.com/abhisrkckl/Vela.jl

Using from Python
-----------------
We need to set the following environment variables first::

   export PYTHON_JULIACALL_HANDLE_SIGNALS=yes
   export PYTHON_JULIACALL_THREADS=4

Here is what an example session looks like::

   > from juliacall import Main as jl
   > jl.seval("Using Vela")
   > vl = jl.Vela
   > mv, tv = vl.read_model_and_toas("J1234+5678.hdf5")
   > lnlike_s = vl.get_lnlike_serial_func(mv, tv)
   > lnlike_p = vl.get_lnlike_parallel_func(mv, tv)