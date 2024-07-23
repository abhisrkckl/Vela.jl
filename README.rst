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
   (vela) $ conda env config vars set PYTHON_JULIACALL_HANDLE_SIGNALS=yes
   (vela) $ conda env config vars set PYTHON_JULIACALL_THREADS=4
   (vela) $ julia
   julia> ] add LocalRegistry
   julia> ] registry add https://github.com/abhisrkckl/julia_registry
   julia> ] add https://github.com/abhisrkckl/Vela.jl
   julia> exit()
   (vela) $ pip install git+https://github.com/abhisrkckl/pint2vela

Using from Python
-----------------
Here is what an example session looks like::

   > from pint2vela import read_model_and_toas, vl
   > 
   > vl = jl.Vela
   > mv, tv = read_model_and_toas("J1234+5678.par", "J1234+5678.tim")
   > lnlike = vl.get_lnlike_func(mv, tv)
