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
a::
   (base) $ conda create -n vela python=3.11
   (base) $ conda activate vela
   (vela) $ conda install -c conda-forge julia pint-pulsar  
   (vela) $ julia
   julia> ] add LocalRegistry
   julia> ] registry add https://github.com/abhisrkckl/julia_registry
   julia> ] add https://github.com/abhisrkckl/Vela.jl
