"""`Vela.jl` is a package for doing Bayesian single-pulsar timing and noise analysis."""
module Vela

using GeometricUnits
using DoubleFloats: Double64
using LinearAlgebra: dot
using .Threads: @spawn, fetch, nthreads
using Unrolled: @unroll
import HDF5
import JSON

include("ephemeris.jl")
include("toa.jl")
include("parameter.jl")
include("component.jl")
include("spindown.jl")
include("phase_offset.jl")
include("solarsystem.jl")
include("troposphere.jl")
include("dispersion.jl")
include("solarwind.jl")
include("measurement_noise.jl")
include("timing_model.jl")
include("residuals.jl")
include("chi2.jl")
include("likelihood.jl")
include("read_model_and_toas.jl")
include("summary_plot.jl")
end
