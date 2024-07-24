"""A package for evaluating the single-pulsar timing & noise model."""
module Vela

using GeometricUnits
using DoubleFloats: Double64
using LinearAlgebra: dot
using .Threads: @spawn, fetch, nthreads
using Unrolled: @unroll
import HDF5
import JSON

include("toa/ephemeris.jl")
include("toa/toa.jl")
include("parameter/parameter.jl")
include("model/component.jl")
include("model/spindown.jl")
include("model/phase_offset.jl")
include("model/jump.jl")
include("model/solarsystem.jl")
include("model/troposphere.jl")
include("model/dispersion.jl")
include("model/solarwind.jl")
include("model/measurement_noise.jl")
include("model/timing_model.jl")
include("residuals.jl")
include("chi2.jl")
include("likelihood.jl")
include("read_model_and_toas.jl")
include("summary_plot.jl")
end
