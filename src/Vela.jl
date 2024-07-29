"""A package for evaluating the single-pulsar timing & noise model."""
module Vela

using GeometricUnits
using DoubleFloats: Double64
using LinearAlgebra: dot
using .Threads: @spawn, fetch, nthreads
using Unrolled: @unroll
using Distributions: Distribution, logpdf, quantile
import JLSO

include("toa/ephemeris.jl")
include("toa/toa.jl")
include("parameter/parameter.jl")
include("model/component.jl")
include("model/basis.jl")
include("model/spindown.jl")
include("model/phase_offset.jl")
include("model/jump.jl")
include("model/solarsystem.jl")
include("model/troposphere.jl")
include("model/dispersion.jl")
include("model/solarwind.jl")
include("model/wavex.jl")
include("model/measurement_noise.jl")
include("model/timing_model.jl")
include("prior/prior.jl")
include("prior/simple_prior.jl")
include("residuals/residuals.jl")
include("likelihood/chi2.jl")
include("likelihood/likelihood.jl")
include("readwrite/readwrite.jl")
end
