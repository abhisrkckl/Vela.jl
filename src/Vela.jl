module Vela

using GeometricUnits
using Quadmath: Float128
using LinearAlgebra: dot
using Distributions
using .Threads: @threads, Atomic, atomic_add!
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
# include("selection.jl")
include("timing_model.jl")
include("priors.jl")
include("read_model_and_toas.jl")
include("summary_plot.jl")

end
