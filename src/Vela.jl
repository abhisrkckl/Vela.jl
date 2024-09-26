"""`Vela.jl` is a package for doing Bayesian single-pulsar timing and noise analysis."""
module Vela

using GeometricUnits
using DoubleFloats: Double64
using LinearAlgebra: dot
using .Threads: @threads, @spawn, fetch, nthreads
using Unrolled: @unroll
import Distributions:
    Distribution,
    ContinuousUnivariateDistribution,
    pdf,
    logpdf,
    cdf,
    logcdf,
    quantile,
    support,
    minimum,
    maximum,
    insupport,
    RealInterval
import JLSO

include("toa/ephemeris.jl")
include("toa/toa.jl")
include("toa/wideband_toa.jl")
include("parameter/parameter.jl")
include("model/component.jl")
include("model/kernel.jl")
include("model/spindown.jl")
include("model/phase_offset.jl")
include("model/glitch.jl")
include("model/jump.jl")
include("model/solarsystem.jl")
include("model/troposphere.jl")
include("model/dispersion.jl")
include("model/fdjumpdm.jl")
include("model/dmjump.jl")
include("model/chromatic.jl")
include("model/solarwind.jl")
include("model/binary/orbit.jl")
include("model/binary/binary_ell1_base.jl")
include("model/binary/binary_ell1.jl")
include("model/binary/binary_ell1h.jl")
include("model/binary/binary_dd_base.jl")
include("model/binary/binary_dd.jl")
include("model/binary/binary_ddh.jl")
include("model/binary/binary_dds.jl")
include("model/binary/binary_ddk.jl")
include("model/frequency_dependent.jl")
include("model/wavex.jl")
include("model/measurement_noise.jl")
include("model/dispersion_measurement_noise.jl")
include("model/timing_model.jl")
include("model/wideband_model.jl")
include("prior/prior.jl")
include("prior/simple_prior.jl")
include("prior/known_priors.jl")
include("residuals/residuals.jl")
include("residuals/wideband_residuals.jl")
include("likelihood/wls_chi2.jl")
include("likelihood/ecorr_chi2.jl")
include("likelihood/wideband_wls_chi2.jl")
include("likelihood/wls_likelihood.jl")
include("likelihood/ecorr_likelihood.jl")
include("likelihood/wideband_wls_likelihood.jl")
include("likelihood/posterior.jl")
include("readwrite/readwrite.jl")

end
