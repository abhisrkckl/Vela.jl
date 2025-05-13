"""`Vela.jl` is a package for doing Bayesian single-pulsar timing and noise analysis."""
module Vela

using GeometricUnits
using DoubleFloats: Double64
using LinearAlgebra: dot, Symmetric, cholesky!, ldiv!, logdet, isposdef
using .Threads: @threads, @spawn, fetch, nthreads
using Unrolled: @unroll
using Distributions
import Distributions:
    pdf, logpdf, cdf, logcdf, quantile, support, minimum, maximum, insupport
import JLSO
import PkgVersion

export GQ

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
include("model/dispersion.jl")
include("model/fdjumpdm.jl")
include("model/dmjump.jl")
include("model/chromatic.jl")
include("model/solarwind.jl")
include("model/binary/orbit.jl")
include("model/binary/binary_ell1_base.jl")
include("model/binary/binary_ell1.jl")
include("model/binary/binary_ell1h.jl")
include("model/binary/binary_ell1k.jl")
include("model/binary/binary_dd_base.jl")
include("model/binary/binary_dd.jl")
include("model/binary/binary_ddh.jl")
include("model/binary/binary_dds.jl")
include("model/binary/binary_ddk.jl")
include("model/frequency_dependent.jl")
include("model/wavex.jl")
include("model/measurement_noise.jl")
include("model/gp_noise.jl")
include("model/dispersion_measurement_noise.jl")
include("model/marginalized_timing_model.jl")
include("model/timing_model.jl")
include("model/wideband_model.jl")
include("prior/prior.jl")
include("prior/prior_scaling.jl")
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
include("likelihood/gls_likelihood.jl")
include("likelihood/gls_ecorr_likelihood.jl")
include("likelihood/posterior.jl")
include("pulsar/pulsar.jl")
include("readwrite/readwrite.jl")

pkg_version() = string(PkgVersion.Version(Vela))

end
