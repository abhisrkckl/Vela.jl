export WhiteNoiseKernel, EcorrKernel

"""Abstract base class of all likelihood kernels"""
abstract type Kernel end

"""A kernel representing only uncorrelated noise.
The covariance matrix is diagonal."""
struct WhiteNoiseKernel <: Kernel end

"""Range of TOAs belonging to an ECORR block."""
struct EcorrGroup
    start::UInt
    stop::UInt
    index::UInt
end

"""A kernel representing white noise and ECORR.
The covariance matrix is block-diagonal.

Assumes that the `TOA`s are sorted in the correct order.
"""
struct EcorrKernel <: Kernel
    ecorr_groups::Vector{EcorrGroup}
end
