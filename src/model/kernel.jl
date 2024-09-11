export WhiteNoiseKernel, EcorrKernel

"""Abstract base class of all likelihood kernels"""
abstract type Kernel end

"""A kernel representing only uncorrelated noise.
The covariance matrix is diagonal."""
struct WhiteNoiseKernel <: Kernel end

"""A kernel representing only time-uncorrelated noise.
The covariance matrix is block-diagonal."""
struct EcorrKernel <: Kernel end
