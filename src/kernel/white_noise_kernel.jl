export Kernel, WhiteNoiseKernel

"""
    Kernel
    
Abstract base class of all likelihood kernels"""
abstract type Kernel end

show(io::IO, ::MIME"text/plain", kernel::Kernel) = show(io, kernel)

"""A kernel representing only uncorrelated noise.
The covariance matrix is diagonal.

Reference:
    [Hobbs+ 2006](http://doi.org/10.1111/j.1365-2966.2006.10302.x),
    [Alam+ 2021](http://doi.org/10.3847/1538-4365/abc6a1)
"""
struct WhiteNoiseKernel <: Kernel end
