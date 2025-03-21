export Kernel, WhiteNoiseKernel, EcorrKernel, EcorrGroup, WoodburyKernel

"""
    Kernel
    
Abstract base class of all likelihood kernels"""
abstract type Kernel end

"""A kernel representing only uncorrelated noise.
The covariance matrix is diagonal.

Reference:
    [Hobbs+ 2006](http://doi.org/10.1111/j.1365-2966.2006.10302.x),
    [Alam+ 2021](http://doi.org/10.3847/1538-4365/abc6a1)
"""
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

Reference:
    [Johnson+ 2024](https://doi.org/10.1103/PhysRevD.109.103012)
"""
struct EcorrKernel <: Kernel
    ecorr_groups::Vector{EcorrGroup}
end

show(io::IO, ::MIME"text/plain", model::Kernel) = show(io, model)
show(io::IO, ek::EcorrKernel) = print(
    io,
    "EcorrKernel($(length(unique(grp.index for grp in ek.ecorr_groups)) - 1) ECORRs, $(length(ek.ecorr_groups)) groups)",
)

struct WoodburyKernel{InnerKernel<:Kernel,GPComponentsTuple<:Tuple} <: Kernel
    inner_kernel::InnerKernel
    gp_components::GPComponentsTuple
    noise_basis::Matrix{Float64}

    function WoodburyKernel(
        inner_kernel::Kernel,
        gp_components::Tuple,
        noise_basis::Matrix{Float64},
    )
        @assert all(is_gp_noise.(gp_components))
        @assert sum(get_gp_npars.(gp_components)) == size(noise_basis)[2]
        return new{typeof(inner_kernel),typeof(gp_components)}(
            inner_kernel,
            gp_components,
            noise_basis,
        )
    end
end

show(io::IO, wk::WoodburyKernel) =
    print(io, "WoodburyKernel($(wk.inner_kernel), $(wk.gp_components))")
