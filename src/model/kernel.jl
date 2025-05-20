export Kernel,
    WhiteNoiseKernel, EcorrKernel, EcorrGroup, WoodburyKernel, get_marginalized_param_names

"""
    Kernel
    
Abstract base class of all likelihood kernels"""
abstract type Kernel end

"""
    WhiteNoiseKernel

A kernel representing only uncorrelated noise.
The covariance matrix is diagonal.

Reference:
    [Hobbs+ 2006](http://doi.org/10.1111/j.1365-2966.2006.10302.x),
    [Alam+ 2021](http://doi.org/10.3847/1538-4365/abc6a1)
"""
struct WhiteNoiseKernel <: Kernel end

get_marginalized_param_names(::WhiteNoiseKernel) = String[]

"""Range of TOAs belonging to an ECORR block."""
struct EcorrGroup
    start::UInt
    stop::UInt
    index::UInt
end

"""
    EcorrKernel

A kernel representing white noise and ECORR.
The covariance matrix is block-diagonal.
Assumes that the `TOA`s are sorted in the correct order.
Not applicable for wideband TOAs.

Reference:
    [Johnson+ 2024](https://doi.org/10.1103/PhysRevD.109.103012)
"""
struct EcorrKernel <: Kernel
    ecorr_groups::Vector{EcorrGroup}
end

get_marginalized_param_names(::EcorrKernel) = String[]

show(io::IO, ::MIME"text/plain", model::Kernel) = show(io, model)
function show(io::IO, ek::EcorrKernel)
    ecorr_idxs = unique(grp.index for grp in ek.ecorr_groups)
    noecorr_flag = 0 in ecorr_idxs
    print(
        io,
        "EcorrKernel($(length(ecorr_idxs) - Int(noecorr_flag)) ECORRs, $(length(ek.ecorr_groups)) groups)",
    )
end

"""
    WoodburyKernel

A kernel representing white noise and correlated noise including ECORR.

This type has an `inner_kernel` attribute which represents the time-uncorrelated
part of the timing noise. It can be a `WhiteNoiseKernel` if the only time-uncorrelated
noise is white noise, or `EcorrKernel` if ECORR noise is also present.

The `gp_components` attribute contains a collection of amplitude-marginalized Gaussian 
noise components. These are treated as part of the covariance matrix. `WoodburyKernel.gp_components`
and `TimingModel.components` must have no common elements.

Reference:
    [Johnson+ 2024](https://doi.org/10.1103/PhysRevD.109.103012)
"""
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

get_marginalized_param_names(wk::WoodburyKernel) =
    vcat([get_marginalized_param_names(gp) for gp in wk.gp_components]...)

show(io::IO, wk::WoodburyKernel) =
    print(io, "WoodburyKernel($(wk.inner_kernel), $(wk.gp_components))")
