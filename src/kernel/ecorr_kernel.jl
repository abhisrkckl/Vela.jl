export EcorrKernel, EcorrGroup

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

show(io::IO, ek::EcorrKernel) = print(
    io,
    "EcorrKernel($(length(unique(grp.index for grp in ek.ecorr_groups)) - 1) ECORRs, $(length(ek.ecorr_groups)) groups)",
)
