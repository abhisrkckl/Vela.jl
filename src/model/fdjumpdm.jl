export DispersionOffset, ExclusiveDispersionOffset

abstract type DispersionOffsetBase <: DispersionComponent end

"""System-dependent DM offsets (`FDJUMPDM`) with non-exclusive selection masks.

Reference:
    [Susobhanan+ 2024](http://doi.org/10.3847/1538-4357/ad59f7)
"""
struct DispersionOffset <: DispersionOffsetBase
    jump_mask::BitMatrix
end

function dispersion_slope(
    dmoff::DispersionOffset,
    ctoa::CorrectedTOA,
    params::NamedTuple,
)::GQ
    @assert length(params.FDJUMPDM) == size(dmoff.jump_mask)[1]
    @assert ctoa.toa.index <= size(dmoff.jump_mask)[2]
    if is_tzr(ctoa.toa)
        return GQ{-1}(0.0)
    else
        return -dot(params.FDJUMPDM, @view(dmoff.jump_mask[:, ctoa.toa.index]))
    end
end

function show(io::IO, dmoff::DispersionOffset)
    num_jumps = size(dmoff.jump_mask)[1]
    print(io, "DispersionOffset($num_jumps FDJUMPDMs)")
end

"""System-dependent DM offsets (`FDJUMPDM`) with exclusive selection masks."""
struct ExclusiveDispersionOffset <: DispersionOffsetBase
    jump_mask::Vector{UInt}
end

function dispersion_slope(
    dmoff::ExclusiveDispersionOffset,
    ctoa::CorrectedTOA,
    params::NamedTuple,
)::GQ
    return if is_tzr(ctoa.toa)
        GQ{-1}(0.0)
    else
        idx = dmoff.jump_mask[ctoa.toa.index]
        if idx == 0
            GQ{-1}(0.0)
        else
            -params.FDJUMPDM[idx]
        end
    end
end

function show(io::IO, dmoff::ExclusiveDispersionOffset)
    num_jumps = length(filter(x -> x > 0, unique(dmoff.jump_mask)))
    print(io, "DispersionOffset($num_jumps FDJUMPDMs, exclusive)")
end
