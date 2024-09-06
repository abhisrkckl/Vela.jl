export DispersionJump, ExclusiveDispersionJump

abstract type DispersionJumpBase <: DispersionComponent end

"""System-dependent wideband DM offsets (`DMJUMP`) with non-exclusive selection masks.
Unlike an `FDJUMPDM`, a `DMJUMP` only provides a DM correction and no delay."""
struct DispersionJump <: DispersionJumpBase
    jump_mask::BitMatrix
end

function dispersion_slope(
    dmjump::DispersionJump,
    ctoa::CorrectedTOA,
    params::NamedTuple,
)::GQ
    @assert length(params.DMJUMP) == size(dmjump.jump_mask)[1]
    @assert ctoa.toa.index <= size(dmjump.jump_mask)[2]
    if ctoa.toa.tzr
        return GQ{-1}(0.0)
    else
        return -dot(params.DMJUMP, @view(dmjump.jump_mask[:, ctoa.toa.index]))
    end
end

function show(io::IO, dmjump::DispersionJump)
    num_jumps = size(dmjump.jump_mask)[1]
    print(io, "DispersionJump($num_jumps DMJUMPs)")
end

"""System-dependent wideband DM offsets (`DMJUMP`) with exclusive selection masks.
Unlike an `FDJUMPDM`, a `DMJUMP` only provides a DM correction and no delay."""
struct ExclusiveDispersionJump <: DispersionJumpBase
    jump_mask::Vector{UInt}
end

function dispersion_slope(
    dmjump::ExclusiveDispersionJump,
    ctoa::CorrectedTOA,
    params::NamedTuple,
)::GQ
    return if (ctoa.toa.index == 0 || ctoa.toa.tzr)
        GQ{-1}(0.0)
    else
        idx = dmjump.jump_mask[ctoa.toa.index]
        if idx == 0
            GQ{-1}(0.0)
        else
            -params.DMJUMP[idx]
        end
    end
end

function show(io::IO, dmjump::ExclusiveDispersionJump)
    num_jumps = length(filter(x -> x > 0, unique(dmjump.jump_mask)))
    print(io, "DispersionJump($num_jumps DMJUMPs, exclusive)")
end
