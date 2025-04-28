export DispersionJump, ExclusiveDispersionJump

abstract type DispersionJumpBase <: DispersionComponent end

delay(::DispersionJumpBase, ::TOA, ::TOACorrection, ::NamedTuple) = time(0.0)

function correct_toa(
    dmjump::DispersionJumpBase,
    wtoa::WidebandTOA,
    wtoacorr::WidebandTOACorrection,
    params::NamedTuple,
)
    dm = dispersion_slope(dmjump, wtoa.toa, wtoacorr.toa_correction, params)
    dmcorr = correct_dminfo(wtoacorr.dm_correction; delta_dm = dm)
    return WidebandTOACorrection(wtoacorr.toa_correction, dmcorr)
end

"""
    DispersionJump
    
System-dependent wideband DM offsets (`DMJUMP`) with non-exclusive selection masks.

Unlike an `FDJUMPDM`, a `DMJUMP` only provides a DM correction and no delay.
Corresponds to `DispersionJump` in `PINT`.

Reference:
    [Alam+ 2021](http://doi.org/10.3847/1538-4365/abc6a1)
"""
struct DispersionJump <: DispersionJumpBase
    jump_mask::BitMatrix
end

function dispersion_slope(
    dmjump::DispersionJump,
    toa::TOA,
    ::TOACorrection,
    params::NamedTuple,
)::GQ
    @assert length(params.DMJUMP) == size(dmjump.jump_mask)[1]
    @assert toa.index <= size(dmjump.jump_mask)[2]
    if is_tzr(toa)
        return GQ{-1}(0.0)
    else
        return -dot(params.DMJUMP, @view(dmjump.jump_mask[:, toa.index]))
    end
end

function show(io::IO, dmjump::DispersionJump)
    num_jumps = size(dmjump.jump_mask)[1]
    print(io, "DispersionJump($num_jumps DMJUMPs)")
end

"""
    ExclusiveDispersionJump

System-dependent wideband DM offsets (`DMJUMP`) with exclusive selection masks.
Unlike an `FDJUMPDM`, a `DMJUMP` only provides a DM correction and no delay."""
struct ExclusiveDispersionJump <: DispersionJumpBase
    jump_mask::Vector{UInt}
end

function dispersion_slope(
    dmjump::ExclusiveDispersionJump,
    toa::TOA,
    ::TOACorrection,
    params::NamedTuple,
)::GQ
    return if (is_tzr(toa))
        GQ{-1}(0.0)
    else
        idx = dmjump.jump_mask[toa.index]
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
