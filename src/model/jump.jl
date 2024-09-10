export PhaseJump, ExclusivePhaseJump

abstract type PhaseJumpBase <: PhaseComponent end

"""System-dependent phase jumps with non-exclusive selection masks.

Reference:
    [Hobbs+ 2006](http://doi.org/10.1111/j.1365-2966.2006.10302.x)
"""
struct PhaseJump <: PhaseJumpBase
    jump_mask::BitMatrix
end

phase(pjmp::PhaseJump, ctoa::CorrectedTOA, params::NamedTuple)::GQ =
    ctoa.toa.tzr ? dimensionless(0.0) :
    (basis_dot(pjmp.jump_mask, params.JUMP, ctoa.toa.index) * (params.F_ + params.F[1]))

function show(io::IO, jmp::PhaseJump)
    num_jumps = size(jmp.jump_mask)[1]
    print(io, "PhaseJump($num_jumps JUMPs)")
end

"""System-dependent phase jumps with exclusive selection masks."""
struct ExclusivePhaseJump <: PhaseJumpBase
    jump_mask::Vector{UInt}
end

function phase(pjmp::ExclusivePhaseJump, ctoa::CorrectedTOA, params::NamedTuple)::GQ
    return if (ctoa.toa.index == 0 || ctoa.toa.tzr)
        dimensionless(0.0)
    else
        idx = pjmp.jump_mask[ctoa.toa.index]
        if idx == 0
            dimensionless(0.0)
        else
            params.JUMP[idx] * (params.F_ + params.F[1])
        end
    end
end

function show(io::IO, jmp::ExclusivePhaseJump)
    num_jumps = length(filter(x -> x > 0, unique(jmp.jump_mask)))
    print(io, "PhaseJump($num_jumps JUMPs, exclusive)")
end
