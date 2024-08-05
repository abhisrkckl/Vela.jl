export PhaseJump, ExclusivePhaseJump

abstract type PhaseJumpBase <: PhaseComponent end

"""System-dependent phase jumps with non-exclusive selection masks."""
struct PhaseJump <: PhaseJumpBase
    jump_mask::BitMatrix
end

function phase(pjmp::PhaseJump, ctoa::CorrectedTOA, params::NamedTuple)::GQ
    @assert length(params.JUMP) == size(pjmp.jump_mask)[1]
    @assert ctoa.toa.index <= size(pjmp.jump_mask)[2]
    if ctoa.toa.tzr
        return dimensionless(0.0)
    else
        return dot(params.JUMP, @view(pjmp.jump_mask[:, ctoa.toa.index])) *
               (params.F_ + params.F[1])
    end
end

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
    num_jumps = length(unique(jmp.jump_mask))
    print(io, "PhaseJump($num_jumps JUMPs, exclusive)")
end
