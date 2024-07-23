export PhaseJump

"""System-dependent phase jumps."""
struct PhaseJump <: PhaseComponent
    jump_mask::BitMatrix
end

function phase(pjmp::PhaseJump, ctoa::CorrectedTOA, params::NamedTuple)::GQ
    @assert length(params.JUMP) == pjmp.jump_mask.dims[1]
    @assert ctoa.toa.index <= pjmp.jump_mask.dims[2]
    if ctoa.toa.tzr
        return dimensionless(0.0)
    else
        return sum(
            jj -> pjmp.jump_mask[jj, ctoa.toa.index] * params.JUMP[jj],
            1:length(params.JUMP),
        ) * (params.F_ + params.F[1])
    end
end

function show(io::IO, jmp::PhaseJump)
    num_jumps = jmp.jump_mask.dims[1]
    print(io, "PhaseJump($num_jumps JUMPs)")
end
show(io::IO, ::MIME"text/plain", jmp::PhaseJump) = show(io, jmp)
