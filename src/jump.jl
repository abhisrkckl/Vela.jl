export PhaseJump

"""System-dependent phase jumps."""
struct PhaseJump <: PhaseComponent
    jump_mask::Matrix{Float64}
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
show(io::IO, ::MIME"text/plain", jmp::PhaseJump) = show(io, jmp)
