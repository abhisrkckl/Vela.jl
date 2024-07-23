export PhaseJump

"""System-dependent phase jumps."""
struct PhaseJump <: PhaseComponent
    jump_mask::BitMatrix
end

phase(pjmp::PhaseJump, ctoa::CorrectedTOA, params::NamedTuple)::GQ =
    ctoa.toa.tzr ? dimensionless(0.0) :
    (basis_dot(pjmp.jump_mask, params.JUMP, ctoa.toa.index) * (params.F_ + params.F[1]))

function show(io::IO, jmp::PhaseJump)
    num_jumps = size(jmp.jump_mask)[1]
    print(io, "PhaseJump($num_jumps JUMPs)")
end
show(io::IO, ::MIME"text/plain", jmp::PhaseJump) = show(io, jmp)
