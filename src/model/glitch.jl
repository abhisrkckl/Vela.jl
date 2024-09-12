export Glitch

"""Glitches in pulsar rotation.

Reference:
    [Hobbs+ 2006](http://doi.org/10.1111/j.1365-2966.2006.10302.x)
"""
struct Glitch <: PhaseComponent end

function glitch_phase(t, t_gl, ϕ_gl, f0_gl, f1_gl, f2_gl, fd_gl, τ_gl)
    dt = t - t_gl
    return (t < t_gl) ? dimensionless(0.0) :
           (
        ϕ_gl +
        dt * (f0_gl + dt * f1_gl / 2 + dt * dt * f2_gl / 6) +
        fd_gl * τ_gl * (1 - exp(-dt / τ_gl))
    )
end

function phase(::Glitch, ctoa::CorrectedTOA, params::NamedTuple)::GQ
    t = corrected_toa_value(ctoa)
    return sum(
        glitch_phase(t, t_gl, ϕ_gl, f0_gl, f1_gl, f2_gl, fd_gl, τ_gl) for
        (t_gl, ϕ_gl, f0_gl, f1_gl, f2_gl, fd_gl, τ_gl) in zip(
            params.GLEP_,
            params.GLPH_,
            params.GLF0_,
            params.GLF1_,
            params.GLF2_,
            params.GLF0D_,
            params.GLTD_,
        )
    )
end
