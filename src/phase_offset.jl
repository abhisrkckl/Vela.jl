using GeometricUnits
using Distributions

export PhaseOffset, phase, lnprior

"""Phase offset between measured TOAs and the TZR TOA."""
struct PhaseOffset <: PhaseComponent end

phase(::PhaseOffset, ctoa::CorrectedTOA, params::NamedTuple)::GQ =
    -Int(!ctoa.toa.tzr) * params.PHOFF

"""Physically motivated prior: PHOFF ∈ (-0.5, 0.5]"""
lnprior(::PhaseOffset, params::NamedTuple) = logpdf(Uniform(-0.5, 0.5), value(params.PHOFF))