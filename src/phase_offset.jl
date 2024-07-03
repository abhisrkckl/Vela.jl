export PhaseOffset, phase, lnprior

"""Phase offset between measured TOAs and the TZR TOA."""
struct PhaseOffset <: PhaseComponent end

phase(::PhaseOffset, ctoa::CorrectedTOA, params::NamedTuple)::GQ =
    -Int(!ctoa.toa.tzr) * params.PHOFF

const prior_PHOFF = Uniform(-0.5, 0.5)

"""Physically motivated prior: PHOFF ∈ (-0.5, 0.5]"""
lnprior(::PhaseOffset, params::NamedTuple) = logpdf(prior_PHOFF, value(params.PHOFF))
