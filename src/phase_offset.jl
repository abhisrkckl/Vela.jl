using GeometricUnits
using Distributions

export PhaseOffset, phase

"""Phase offset between measured TOAs and the TZR TOA."""
struct PhaseOffset <: PhaseComponent end

phase(::PhaseOffset, ctoa::CorrectedTOA, params::NamedTuple)::GQ =
    -Int(!ctoa.toa.tzr) * params.PHOFF

prior_distributions(::PhaseOffset) = (
    SimplePrior(:PHOFF, Uniform(-0.5, 0.5))
)