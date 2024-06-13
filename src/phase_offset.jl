using GeometricUnits

export PhaseOffset, phase

"""Phase offset between measured TOAs and the TZR TOA."""
struct PhaseOffset <: PhaseComponent end

phase(::PhaseOffset, ctoa::CorrectedTOA, params::NamedTuple)::GQ =
    -Int(!ctoa.toa.tzr) * params.PHOFF
