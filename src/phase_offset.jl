using GeometricUnits

export PhaseOffset, phase

struct PhaseOffset <: PhaseComponent end

phase(::PhaseOffset, toa::TOA, params::NamedTuple)::GQ = -Int(!toa.tzr) * params.PHOFF
