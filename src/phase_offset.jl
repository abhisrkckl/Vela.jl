using GeometricUnits

export PhaseOffset, phase

struct PhaseOffset <: PhaseComponent end

read_params_from_dict(::PhaseOffset, params::Dict) = (PHOFF = params["PHOFF"][1],)

phase(::PhaseOffset, toa::TOA, params::NamedTuple)::GQ = -Int(toa.tzr) * params.PHOFF
