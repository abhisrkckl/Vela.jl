using GeometricUnits

export PhaseOffset, phase

"""Phase offset between measured TOAs and the TZR TOA."""
struct PhaseOffset <: PhaseComponent end

read_params_from_dict(::PhaseOffset, params::Dict) = (PHOFF = params[:PHOFF][1],)

phase(::PhaseOffset, ctoa::CorrectedTOA, params::NamedTuple)::GQ =
    -Int(!ctoa.toa.tzr) * params.PHOFF
