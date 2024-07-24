export PhaseOffset, phase

"""Phase offset between measured TOAs and the TZR TOA.

Corresponds to `PhaseOffset` in `PINT`."""
struct PhaseOffset <: PhaseComponent end

"""Phase correction due to the overall phase offset."""
phase(::PhaseOffset, ctoa::CorrectedTOA, params::NamedTuple)::GQ =
    -Int(!ctoa.toa.tzr) * params.PHOFF
