export PhaseOffset, phase

"""Phase offset between measured TOAs and the TZR TOA.

Corresponds to `PhaseOffset` in `PINT`.

Reference:
    [Susobhanan+ 2024](http://doi.org/10.3847/1538-4357/ad59f7)
"""
struct PhaseOffset <: PhaseComponent end

"""Phase correction due to the overall phase offset."""
phase(::PhaseOffset, ctoa::CorrectedTOA, params::NamedTuple)::GQ =
    -Int(!ctoa.toa.tzr) * params.PHOFF
