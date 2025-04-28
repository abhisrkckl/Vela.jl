export PhaseOffset, phase

"""Phase offset between measured TOAs and the TZR TOA.
Corresponds to `PhaseOffset` in `PINT`.

Reference:
    [Susobhanan+ 2024](http://doi.org/10.3847/1538-4357/ad59f7)
"""
struct PhaseOffset <: PhaseComponent end

"""Phase correction due to the overall phase offset."""
phase(::PhaseOffset, toa::TOA, ::TOACorrection, params::NamedTuple) =
    -Int(!is_tzr(toa)) * params.PHOFF
