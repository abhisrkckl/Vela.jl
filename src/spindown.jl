using GeometricUnits
using Quadmath

export Spindown, phase

struct Spindown <: PhaseComponent
    number_of_terms::UInt
end

phase(spindown::Spindown, toa::TOA, params) = dimensionless(0.0)
