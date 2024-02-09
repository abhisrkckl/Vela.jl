using GeometricUnits
using Quadmath

struct Spindown <: PhaseComponent
    number_of_terms::UInt
end

phase(spindown::Spindown, toa::TOA, params) = dimensionless(0)
