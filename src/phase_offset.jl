using GeometricUnits

struct PhaseOffset <: PhaseComponent end

phase(::PhaseOffset, toa::TOA, params::Dict)::GQ =
    toa.tzr ? dimensionless(0.0) : params["PHOFF"]
