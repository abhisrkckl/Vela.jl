export DispersionTaylor, dispersion_slope

struct DispersionTaylor <: DispersionComponent
    number_of_terms::UInt
end

dispersion_slope(::DispersionTaylor, ::TOA, params) = GQ(0.0, -1)
