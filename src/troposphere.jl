export Troposphere, delay

struct Troposphere <: DelayComponent end

delay(::Troposphere, toa::TOA, params) = time(0.0)