export MeasurementNoise, efac, equad2

"""Modifications to the measured TOA uncertainties: 
EFACs (error factors) and EQUADs (errors added in quadrature).
Assumes that the EFACs & EQUADs are exclusive.

Reference:
    [Lentati+ 2014](http://doi.org/10.1093/mnras/stt2122)
"""
struct MeasurementNoise <: WhiteNoiseComponent
    efac_index_mask::Vector{UInt}
    equad_index_mask::Vector{UInt}
end

"""The EFAC corresponding to a TOA (assumes exclusivity)."""
function efac(wn::MeasurementNoise, ctoa::CorrectedTOA, params::NamedTuple)
    idx = wn.efac_index_mask[ctoa.toa.index]
    return (idx == 0) ? dimensionless(1.0) : params.EFAC[idx]
end

"""The EQUAD corresponding to a TOA (assumes exclusivity)."""
function equad2(wn::MeasurementNoise, ctoa::CorrectedTOA, params::NamedTuple)
    idx = wn.equad_index_mask[ctoa.toa.index]
    return (idx == 0) ? GQ{2}(0.0) : params.EQUAD[idx]^Val(2)
end

"""Apply EFAC and EQUAD to a TOA."""
correct_toa(wn::MeasurementNoise, ctoa::CorrectedTOA, params::NamedTuple) =
    is_tzr(ctoa.toa) ? correct_toa(ctoa) :
    correct_toa(ctoa; efac = efac(wn, ctoa, params), equad2 = equad2(wn, ctoa, params))

function show(io::IO, wn::MeasurementNoise)
    num_efacs = length(filter(x -> x > 0, unique(wn.efac_index_mask)))
    num_equads = length(filter(x -> x > 0, unique(wn.equad_index_mask)))
    print(io, "MeasurementNoise($num_efacs EFACs, $num_equads EQUADs)")
end
