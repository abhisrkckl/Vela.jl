export MeasurementNoise, efac, equad2

"""Modifications to the measured TOA uncertainties: 
EFACs (error factors) and EQUADs (errors added in quadrature).
Assumes that the EFACs & EQUADs are are exclusive."""
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
    return (idx == 0) ? GQ(0.0, 2) : params.EQUAD[idx]^2
end

correct_toa(wn::MeasurementNoise, ctoa::CorrectedTOA, params::NamedTuple) =
    correct_toa(ctoa; efac = efac(wn, ctoa, params), equad2 = equad2(wn, ctoa, params))
