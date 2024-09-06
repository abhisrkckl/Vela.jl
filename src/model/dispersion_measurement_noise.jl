export DispersionMeasurementNoise, dmefac, dmequad2

"""Modifications to the measured wideband DM uncertainties: DMEFACs and DMEQUADs.
Assumes that the DMEFACs & DMEQUADs are exclusive."""
struct DispersionMeasurementNoise <: WhiteNoiseComponent
    dmefac_index_mask::Vector{UInt}
    dmequad_index_mask::Vector{UInt}
end

"""The DMEFAC corresponding to a TOA (assumes exclusivity)."""
function dmefac(dwn::DispersionMeasurementNoise, cwtoa::CorrectedWidebandTOA, params::NamedTuple)
    idx = dwn.dmefac_index_mask[cwtoa.corrected_toa.toa.index]
    return (idx == 0) ? dimensionless(1.0) : params.DMEFAC[idx]
end

"""The DMEQUAD corresponding to a TOA (assumes exclusivity)."""
function dmequad2(dwn::DispersionMeasurementNoise, cwtoa::CorrectedWidebandTOA, params::NamedTuple)
    idx = dwn.dmequad_index_mask[cwtoa.corrected_toa.toa.index]
    return (idx == 0) ? GQ{2}(0.0) : params.DMEQUAD[idx]^Val(2)
end

"""Apply DMEFAC and DMEQUAD to a TOA."""
function correct_toa(dwn::DispersionMeasurementNoise, cwtoa::CorrectedWidebandTOA, params::NamedTuple)
    cdminfo = correct_dminfo(
        cwtoa.corrected_dminfo; 
        dmefac=dmefac(dwn, cwtoa, params),
        dmequad2=dmequad2(dwn, cwtoa, params),
    )
    return CorrectedWidebandTOA(cwtoa.corrected_toa, cdminfo)
end

function show(io::IO, dwn::DispersionMeasurementNoise)
    num_dmefacs = length(filter(x -> x > 0, unique(dwn.dmefac_index_mask)))
    num_dmequads = length(filter(x -> x > 0, unique(dwn.dmequad_index_mask)))
    print(io, "DispersionMeasurementNoise($num_dmefacs DMEFACs, $num_dmequads DMEQUADs)")
end
