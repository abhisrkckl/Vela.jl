export DispersionMeasurementNoise, dmefac, dmequad2

"""Modifications to the measured wideband DM uncertainties: DMEFACs and DMEQUADs.
Assumes that the DMEFACs & DMEQUADs are exclusive.

Reference:
    [Alam+ 2021](http://doi.org/10.3847/1538-4365/abc6a1)
"""
struct DispersionMeasurementNoise <: WhiteNoiseComponent
    dmefac_index_mask::Vector{UInt}
    dmequad_index_mask::Vector{UInt}
end

"""The DMEFAC corresponding to a TOA (assumes exclusivity)."""
function dmefac(dwn::DispersionMeasurementNoise, wtoa::WidebandTOA, params::NamedTuple)
    idx = dwn.dmefac_index_mask[wtoa.toa.index]
    return (idx == 0) ? dimensionless(1.0) : params.DMEFAC[idx]
end

"""The DMEQUAD corresponding to a TOA (assumes exclusivity)."""
function dmequad2(dwn::DispersionMeasurementNoise, wtoa::WidebandTOA, params::NamedTuple)
    idx = dwn.dmequad_index_mask[wtoa.toa.index]
    return (idx == 0) ? GQ{-2}(0.0) : params.DMEQUAD[idx]^Val(2)
end

"""Apply DMEFAC and DMEQUAD to a TOA."""
function correct_toa(
    dwn::DispersionMeasurementNoise,
    wtoa::WidebandTOA,
    wtoacorr::WidebandTOACorrection,
    params::NamedTuple,
)
    dmcorr = correct_dminfo(
        wtoacorr.dm_correction;
        dmefac = dmefac(dwn, wtoa, params),
        dmequad2 = dmequad2(dwn, wtoa, params),
    )
    return WidebandTOACorrection(wtoacorr.toa_correction, dmcorr)
end

correct_toa(::DispersionMeasurementNoise, ::TOA, toacorr::TOACorrection, ::NamedTuple) =
    correct_toa(toacorr)

function show(io::IO, dwn::DispersionMeasurementNoise)
    num_dmefacs = length(filter(x -> x > 0, unique(dwn.dmefac_index_mask)))
    num_dmequads = length(filter(x -> x > 0, unique(dwn.dmequad_index_mask)))
    print(io, "DispersionMeasurementNoise($num_dmefacs DMEFACs, $num_dmequads DMEQUADs)")
end
