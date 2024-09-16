export DMInfo,
    CorrectedDMInfo,
    correct_dminfo,
    dm_residual,
    scaled_dm_error_sqr,
    WidebandTOA,
    CorrectedWidebandTOA

"""
    DMInfo

DM information associated with a wideband TOA.

References:
    [Pennucci+ 2014](http://doi.org/10.1088/0004-637X/790/2/93),
    [Pennucci 2019](http://doi.org/10.3847/1538-4357/aaf6ef)
"""
struct DMInfo
    value::GQ{-1,Float64}
    error::GQ{-1,Float64}
end

"""
    CorrectedDMInfo
The accumulated timing & noise model corrections applied to wideband DM.
"""
struct CorrectedDMInfo
    dminfo::DMInfo
    model_dm::GQ{-1,Float64}
    dmefac::GQ{0,Float64}
    dmequad2::GQ{-2,Float64}
end

CorrectedDMInfo(dminfo::DMInfo) =
    CorrectedDMInfo(dminfo, GQ{-1}(0.0), dimensionless(1.0), GQ{-2}(0.0))

"""
    correct_dminfo(
        cdminfo::CorrectedDMInfo;
        delta_dm = GQ{-1}(0.0),
        dmefac = dimensionless(1.0),
        dmequad2 = GQ{-2}(0.0),
    )
    
Apply a correction to a `DMInfo`.
"""
correct_dminfo(
    cdminfo::CorrectedDMInfo;
    delta_dm = GQ{-1}(0.0),
    dmefac = dimensionless(1.0),
    dmequad2 = GQ{-2}(0.0),
) = CorrectedDMInfo(
    cdminfo.dminfo,
    cdminfo.model_dm + delta_dm,
    cdminfo.dmefac * dmefac,
    cdminfo.dmequad2 + dmequad2,
)

"""
    dm_residual(cdminfo::CorrectedDMInfo)

Compute the DM residual.
"""
dm_residual(cdminfo::CorrectedDMInfo) = cdminfo.dminfo.value - cdminfo.model_dm

"""
    scaled_dm_error_sqr(cdminfo::CorrectedDMInfo)

Compute the scaled DM uncertainty incorporating the DMEFAC and the DMEQUAD.
"""
scaled_dm_error_sqr(cdminfo::CorrectedDMInfo) =
    (cdminfo.dminfo.error * cdminfo.dminfo.error + cdminfo.dmequad2) *
    cdminfo.dmefac *
    cdminfo.dmefac

"""
    WidebandTOA

A single wideband TOA observation.

`toa.value` is the wideband TOA measurement in the TDB frame incorporating the clock 
corrections. `toa.ephem` contains the solar system ephemerides. These are computed using 
`PINT`.

References:
    [Pennucci+ 2014](http://doi.org/10.1088/0004-637X/790/2/93),
    [Pennucci 2019](http://doi.org/10.3847/1538-4357/aaf6ef),
    [Luo+ 2021](http://doi.org/10.3847/1538-4357/abe62f)
"""
struct WidebandTOA <: TOABase
    toa::TOA
    dminfo::DMInfo
end

show(io::IO, wtoa::WidebandTOA) = print(
    io,
    "WidebandTOA[MJD:$(trunc(Int, wtoa.toa.value.x/day_to_s)), Freq(MHz):$(trunc(Int, wtoa.toa.observing_frequency.x/1e6))]",
)
show(io::IO, toas::Vector{WidebandTOA}) =
    print(io, "[Vector containing $(length(toas)) wideband TOAs.]")

"""
    CorrectedWidebandTOA

The accumulated timing & noise model corrections applied to a wideband TOA.
"""
struct CorrectedWidebandTOA <: CorrectedTOABase
    corrected_toa::CorrectedTOA
    corrected_dminfo::CorrectedDMInfo
end

CorrectedWidebandTOA(wtoa::WidebandTOA) =
    CorrectedWidebandTOA(CorrectedTOA(wtoa.toa), CorrectedDMInfo(wtoa.dminfo))
