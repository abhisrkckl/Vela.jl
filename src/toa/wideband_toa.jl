export DMInfo,
    DMInfoCorrection,
    correct_dminfo,
    dm_residual,
    scaled_dm_error_sqr,
    WidebandTOA,
    WidebandTOACorrection

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

"""The accumulated timing & noise model corrections applied to wideband DM."""
struct DMInfoCorrection
    model_dm::GQ{-1,Float64}
    dmefac::GQ{0,Float64}
    dmequad2::GQ{-2,Float64}
end

DMInfoCorrection() = DMInfoCorrection(GQ{-1}(0.0), dimensionless(1.0), GQ{-2}(0.0))

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
    dmcorr::DMInfoCorrection;
    delta_dm = GQ{-1}(0.0),
    dmefac = dimensionless(1.0),
    dmequad2 = GQ{-2}(0.0),
) = DMInfoCorrection(
    dmcorr.model_dm + delta_dm,
    dmcorr.dmefac * dmefac,
    dmcorr.dmequad2 + dmequad2,
)

dm_residual(dminfo::DMInfo, dmcorr::DMInfoCorrection) = dminfo.value - dmcorr.model_dm

scaled_dm_error_sqr(dminfo::DMInfo, dmcorr::DMInfoCorrection) =
    (dminfo.error * dminfo.error + dmcorr.dmequad2) * dmcorr.dmefac * dmcorr.dmefac

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

"""The accumulated timing & noise model corrections applied to a wideband TOA."""
struct WidebandTOACorrection <: TOACorrectionBase
    toa_correction::TOACorrection
    dm_correction::DMInfoCorrection
end

WidebandTOACorrection() = WidebandTOACorrection(TOACorrection(), DMInfoCorrection())
