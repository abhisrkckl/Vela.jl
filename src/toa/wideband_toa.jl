"""DM information associated with a wideband TOA."""
struct DMInfo
    value::GQ{-1,Float64}
    error::GQ{-1,Float64}
end

"""The accumulated timing & noise model corrections applied to wideband DM."""
struct CorrectedDMInfo
    dminfo::DMInfo
    model_dm::GQ{-1,Float64}
    dmefac::GQ{0,Float64}
    dmequad2::GQ{-2,Float64}
end

CorrectedDMInfo(dminfo::DMInfo) =
    CorrectedDMInfo(dminfo, GQ{-1}(0.0), dimensionless(1.0), GQ{-2}(0.0))

correct_dminfo(
    cdminfo::CorrectedDMInfo;
    model_dm = GQ{-1}(0.0),
    dmefac = dimensionless(1.0),
    dmequad2 = GQ{-2}(0.0),
) = CorrectedDMInfo(
    cdminfo.dminfo,
    cdminfo.model_dm + model_dm,
    cdminfo.dmefac * dmefac,
    cdminfo.dmequad2 + dmequad2,
)

dm_residual(cdminfo::CorrectedDMInfo) = cdminfo.dminfo.value - cdminfo.model_dm

scaled_dm_error_sqr(cdminfo::CorrectedDMInfo) =
    (cdminfo.dminfo.error * cdminfo.dminfo.error + cdminfo.dmequad2) *
    cdminfo.dmefac *
    cdminfo.dmefac

"""A single wideband TOA observation."""
struct WidebandTOA <: TOABase
    toa::TOA
    dminfo::DMInfo
end

"""The accumulated timing & noise model corrections applied to a wideband TOA."""
struct CorrectedWidebandTOA <: CorrectedTOABase
    corrected_toa::CorrectedTOA
    corrected_dminfo::CorrectedDMInfo
end

CorrectedWidebandTOA(wtoa::WidebandTOA) =
    CorrectedWidebandTOA(CorrectedTOA(wtoa.toa), CorrectedDMInfo(wtoa.dminfo))
