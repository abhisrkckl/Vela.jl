export ResidBase,
    ResidCorrectionBase,
    NarrowbandResid,
    NarrowbandResidCorrection,
    WidebandResid,
    WidebandResidCorrection,
    corrected_time_residual,
    corrected_dm_residual

abstract type ResidBase end
abstract type ResidCorrectionBase end

struct NarrowbandResid <: ResidBase
    tres::Float64
    terr2::Float64
end

function calc_resids(model::TimingModel, toas::Vector{TOA})
    tresids =
        map(value, form_residuals(model, toas, model.param_handler._default_params_tuple))
    terrs2 = map(toa -> value(toa.error)^2, toas)
    return map(NarrowbandResid, tresids, terrs2)
end

struct NarrowbandResidCorrection <: ResidCorrectionBase
    delay::Float64
    efac::Float64
    equad2::Float64
end

NarrowbandResidCorrection() = NarrowbandResidCorrection(0.0, 1.0, 0.0)

corrected_time_residual(res::ResidBase, cres::ResidCorrectionBase) = res.tres - cres.delay

scaled_toa_error_sqr(res::ResidBase, cres::ResidCorrectionBase) =
    cres.efac * cres.efac * (res.terr2 + cres.equad2)

struct WidebandResid <: ResidBase
    tres::Float64
    dres::Float64
    terr2::Float64
    derr2::Float64
end

function calc_resids(model::TimingModel, wtoas::Vector{WidebandTOA})
    wres = form_residuals(model, wtoas, model.param_handler._default_params_tuple)

    tresids = map(wr -> value(wr[1]), wres)
    dresids = map(wr -> value(wr[2]), wres)

    terrs2 = map(toa -> value(wtoa.toa.error)^2, wtoas)
    derrs2 = map(toa -> value(wtoa.dminfo.error)^2, wtoas)
    return map(WidebandResid, tresids, dresids, terrs2, derrs2)
end

struct WidebandResidCorrection <: ResidCorrectionBase
    delay::Float64
    deltadm::Float64
    efac::Float64
    equad2::Float64
    dmefac::Float64
    dmequad2::Float64
end

WidebandResidCorrection() = WidebandResidCorrection(0.0, 0.0, 1.0, 0.0, 1.0, 0.0)

corrected_dm_residual(res::WidebandResid, cres::WidebandResidCorrection) =
    res.dres - cres.deltadm

scaled_dm_error_sqr(res::WidebandResid, cres::WidebandResidCorrection) =
    cres.dmefac * cres.dmefac * (res.derr2 + cres.dmequad2)
