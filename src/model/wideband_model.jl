"""Update a `CorrectedWidebandTOA` object. 
Most components only affect the TOA part and not the DM part."""
correct_toa(
    component::Component,
    wtoa::WidebandTOA,
    wtoacorr::WidebandTOACorrection,
    params::NamedTuple,
) = WidebandTOACorrection(
    correct_toa(component, wtoa.toa, wtoacorr.toa_correction, params),
    wtoacorr.dm_correction,
)

"""Update a `CorrectedWidebandTOA` object. 
`DispersionComponent`s affect both the TOA part and the DM part."""
function correct_toa(
    dmcomp::DispersionComponent,
    wtoa::WidebandTOA,
    wtoacorr::WidebandTOACorrection,
    params::NamedTuple,
)
    dm = dispersion_slope(dmcomp, wtoa.toa, wtoacorr.toa_correction, params)
    delay =
        dm / doppler_corrected_observing_frequency(wtoa.toa, wtoacorr.toa_correction)^Val(2)
    toacorr = correct_toa(wtoacorr.toa_correction; delay = delay)
    dmcorr = correct_dminfo(wtoacorr.dm_correction; delta_dm = dm)
    return WidebandTOACorrection(toacorr, dmcorr)
end

"""Update a `WidebandTOA` object using a timing model."""
correct_toa(model::TimingModel, wtoa::WidebandTOA, params::NamedTuple) =
    correct_toa(model, wtoa, WidebandTOACorrection(), params)
