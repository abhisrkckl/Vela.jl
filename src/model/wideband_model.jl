"""Update a `CorrectedWidebandTOA` object. 
Most components only affect the TOA part and not the DM part."""
correct_toa(component::Component, cwtoa::CorrectedWidebandTOA, params::NamedTuple) =
    CorrectedWidebandTOA(
        correct_toa(component, cwtoa.corrected_toa, params),
        cwtoa.corrected_dminfo,
    )

"""Update a `CorrectedWidebandTOA` object. 
`DispersionComponent`s affect both the TOA part and the DM part."""
function correct_toa(
    dmcomp::DispersionComponent,
    cwtoa::CorrectedWidebandTOA,
    params::NamedTuple,
)
    dm = dispersion_slope(dmcomp, cwtoa.corrected_toa, params)
    delay = dm / doppler_corrected_observing_frequency(cwtoa.corrected_toa)^Val(2)
    ctoa = correct_toa(cwtoa.corrected_toa; delay = delay)
    cdminfo = correct_dminfo(cwtoa.corrected_dminfo; delta_dm = dm)
    return CorrectedWidebandTOA(ctoa, cdminfo)
end

"""Update a `WidebandTOA` object using a timing model."""
correct_toa(model::TimingModel, toa::WidebandTOA, params::NamedTuple) =
    correct_toa(model, CorrectedWidebandTOA(toa), params)
