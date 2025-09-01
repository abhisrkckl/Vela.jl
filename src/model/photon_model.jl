"""Update a `CorrectedWidebandTOA` object. 
Most components only affect the TOA part and not the DM part."""
correct_toa(
    component::Component,
    photon::Photon,
    photcorr::PhotonCorrection,
    params::NamedTuple,
) = PhotonCorrection(
    correct_toa(component, photon.toa, photcorr.arrival_time_correction, params),
)

"""Update a `WidebandTOA` object using a timing model."""
correct_toa(model::TimingModel, wtoa::Photon, params::NamedTuple) =
    correct_toa(model, wtoa, PhotonCorrection(), params)
