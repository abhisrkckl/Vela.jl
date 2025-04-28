export SolarWind, SolarWindDispersion

"""
    sun_angle_and_distance

Compute the pulsar-Sun-observatory angle (ρ) and the observatory-sun distance.
This is different from the sun angle defined in PINT; that angle is π-ρ in our
notation. This function assumes that the ssb_psr_pos has already been computed and
is available in the CorrectedTOA object."""
function sun_angle_and_distance(toa::TOA, toacorr::TOACorrection)
    Lhat = toacorr.ssb_psr_pos
    @assert !all(iszero, Lhat) "ssb_psr_pos is not available!"
    rvec = toa.ephem.obs_sun_pos
    r = sqrt(dot(rvec, rvec))
    cos_ρ = -dot(Lhat, rvec) / r
    return acos(cos_ρ), r
end

"""
    SolarWind

Abstract base type for solar wind models."""
abstract type SolarWind <: DispersionComponent end

"""
    SolarWindDispersion

Solar wind model assuming a spherically symmetric distribution of electrons which
falls off as an inverse-square function of the distance from the Sun. Corresponds
to the `SolarWindDispersion` model in `PINT`.

Reference:
    [Edwards+ 2006](http://doi.org/10.1111/j.1365-2966.2006.10870.x)
"""
struct SolarWindDispersion <: SolarWind end

function dispersion_slope(
    ::SolarWindDispersion,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
)::GQ
    ρ, r = sun_angle_and_distance(toa, toacorr)
    t0 = params.SWEPOCH
    t = corrected_toa_value(toa, toacorr, Float64)
    ne_sw = taylor_horner(t - t0, params.NE_SW)
    return ne_sw * AU * AU * ρ / (r * sin(ρ))
end
