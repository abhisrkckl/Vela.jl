export SolarWind, SolarWindDispersion

"""Compute the pulsar-Sun-observatory angle (ρ) and the observatory-sun distance.
This is different from the sun angle defined in PINT; that angle is π-ρ in our
notation. This function assumes that the ssb_psr_pos has already been computed and
is available in the CorrectedTOA object."""
function sun_angle_and_distance(ctoa::CorrectedTOA)
    Lhat = ctoa.ssb_psr_pos
    @assert !all(iszero.(Lhat)) "ssb_psr_pos is not available!"
    rvec = ctoa.toa.ephem.obs_sun_pos
    r = sqrt(dot(rvec, rvec))
    cos_ρ = -dot(Lhat, rvec) / r
    return acos(cos_ρ), r
end

"""Abstract base type for solar wind models."""
abstract type SolarWind <: DispersionComponent end

"""Solar wind model assuming a spherically symmetric distribution of electrons which
falls off as an inverse-square function of the distance from the Sun."""
struct SolarWindDispersion <: SolarWind end

function dispersion_slope(::SolarWindDispersion, ctoa::CorrectedTOA, params)
    ρ, r = sun_angle_and_distance(ctoa)
    return params.NE_SW * AU * AU * ρ / (r * sin(ρ))
end
