export Photon

struct Photon <: TOABase
    arrival_time::TOA
    weight::GQ{0,Float64}

    function Photon(arrival_time, frequency, weight, ephem, index)
        toa = TOA(arrival_time, time(0.0), frequency, dimensionless(0.0), ephem, index)
        return new(toa, weight)
    end
end

struct PhotonCorrection <: TOACorrectionBase
    arrival_time_correction::TOACorrection
end

function PhotonCorrection(delay, phase, ssb_psr_pos)
    @assert all(iszero.(ssb_psr_pos)) ||
            dot(ssb_psr_pos, ssb_psr_pos) ≈ dimensionless(1.0) "ssb_psr_pos must be a zero vector (representing pending computation) or a unit vector."

    return PhotonCorrection(TOACorrection(delay, phase, dimensionless(1.0), GQ{2}(0.0), frequency(0.0), dimensionless(1.1), ssb_psr_pos))
end

PhotonCorrection() = PhotonCorrection(time(0.0), dimensionless(0.0), dimensionless.((0.0, 0.0, 0.0)),)

fractional_phase(photcorr::PhotonCorrection) = photcorr.arrival_time_correction.phase - round(photcorr.arrival_time_correction.phase)

