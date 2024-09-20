export Spindown, phase, spin_frequency, read_params_from_dict

"""Rotation or the pulsar represented as a Taylor series in spin frequency.

Corresponds to `Spindown` in `PINT`.

Reference:
    [Backer & Hellings 1986](http://doi.org/10.1146/annurev.aa.24.090186.002541)
"""
struct Spindown <: PhaseComponent end

"""Rotational phase of the pulsar."""
spindown_phase(Δt, f_, fs)::GQ{0,Double64} = f_ * Δt + taylor_horner_integral(Δt, fs)

"""Instantaneous rotational frequency of the pulsar."""
spindown_frequency(Δt, f_, fs)::GQ{-1,Float64} = f_ + taylor_horner(Δt, fs)

"""Update the `CorrectedTOA` object with the rotational phase and spin frequency"""
function correct_toa(::Spindown, toa::TOA, toacorr::TOACorrection, params::NamedTuple)
    # t0 = params.PEPOCH
    # t = corrected_toa_value(toa, toacorr)
    # Δt = t - t0
    Δt = corrected_toa_value(toa, toacorr)
    fs = params.F
    f_ = params.F_
    return correct_toa_phase(
        toacorr;
        phase = spindown_phase(Δt, f_, fs),
        delta_spin_frequency = spindown_frequency(GQ{Float64}(Δt), f_, fs),
    )
end
