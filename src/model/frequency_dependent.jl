export FrequencyDependent

"""A frequency-dependent delay to account for frequency-dependent profile evolution."""
struct FrequencyDependent <: DelayComponent end

function delay(::FrequencyDependent, ctoa::CorrectedTOA, params::NamedTuple)::GQ
    fds = params.FD
    ν = doppler_corrected_observing_frequency(ctoa)
    νref = frequency(1e9) # 1 GHz
    return sum(fd * (log(ν / νref))^p for (p, fd) in enumerate(fds))
end
