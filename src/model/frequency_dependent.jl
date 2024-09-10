export FrequencyDependent

"""A frequency-dependent delay to account for frequency-dependent profile evolution.

Reference:
    [Arzoumanian+ 2015](http://doi.org/10.1088/0004-637X/813/1/65)
"""
struct FrequencyDependent <: DelayComponent end

function delay(::FrequencyDependent, ctoa::CorrectedTOA, params::NamedTuple)::GQ
    fds = params.FD
    ν = doppler_corrected_observing_frequency(ctoa)
    νref = frequency(1e9) # 1 GHz
    return sum(fd * (log(ν / νref))^p for (p, fd) in enumerate(fds))
end
