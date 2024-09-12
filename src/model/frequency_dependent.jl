export FrequencyDependent, FrequencyDependentJump

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

"""A frequency-dependent delay to account for frequency-dependent profile evolution.

Reference:
    [Susobhanan+ 2015](http://doi.org/10.3847/1538-4357/ad59f7)
"""
struct FrequencyDependentJump <: DelayComponent
    jump_mask::BitMatrix
    exponents::Vector{UInt}

    function FrequencyDependentJump(mask, idxs)
        @assert size(mask)[1] == length(idxs)
        return new(jump_mask, idxs)
    end
end

function delay(fdj::FrequencyDependentJump, ctoa::CorrectedTOA, params::NamedTuple)::GQ
    if ctoa.toa.tzr
        return time(0.0)
    end

    fdjumps = params.FDJUMP
    ν = doppler_corrected_observing_frequency(ctoa)
    νref = frequency(1e9)

    mask = @view(fdj.jump_mask[:, ctoa.toa.index])

    delay = time(0.0)
    for (fdjump, p, m) in zip(fdjumps, fdj.exponents, mask)
        delay += m * fdjump * (log(ν / νref))^p
    end

    return delay
end