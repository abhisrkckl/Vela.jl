export FrequencyDependent, FrequencyDependentJump

"""A frequency-dependent delay to account for frequency-dependent profile evolution.

Reference:
    [Arzoumanian+ 2015](http://doi.org/10.1088/0004-637X/813/1/65)
"""
struct FrequencyDependent <: DelayComponent end

function delay(
    ::FrequencyDependent,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
)::GQ
    fds = params.FD
    ν = doppler_corrected_observing_frequency(toa, toacorr)
    νref = frequency(1e9) # 1 GHz
    λ = log(ν / νref)
    return sum(fd * λ^p for (p, fd) in enumerate(fds))
end

"""A frequency-dependent delay to account for frequency-dependent profile evolution.

Reference:
    [Susobhanan+ 2015](http://doi.org/10.3847/1538-4357/ad59f7)
"""
struct FrequencyDependentJump <: DelayComponent
    jump_mask::BitMatrix
    exponents::Vector{UInt}

    function FrequencyDependentJump(mask, exps)
        @assert size(mask)[1] == length(exps)
        return new(mask, exps)
    end
end

function delay(
    fdj::FrequencyDependentJump,
    toa::TOA,
    toacorr::TOACorrection,
    params::NamedTuple,
)::GQ
    if is_tzr(toa)
        return time(0.0)
    end

    fdjumps = params.FDJUMP
    ν = doppler_corrected_observing_frequency(toa, toacorr)
    νref = frequency(1e9)

    mask = @view(fdj.jump_mask[:, toa.index])

    λ = log(ν / νref)

    delay = time(0.0)
    for (fdjump, p, m) in zip(fdjumps, fdj.exponents, mask)
        if m
            delay += fdjump * λ^p
        end
    end

    return delay
end

function show(io::IO, fdjmp::FrequencyDependentJump)
    num_jumps = length(fdjmp.exponents)
    print(io, "FrequencyDependentJump($num_jumps FDJUMPs)")
end
