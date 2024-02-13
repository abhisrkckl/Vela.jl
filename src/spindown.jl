using GeometricUnits
using Quadmath

export Spindown, phase, rotational_frequency

struct Spindown <: PhaseComponent
    number_of_terms::UInt
end

_to_f128(q::GQ) = GQ(Float128(q.x), q.d)

function phase(spindown::Spindown, toa::TOA, params)::GQ{Float128}
    t0 = _to_f128(params["PEPOCH"])
    t = toa.value
    cs = [_to_f128(params["F$i"]) for i = 0:(spindown.number_of_terms-1)]
    c0 = dimensionless(Float128(0.0))
    th = TaylorSeries(t0, c0, cs)
    return th(t)
end

function spin_frequency(spindown::Spindown, toa::TOA, params)
    if spindown.number_of_terms == 1
        return params["F0"]
    end

    t0 = params["PEPOCH"]
    t = toa.value
    f0 = params["F0"]
    cs = [params["F$i"] for i = 1:(spindown.number_of_terms-1)]
    th = TaylorSeries(t0, f0, cs)
    return th(t)
end

correct_toa(spindown::Spindown, toa::TOA, params) = TOA(
    toa.value,
    toa.error,
    toa.observing_frequency,
    toa.phase + phase(spindown, toa, params),
    spin_frequency(spindown, toa, params),
    toa.ephem_vectors,
    toa.level + 1,
    toa.tzr,
)
