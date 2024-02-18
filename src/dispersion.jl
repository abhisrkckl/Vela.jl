export DispersionTaylor, dispersion_slope

struct DispersionTaylor <: DispersionComponent
    number_of_terms::UInt
end

function dispersion_slope(dmt::DispersionTaylor, toa::TOA, params)
    if dmt.number_of_terms == 1
        return params["DM"]
    end

    t0 = params["DMEPOCH"]
    t = toa.value
    cs = [params["DM$i"] for i = 1:(dmt.number_of_terms-1)]
    c0 = params["DM"]
    th = TaylorSeries(t0, c0, cs)
    return th(t)
end
