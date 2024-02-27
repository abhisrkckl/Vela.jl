export DispersionTaylor, dispersion_slope

struct DispersionTaylor <: DispersionComponent end

function dispersion_slope(::DispersionTaylor, toa::TOA, params)
    t0 = params["DMEPOCH"][1]
    t = toa.value
    dms = params["DM"]
    dm = taylor_horner(t - t0, dms)
    return quantity_like(dms[1], dm.x)
end
