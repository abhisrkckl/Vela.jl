export DispersionTaylor, dispersion_slope

struct DispersionTaylor <: DispersionComponent end

read_params_from_dict(::DispersionTaylor, params::Dict) =
    (DMEPOCH = params[:DMEPOCH][1], DM = params[:DM])

function dispersion_slope(::DispersionTaylor, toa::TOA, params)
    t0 = params.DMEPOCH
    t = toa.value
    dms = params.DM
    dm = taylor_horner(t - t0, dms)
    return quantity_like(dms[1], dm.x)
end
