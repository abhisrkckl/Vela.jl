export get_lnpost_func

function get_lnpost_func(model::TimingModel, toas::Vector{TOA})
    lnlike_func = get_lnlike_func(model, toas)
    lnprior_func = get_lnprior_func(model)

    function lnpost(params)
        lnpr = lnprior_func(params)
        return isfinite(lnpr) ? lnpr + lnlike_func(params) : -Inf
    end

    return lnpost
end
