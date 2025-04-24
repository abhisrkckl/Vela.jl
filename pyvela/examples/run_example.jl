#!/usr/bin/env julia

# This script demonstrates the use of `Vela` in Julia.
# It needs several packages that are not dependencies of `Vela`.
# Please install them separately.

using Vela
using AffineInvariantMCMC
using PairPlots
using CairoMakie
using DataFrames
using CSV

function main()
    parfile, timfile = ARGS

    # Vela doesn't know how to read par and tim files. 
    # They must be converted into a JLSO file first. 
    run(`pyvela-jlso $parfile $timfile -o __$parfile.jlso`)

    m, t = Vela.load_pulsar_data("__$parfile.jlso")

    lnpost = get_lnpost_func(m, t)
    prior_transform = get_prior_transform_func(m)

    param_names = get_free_param_names(m)

    ndim = length(param_names)
    thin = 10
    nsamples_perwalker = 5000
    nwalker = ndim * 5
    burnin = 1500

    x0 = mapslices(prior_transform, rand(ndim, nwalker); dims = 1)

    print("Burning in...")
    burnin_chain, _ = @time AffineInvariantMCMC.sample(lnpost, nwalker, x0, burnin, 1)
    x1 = burnin_chain[:, :, end]

    print("Running actual MCMC...")
    chain, lnpost_vals =
        @time AffineInvariantMCMC.sample(lnpost, nwalker, x1, nsamples_perwalker, thin)
    flat_chain, flat_lnpost_vals = AffineInvariantMCMC.flattenmcmcarray(chain, lnpost_vals)

    print("Plotting...")
    scale_factors = get_scale_factors(m)
    chain_df = DataFrame((flat_chain ./ scale_factors)', param_names)
    fig = pairplot(chain_df)
    save("$(m.pulsar_name)_chain_emcee_jl.png", fig)

    print("Saving chain...")
    CSV.write("$(m.pulsar_name)_chain_emcee_jl.txt", chain_df; delim = "    ")
end

main()
