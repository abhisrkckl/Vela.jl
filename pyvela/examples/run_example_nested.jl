using Vela
using NestedSamplers
using PairPlots
using CairoMakie
using DataFrames
using CSV
using StatsBase

function main()
    jlsofile = ARGS[0]
    m, t = Vela.load_pulsar_data(jlsofile)

    lnlike = get_lnlike_func(m, t)
    prior_transform = get_prior_transform_func(m)

    param_names = get_free_param_names(m)
    ndim = length(param_names)

    nested_model = NestedModel(lnlike, prior_transform)

    bounds = Bounds.MultiEllipsoid
    # prop = Proposals.Slice(slices=10)
    prop = Proposals.Uniform()
    nlive = 1000
    sampler = Nested(ndim, nlive; bounds = bounds, proposal = prop)

    chain, state = sample(nested_model, sampler; dlogz = 0.1)

    # Resample with equal weights
    chain_resampled = sample(chain, Weights(vec(chain["weights"])), length(chain))

    scale_factors = get_scale_factors(m)
    chain_scaled = (chain_resampled.value.data[:, 1:(end-1), 1]' ./ scale_factors)'
    chain_df = DataFrame(chain_scaled, param_names)

    print("logZ = $(state.logz) ± $(state.logzerr)")

    fig = pairplot(chain_df)
    save("$(m.pulsar_name)_chain_nested_jl.png", fig)

    print("Saving chain...")
    CSV.write("$(m.pulsar_name)_chain_nested_jl.txt", chain_df; delim = "    ")
end

main()
