using Vela
using AffineInvariantMCMC
using PythonCall
using PairPlots
using CairoMakie
using DataFrames
using CSV

p2v = pyimport("pint2vela")


parfile, timfile = ARGS
const m, t =
    map(pyconvert, (TimingModel, Vector{TOA}), p2v.read_model_and_toas(parfile, timfile))

const lnpost = get_lnpost_func(m, t)
const prior_transform = get_prior_transform_func(m)

const param_names = get_free_param_names(m)

const ndim = length(param_names)
const nwalker = 100
const thin = 10
const nsamples_perwalker = 1000
const burnin = 100

# Initial points
const x0 = mapslices(prior_transform, rand(ndim, nwalker); dims = 1)

# Burn-in
print("Burning in...")
const burnin_chain, _ = @time AffineInvariantMCMC.sample(lnpost, nwalker, x0, burnin, 1)
const x1 = burnin_chain[:, :, end]

# Do the actual MCMC run
print("Running actual MCMC...")
const chain, lnpost_vals =
    @time AffineInvariantMCMC.sample(lnpost, nwalker, x1, nsamples_perwalker, thin)
const flat_chain, flat_lnpost_vals =
    AffineInvariantMCMC.flattenmcmcarray(chain, lnpost_vals)

print("Plotting...")
const scale_factors = get_scale_factors(m)
const chain_df = DataFrame((flat_chain ./ scale_factors)', param_names)
fig = pairplot(chain_df)
save("$(m.pulsar_name)_chain_emcee_jl.png", fig)

print("Saving chain...")
CSV.write("$(m.pulsar_name)_chain_emcee_jl.txt", chain_df; delim = "    ")
