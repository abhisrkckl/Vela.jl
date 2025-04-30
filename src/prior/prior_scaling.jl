export scale_prior_args

prior_scaling(::Type{Arcsine}, args) = map(oneunit, args)

prior_scaling(::Type{Beta}, args) = map(zero, args)

prior_scaling(::Type{BetaPrime}, args) = map(zero, args)

prior_scaling(::Type{Biweight}, args) = map(oneunit, args)

prior_scaling(::Type{Cauchy}, args) = map(oneunit, args)

prior_scaling(::Type{Chi}, args) = map(zero, args)

prior_scaling(::Type{Chisq}, args) = map(zero, args)

prior_scaling(::Type{Cosine}, args) = map(oneunit, args)

prior_scaling(::Type{Epanechnikov}, args) = map(oneunit, args)

prior_scaling(::Type{Erlang}, args) = ((), (0,), (0, 1))[length(args)+1]

prior_scaling(::Type{Exponential}, args) = map(oneunit, args)

prior_scaling(::Type{FDist}, args) = map(zero, args)

prior_scaling(::Type{Frechet}, args) = length(args) == ((), (0,), (0, 1))[length(args)+1]

prior_scaling(::Type{Gamma}, args) = length(args) == ((), (0,), (0, 1))[length(args)+1]

prior_scaling(::Type{GeneralizedExtremeValue}, args) = (1, 1, 0)

prior_scaling(::Type{GeneralizedPareto}, args) =
    ((), (0,), (1, 0), (1, 1, 0))[length(args)+1]

prior_scaling(::Type{Gumbel}, args) = map(oneunit, args)

prior_scaling(::Type{InverseGamma}, args) = ((), (0,), (0, 1))[length(args)+1]

prior_scaling(::Type{InverseGaussian}, args) = map(oneunit, args)

prior_scaling(::Type{Kolmogorov}, args) = ()

prior_scaling(::Type{Kumaraswamy}, args) = (0, 0)

prior_scaling(::Type{Laplace}, args) = map(oneunit, args)

prior_scaling(::Type{Levy}, args) = map(oneunit, args)

prior_scaling(::Type{Logistic}, args) = map(oneunit, args)

prior_scaling(::Type{LogNormal}, args) = map(zero, args)

prior_scaling(::Type{LogUniform}, args) = map(oneunit, args)

prior_scaling(::Type{Normal}, args) = map(oneunit, args)

prior_scaling(::Type{Pareto}, args) = ((), (0,), (0, 1))[length(args)+1]

prior_scaling(::Type{Rayleigh}, args) = map(oneunit, args)

prior_scaling(::Type{Semicircle}, args) = map(oneunit, args)

prior_scaling(::Type{SymTriangularDist}, args) = map(oneunit, args)

prior_scaling(::Type{TDist}, args) = map(zero, args)

prior_scaling(::Type{TriangularDist}, args) = map(oneunit, args)

prior_scaling(::Type{Triweight}, args) = map(oneunit, args)

prior_scaling(::Type{Uniform}, args) = map(oneunit, args)

prior_scaling(::Type{VonMises}, args) = ((), (0,), (0, 1))[length(args)+1]

prior_scaling(::Type{Weibull}, args) = ((), (0,), (0, 1))[length(args)+1]

scale_prior_args(distr_type, args, unit_conversion_factor) =
    args .* (unit_conversion_factor .^ prior_scaling(distr_type, args))

unscale_prior_args(distr_type, args, unit_conversion_factor) =
    args .* ((1 / unit_conversion_factor) .^ prior_scaling(distr_type, args))
