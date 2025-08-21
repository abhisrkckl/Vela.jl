@testset "priors" begin
    @testset "SimplePrior" begin
        params = (PHOFF = dimensionless(0.1), F = (frequency(100.0), GQ{-2}(-1e-14)))
        bad_params = (PHOFF = dimensionless(1.1), F = (frequency(-200.0), GQ{-2}(-2e-14)))

        phoff_prior = SimplePrior{:PHOFF}(Uniform(-0.5, 0.5), Vela.USER_DEFINED_PRIOR)
        @test isfinite(lnprior(phoff_prior, params))
        @test @ballocated(lnprior($phoff_prior, $params)) == 0
        @test !isfinite(lnprior(phoff_prior, bad_params))
        @test prior_transform(phoff_prior, 0.5) == 0
        @test Vela.param_name(phoff_prior) == :PHOFF
        @test Vela.distr_args(phoff_prior) == (-0.5, 0.5)
        @test Vela.distr_name(phoff_prior) == :Uniform

        f0_prior = SimplePriorMulti{:F,UInt(1)}(
            truncated(Normal(100.0, 1e-7); lower = 0.0),
            Vela.USER_DEFINED_PRIOR,
        )
        @test isfinite(lnprior(f0_prior, params))
        @test @ballocated(lnprior($f0_prior, $params)) == 0
        @test !isfinite(lnprior(f0_prior, bad_params))
        @test prior_transform(f0_prior, 0.5) ≈ 100.0
        @test Vela.param_name(f0_prior) == :F
        @test Vela.param_index(f0_prior) == 1
        @test Vela.distr_args(f0_prior) == (100.0, 1e-7)
        @test Vela.distr_name(f0_prior) == :Normal

        f1_prior = SimplePriorMulti{:F,UInt(2)}(
            Uniform(-1.01e-14, -0.9e-14),
            Vela.USER_DEFINED_PRIOR,
        )
        @test isfinite(lnprior(f1_prior, params))
        @test @ballocated(lnprior($f1_prior, $params)) == 0
        @test !isfinite(lnprior(f1_prior, bad_params))
        @test prior_transform(f1_prior, 0.0) == -1.01e-14
        @test Vela.param_name(f1_prior) == :F
        @test Vela.param_index(f1_prior) == 2

        priors = (phoff_prior, f0_prior, f1_prior)
        @test isfinite(lnprior(priors, params))
        @test @ballocated(lnprior($priors, $params)) == 0
        @test !isfinite(lnprior(priors, bad_params))
        @test prior_transform(priors, [0.5, 0.5, 0.0]) ≈ [0.0, 100.0, -1.01e-14]
    end

    @testset "custom priors" begin
        Rmax = 1.66778354e12 # Radius of the Galaxy in s
        dists = [
            KINPriorDistribution(),
            SINIPriorDistribution(),
            STIGMAPriorDistribution(),
            SHAPMAXPriorDistribution(),
            PXPriorDistribution(Rmax),
        ]
        for d in dists
            xl = minimum(d) - 1
            if isfinite(xl)
                @test pdf(d, xl) == 0
                @test cdf(d, xl) == 0
            end

            xh = maximum(d) + 1
            if isfinite(xh)
                @test pdf(d, xh) == 0
                @test cdf(d, xh) == 1
            end

            xin = quantile(d, 0.5)
            @test isfinite(xin)

            @test pdf(d, xin) > 0
            @test 0 < cdf(d, xin) < 1
            @test isfinite(logpdf(d, xin))
            @test isfinite(logcdf(d, xin))

            @test quantile(d, 0.0) == minimum(d)
            @test quantile(d, 1.0) == maximum(d)
        end
    end

    @testset "prior scaling" begin
        args = (1.0, 2.0)
        distr_types = [
            Arcsine,
            Beta,
            BetaPrime,
            Biweight,
            Cauchy,
            Chi,
            Chisq,
            Cosine,
            Epanechnikov,
            Erlang,
            Exponential,
            FDist,
            Frechet,
            Gamma,
            GeneralizedExtremeValue,
            GeneralizedPareto,
            Gumbel,
            InverseGamma,
            InverseGaussian,
            Kolmogorov,
            Kumaraswamy,
            Laplace,
            Levy,
            Logistic,
            LogNormal,
            LogUniform,
            Normal,
            Pareto,
            Rayleigh,
            Semicircle,
            SymTriangularDist,
            TDist,
            TriangularDist,
            Triweight,
            Uniform,
            VonMises,
            Weibull,
        ]

        for D in distr_types
            @test all(map(isfinite, Vela.prior_scaling(D, args)))
        end

        @test all(isfinite.(Vela.prior_scaling(PGeneralizedGaussian, (0.0, 1.5, 0.4))))

        @test Vela.scale_prior_args(Uniform, args, 2.0) == (2.0, 4.0)
        @test Vela.unscale_prior_args(Uniform, args, 1 / 2.0) == (2.0, 4.0)
    end
end
