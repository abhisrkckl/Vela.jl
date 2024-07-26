@testset "priors" begin
    @testset "SimplePrior" begin
        params = (PHOFF = dimensionless(0.1), F = (frequency(100.0), GQ(-1e-14, -2)))
        bad_params = (PHOFF = dimensionless(1.1), F = (frequency(-200.0), GQ(-2e-14, -2)))

        phoff_prior = SimplePrior{:PHOFF}(Uniform(-0.5, 0.5))
        @test isfinite(lnprior(phoff_prior, params))
        @test @ballocated(lnprior($phoff_prior, $params)) == 0
        @test !isfinite(lnprior(phoff_prior, bad_params))
        @test prior_transform(phoff_prior, 0.5) == 0

        f0_prior = SimplePriorMulti{:F,UInt(1)}(truncated(Normal(100.0, 1e-7); lower = 0.0))
        @test isfinite(lnprior(f0_prior, params))
        @test @ballocated(lnprior($f0_prior, $params)) == 0
        @test !isfinite(lnprior(f0_prior, bad_params))
        @test prior_transform(f0_prior, 0.5) ≈ 100.0

        f1_prior = SimplePriorMulti{:F,UInt(2)}(Uniform(-1.01e-14, -0.9e-14))
        @test isfinite(lnprior(f1_prior, params))
        @test @ballocated(lnprior($f1_prior, $params)) == 0
        @test !isfinite(lnprior(f1_prior, bad_params))
        @test prior_transform(f1_prior, 0.0) == -1.01e-14

        priors = (phoff_prior, f0_prior, f1_prior)
        @test isfinite(lnprior(priors, params))
        @test @ballocated(lnprior($priors, $params)) == 0
        @test !isfinite(lnprior(priors, bad_params))
        @test prior_transform(priors, [0.5, 0.5, 0.0]) ≈ [0.0, 100.0, -1.01e-14]
    end
end
