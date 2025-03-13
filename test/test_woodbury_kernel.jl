@testset "WoodburyKernel" begin
    kernel = WoodburyKernel(
        WhiteNoiseKernel(),
        (
            PowerlawRedNoiseGP(10),
            PowerlawDispersionNoiseGP(15),
            PowerlawChromaticNoiseGP(5),
        ),
        randn(1000, 60),
    )

    params = (
        TNREDAMP = dimensionless(-13.5),
        TNREDGAM = dimensionless(3.5),
        PLREDFREQ = frequency(1e-8),
        TNDMAMP = dimensionless(-13.5),
        TNDMGAM = dimensionless(3.5),
        PLDMFREQ = frequency(1e-8),
        TNCHROMIDX = dimensionless(4.0),
        TNCHROMAMP = dimensionless(-13.5),
        TNCHROMGAM = dimensionless(3.5),
        PLCHROMFREQ = frequency(1e-8),
    )

    @test length(calc_noise_weights_inv(kernel, params)) == 60
    @test all(calc_noise_weights_inv(kernel, params) .> 0)

    @testset "gls likelihood utils" begin
        y = randn(100)
        Ndiag = 1 .+ rand(100)
        M = randn(100, 5)
        Phiinv = 1 .+ rand(5)
        @testset "_calc_y_Ninv_y" begin
            y_Ninv_y = Vela._calc_y_Ninv_y(Ndiag, y)
            @test y_Ninv_y ≈ dot(y, y ./ Ndiag)
        end

        @testset "_calc_Σinv__and__MT_Ninv_y" begin
            Σinv, MT_Ninv_y = Vela._calc_Σinv__and__MT_Ninv_y(M, Ndiag, Phiinv, y)
            @test MT_Ninv_y ≈ transpose(M) * (y ./ Ndiag)
            @test Σinv ≈ Diagonal(Phiinv) + transpose(M) * (M ./ Ndiag)
        end
    end
end
