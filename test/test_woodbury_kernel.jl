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

    @test calc_noise_weights_inv(kernel, params)[1:10] ==
          calc_noise_weights_inv(kernel, params)[11:20]
    @test calc_noise_weights_inv(kernel, params)[21:35] ==
          calc_noise_weights_inv(kernel, params)[36:50]
    @test calc_noise_weights_inv(kernel, params)[51:55] ==
          calc_noise_weights_inv(kernel, params)[56:60]

    @testset "gls likelihood utils" begin
        y = randn(100)
        Ndiag = 1 .+ rand(100)
        M = randn(100, 5)
        Phiinv = 1 .+ rand(5)

        @testset "_calc_y_Ninv_y" begin
            y_Ninv_y = Vela._calc_y_Ninv_y(WhiteNoiseKernel(), Ndiag, y)
            @test y_Ninv_y ≈ dot(y, y ./ Ndiag)
        end

        @testset "_calc_Σinv__and__MT_Ninv_y" begin
            Σinv, MT_Ninv_y =
                Vela._calc_Σinv__and__MT_Ninv_y(WhiteNoiseKernel(), M, Ndiag, Phiinv, y)
            @test MT_Ninv_y ≈ transpose(M) * (y ./ Ndiag)
            @test Σinv ≈ Diagonal(Phiinv) + transpose(M) * (M ./ Ndiag)
        end

        @testset "_gls_lnlike_serial" begin
            lnlike = Vela._gls_lnlike_serial(WhiteNoiseKernel(), M, Ndiag, Phiinv, y)

            Ninv_y = y ./ Ndiag
            y_Ninv_y = dot(y, Ninv_y)
            MT_Ninv_y = transpose(M) * Ninv_y
            MT_Ninv_M = transpose(M) * (M ./ Ndiag)
            Sigmainv = Diagonal(Phiinv) + MT_Ninv_M
            Sigma_MT_Ninv_y = Sigmainv \ MT_Ninv_y
            y_Ninv_M_Sigma_MT_Ninv_y = dot(MT_Ninv_y, Sigma_MT_Ninv_y)
            logdet_N = sum(log.(Ndiag))
            logdet_Phi = -sum(log.(Phiinv))
            logdet_Sigmainv = logdet(Sigmainv)
            lnlike_brute =
                -0.5 * (
                    y_Ninv_y - y_Ninv_M_Sigma_MT_Ninv_y +
                    logdet_N +
                    logdet_Phi +
                    logdet_Sigmainv
                )
            @test lnlike_brute ≈ lnlike

            C = Diagonal(Ndiag) + M * Diagonal(1 ./ Phiinv) * transpose(M)
            @test lnlike ≈ -0.5 * (dot(y, C \ y) + logdet(C))
        end
    end
end
