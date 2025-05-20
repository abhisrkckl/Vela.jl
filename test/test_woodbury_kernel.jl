@testset "WoodburyKernel" begin
    kernel = WoodburyKernel(
        WhiteNoiseKernel(),
        (
            PowerlawRedNoiseGP(10, 4, 2.0),
            PowerlawDispersionNoiseGP(15, 4, 2.0),
            PowerlawChromaticNoiseGP(5, 4, 2.0),
        ),
        randn(1000, 84),
    )

    @test length(get_marginalized_param_names(kernel)) == 84

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

    @test length(calc_noise_weights_inv(kernel, params)) == 60 + 24
    @test all(calc_noise_weights_inv(kernel, params) .> 0)

    @test calc_noise_weights_inv(kernel, params)[1:14] ==
          calc_noise_weights_inv(kernel, params)[15:28]
    @test calc_noise_weights_inv(kernel, params)[29:47] ==
          calc_noise_weights_inv(kernel, params)[48:66]
    @test calc_noise_weights_inv(kernel, params)[67:75] ==
          calc_noise_weights_inv(kernel, params)[76:84]

    @testset "gls likelihood utils" begin
        y = randn(100)
        Ndiag = 1 .+ rand(100)
        M = randn(100, 5)
        Phiinv = 1 .+ rand(5)

        @testset "_calc_y_Ninv_y__and__logdet_N" begin
            y_Ninv_y, logdet_N =
                Vela._calc_y_Ninv_y__and__logdet_N(WhiteNoiseKernel(), Ndiag, y, (;))
            @test y_Ninv_y ≈ dot(y, y ./ Ndiag)
            @test logdet_N ≈ sum(log.(Ndiag))
        end

        @testset "_calc_Σinv__and__MT_Ninv_y" begin
            Σinv, MT_Ninv_y = Vela._calc_Σinv__and__MT_Ninv_y(
                WhiteNoiseKernel(),
                M,
                Ndiag,
                Phiinv,
                y,
                (;),
            )
            @test MT_Ninv_y ≈ transpose(M) * (y ./ Ndiag)
            @test Σinv ≈ Diagonal(Phiinv) + transpose(M) * (M ./ Ndiag)
        end

        @testset "_gls_lnlike_serial" begin
            lnlike = Vela._gls_lnlike_serial(WhiteNoiseKernel(), M, Ndiag, Phiinv, y, (;))

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

    @testset "gls ecorr likelihood utils" begin
        y = randn(9)
        Ndiag = 1 .+ rand(9)
        M = randn(9, 5)
        Phiinv = 1 .+ rand(5)

        inner_kernel =
            EcorrKernel([EcorrGroup(1, 3, 1), EcorrGroup(4, 6, 2), EcorrGroup(7, 9, 3)])
        params = (; ECORR = (time(0.5), time(1.1), time(0.9)))

        U = [1 1 1 0 0 0 0 0 0; 0 0 0 1 1 1 0 0 0; 0 0 0 0 0 0 1 1 1]
        Psi = collect(value.(params.ECORR) .^ 2)
        Nc = Diagonal(Ndiag) + transpose(U) * Diagonal(Psi) * U

        @testset "_calc_y_Ninv_y__and__logdet_N" begin
            y_Ninv_y, logdet_N =
                Vela._calc_y_Ninv_y__and__logdet_N(inner_kernel, Ndiag, y, params)
            @test y_Ninv_y ≈ dot(y, Nc \ y)
            @test logdet_N ≈ logdet(Nc)
        end

        @testset "_calc_Σinv__and__MT_Ninv_y" begin
            Σinv, MT_Ninv_y =
                Vela._calc_Σinv__and__MT_Ninv_y(inner_kernel, M, Ndiag, Phiinv, y, params)
            @test MT_Ninv_y ≈ transpose(M) * (Nc \ y)
            @test Σinv ≈ Diagonal(Phiinv) + transpose(M) * (Nc \ M)
        end

        @testset "_gls_lnlike_serial" begin
            lnlike = Vela._gls_lnlike_serial(inner_kernel, M, Ndiag, Phiinv, y, params)

            Ninv_y = Nc \ y
            y_Ninv_y = dot(y, Ninv_y)
            MT_Ninv_y = transpose(M) * Ninv_y
            MT_Ninv_M = transpose(M) * (Nc \ M)
            Sigmainv = Diagonal(Phiinv) + MT_Ninv_M
            Sigma_MT_Ninv_y = Sigmainv \ MT_Ninv_y
            y_Ninv_M_Sigma_MT_Ninv_y = dot(MT_Ninv_y, Sigma_MT_Ninv_y)
            logdet_N = logdet(Nc)
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

            C = Nc + M * Diagonal(1 ./ Phiinv) * transpose(M)
            @test lnlike ≈ -0.5 * (dot(y, C \ y) + logdet(C))
        end
    end
end
