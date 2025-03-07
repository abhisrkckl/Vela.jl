@testset "ncp" begin
    @testset "matrix_ops" begin
        ntoa = 1000
        npar = 5

        M = randn(ntoa, npar)
        Ndiag = 0.1 .+ rand(ntoa)
        Phidiag = 0.1 .+ rand(npar)
        y = randn(ntoa)
        alpha = randn(npar)

        Sigmainv_brute = Diagonal(1 ./ Phidiag) + transpose(M) * (Diagonal(1 ./ Ndiag) * M)
        MT_Ninv_y_brute = transpose(M) * (Diagonal(1 ./ Ndiag) * y)

        Sigmainv, MT_Ninv_y = Vela.calc_Sigmainv_and_MT_Ninv_y(M, Ndiag, Phidiag, y)

        @test isapprox(Sigmainv, Sigmainv_brute)
        @test isapprox(MT_Ninv_y, MT_Ninv_y_brute)
    end
end
