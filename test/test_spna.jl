@testset "SPNA" begin
    @testset "Resids" begin
        nbres = NarrowbandResid(-1.1e-6, 1e-6)
        cnbres = NarrowbandResidCorrection()
        @test corrected_time_residual(nbres, cnbres) ≈ nbres.tres
        @test scaled_toa_error_sqr(nbres, cnbres) ≈ nbres.terr2

        wbres = WidebandResid(-1.1e-6, 1e+4, 1e-6, 1e+5)
        cwbres = WidebandResidCorrection()
        @test corrected_time_residual(wbres, cwbres) ≈ wbres.tres
        @test scaled_toa_error_sqr(wbres, cwbres) ≈ wbres.terr2
        @test corrected_dm_residual(wbres, cwbres) ≈ wbres.dres
        @test scaled_dm_error_sqr(wbres, cwbres) ≈ wbres.derr2
    end
end
