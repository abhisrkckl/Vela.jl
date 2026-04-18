@testset "ChromaticExponentialDip" begin
    toa = default_toa()
    ctoa = TOACorrection()

    toa2 = TOA(
        time(Double64((54954.0 - epoch_mjd) * day_to_s)),
        time(1e-6),
        frequency(2.5e9),
        dimensionless(Double64(0.0)),
        default_ephem(),
        1,
    )
    ctoa2 = TOACorrection()

    params = (
        EXPDIPFREF = frequency(1.0e9),
        EXPDIPEPS = time(1e-4),
        EXPDIPEP = (time(day_to_s * (54952.92239 - epoch_mjd)),),
        EXPDIPAMP_ = (time(1e-6),),
        EXPDIPIDX_ = (dimensionless(4.0),),
        EXPDIPTAU_ = (time(10.0 * day_to_s),),
    )

    expdip = ChromaticExponentialDip()

    @test delay(expdip, toa, ctoa, params) == dimensionless(0.0)
    @test delay(expdip, toa2, ctoa2, params) != dimensionless(0.0)

    @test @ballocated(delay($expdip, $toa2, $ctoa2, $params)) == 0
end
