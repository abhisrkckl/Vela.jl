@testset "MeasurementNoise" begin
    toa = default_toa()
    ctoa = CorrectedTOA(toa)

    params = (EFAC = (dimensionless(1.1),),)

    wn = MeasurementNoise(UInt[1], UInt[0])
    @test efac(wn, ctoa, params) == params.EFAC[1]
    @test equad2(wn, ctoa, params) == time(0.0)^Val(2)
    ctoa1 = correct_toa(wn, ctoa, params)

    @ballocated(correct_toa($wn, $ctoa, $params)) == 0

    @test ctoa1.delay == ctoa.delay
    @test ctoa1.phase == ctoa.phase
    @test ctoa1.doppler == ctoa.doppler
    @test ctoa1.spin_frequency == ctoa.spin_frequency
    @test scaled_toa_error_sqr(ctoa1) > ctoa1.toa.error^Val(2)
    display(wn)
end
