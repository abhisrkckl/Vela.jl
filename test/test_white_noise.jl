@testset "MeasurementNoise" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (EFAC = (dimensionless(1.1),),)

    wn = MeasurementNoise(UInt[1], UInt[0])
    @test efac(wn, toa, params) == params.EFAC[1]
    @test equad2(wn, toa, params) == time(0.0)^Val(2)
    ctoa1 = correct_toa(wn, toa, ctoa, params)

    @ballocated(correct_toa($wn, $toa, $ctoa, $params)) == 0

    @test ctoa1.delay == ctoa.delay
    @test ctoa1.phase == ctoa.phase
    @test ctoa1.doppler == ctoa.doppler
    @test ctoa1.spin_frequency == ctoa.spin_frequency
    @test scaled_toa_error_sqr(toa, ctoa1) > toa.error^Val(2)
    display(wn)
end

@testset "DispersionMeasurementNoise" begin
    wtoa = default_wbtoa()
    cwtoa = WidebandTOACorrection()

    params = (DMEFAC = (dimensionless(1.1),), DMEQUAD = (GQ{-1}(1e9),))

    dmwn = DispersionMeasurementNoise(UInt[1], UInt[1])

    @test dmefac(dmwn, wtoa, params) == params.DMEFAC[1]
    @test dmequad2(dmwn, wtoa, params) == params.DMEQUAD[1]^Val(2)
    cwtoa1 = correct_toa(dmwn, wtoa, cwtoa, params)

    @ballocated(correct_toa($dmwn, $wtoa, $cwtoa, $params)) == 0

    @test correct_toa(dmwn, wtoa.toa, cwtoa.toa_correction, params) ==
          correct_toa(cwtoa.toa_correction)

    @test cwtoa1.toa_correction.delay == cwtoa.toa_correction.delay
    @test cwtoa1.toa_correction.phase == cwtoa.toa_correction.phase
    @test cwtoa1.toa_correction.doppler == cwtoa.toa_correction.doppler
    @test cwtoa1.toa_correction.spin_frequency == cwtoa.toa_correction.spin_frequency
    @test scaled_dm_error_sqr(wtoa.dminfo, cwtoa1.dm_correction) > wtoa.dminfo.error^Val(2)

    display(dmwn)
end
