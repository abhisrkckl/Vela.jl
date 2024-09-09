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

@testset "DispersionMeasurementNoise" begin
    wtoa = default_wbtoa()
    cwtoa = CorrectedWidebandTOA(wtoa)

    params = (DMEFAC = (dimensionless(1.1),), DMEQUAD = (GQ{-1}(1e9),))

    dmwn = DispersionMeasurementNoise(UInt[1], UInt[1])

    @test dmefac(dmwn, cwtoa, params) == params.DMEFAC[1]
    @test dmequad2(dmwn, cwtoa, params) == params.DMEQUAD[1]^Val(2)
    cwtoa1 = correct_toa(dmwn, cwtoa, params)

    @ballocated(correct_toa($dmwn, $cwtoa, $params)) == 0

    @test cwtoa1.corrected_toa.delay == cwtoa.corrected_toa.delay
    @test cwtoa1.corrected_toa.phase == cwtoa.corrected_toa.phase
    @test cwtoa1.corrected_toa.doppler == cwtoa.corrected_toa.doppler
    @test cwtoa1.corrected_toa.spin_frequency == cwtoa.corrected_toa.spin_frequency
    @test scaled_dm_error_sqr(cwtoa1.corrected_dminfo) >
          cwtoa1.corrected_dminfo.dminfo.error^Val(2)

    display(dmwn)
end
