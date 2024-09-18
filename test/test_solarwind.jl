@testset "SolarWindDispersion" begin
    toa = default_toa()
    ctoa = TOACorrection()

    tzrtoa = default_tzrtoa()
    ctzrtoa = TOACorrection()

    wtoa = default_wbtoa()
    cwtoa = WidebandTOACorrection()

    params = (
        POSEPOCH = time(53470.0 * day_to_s),
        ELAT = dimensionless(1.2),
        ELONG = dimensionless(1.25),
        PX = GQ{-1}(3e-12),
        PMELAT = GQ{-1}(-7e-16),
        PMELONG = GQ{-1}(-5e-16),
        NE_SW = GQ{-2}(1.6e8),
    )

    swd = SolarWindDispersion()
    @test_throws AssertionError dispersion_slope(swd, toa, ctoa, params)

    ss = SolarSystem(true, true)
    ctoa1 = correct_toa(ss, toa, ctoa, params)
    @test dispersion_slope(swd, toa, ctoa1, params) != GQ{-1}(0.0)

    @test @ballocated(delay($swd, $toa, $ctoa1, $params)) == 0

    cwtoa1 = correct_toa(ss, wtoa, cwtoa, params)
    cwtoa2 = correct_toa(swd, wtoa, cwtoa1, params)
    @test cwtoa2.dm_correction != cwtoa.dm_correction

    @test @ballocated(correct_toa($swd, $wtoa, $cwtoa1, $params)) == 0

    display(swd)
end
