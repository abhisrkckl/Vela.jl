@testset "SolarWindDispersion" begin
    ephem = SolarSystemEphemeris(
        ssb_obs_pos,
        ssb_obs_vel,
        obs_sun_pos,
        obs_jupiter_pos,
        obs_saturn_pos,
        obs_venus_pos,
        obs_uranus_pos,
        obs_neptune_pos,
        obs_earth_pos,
    )

    toa = TOA(
        time(Double64(53470.0 * day_to_s)),
        time(1e-6),
        frequency(2.5e9),
        dimensionless(Double64(0.0)),
        false,
        ephem,
        1,
    )
    ctoa = CorrectedTOA(toa)

    tzrtoa =
        make_tzr_toa(time(Double64(53475.0 * day_to_s)), frequency(2.5e9), false, ephem)
    ctzrtoa = CorrectedTOA(tzrtoa)

    dminfo = DMInfo(GQ{-1}(1e16), GQ{-1}(1e11))
    wtoa = WidebandTOA(toa, dminfo)
    cwtoa = CorrectedWidebandTOA(wtoa)

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
    @test_throws AssertionError dispersion_slope(swd, ctoa, params)

    ss = SolarSystem(true, true)
    ctoa1 = correct_toa(ss, ctoa, params)
    @test dispersion_slope(swd, ctoa1, params) != GQ{-1}(0.0)

    @test @ballocated(delay($swd, $ctoa1, $params)) == 0

    cwtoa1 = correct_toa(ss, cwtoa, params)
    cwtoa2 = correct_toa(swd, cwtoa1, params)
    @test cwtoa2.corrected_dminfo != cwtoa.corrected_dminfo

    @test @ballocated(correct_toa($swd, $cwtoa1, $params)) == 0
end