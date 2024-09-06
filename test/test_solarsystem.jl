@testset "SolarSystem" begin
    toa = default_toa()
    ctoa = CorrectedTOA(toa)

    tzrtoa = default_tzrtoa()
    ctzrtoa = CorrectedTOA(tzrtoa)

    wtoa = default_wbtoa()
    cwtoa = CorrectedWidebandTOA(wtoa)

    params_ecl = (
        POSEPOCH = time(53470.0 * day_to_s),
        ELAT = dimensionless(1.2),
        ELONG = dimensionless(1.25),
        PX = GQ{-1}(3e-12),
        PMELAT = GQ{-1}(-7e-16),
        PMELONG = GQ{-1}(-5e-16),
    )

    params_eql = (
        POSEPOCH = time(53470.0 * day_to_s),
        RAJ = dimensionless(1.2),
        DECJ = dimensionless(1.25),
        PX = GQ{-1}(3e-12),
        PMRA = GQ{-1}(-7e-16),
        PMDEC = GQ{-1}(-5e-16),
    )

    for (ecl, params) in zip((true, false), (params_ecl, params_eql))
        for planet_shapiro in (true, false)
            ss = SolarSystem(ecl, planet_shapiro)
            @test ss.ecliptic_coordinates == ecl && ss.planet_shapiro == planet_shapiro

            ctoa1 = correct_toa(ss, ctoa, params)

            @test @ballocated(correct_toa($ss, $ctoa, $params)) == 0

            @test ctoa1.phase == ctoa.phase
            @test ctoa1.delay != ctoa.delay
            @test ctoa.doppler == 0 && ctoa1.doppler != 0
            @test !ctoa.barycentered && ctoa1.barycentered
            @test ctoa1.level == ctoa.level + 1

            ctoa2 = correct_toa(ss, ctoa1, params)
            @test (ctoa2.delay == ctoa1.delay) && (ctoa2.doppler == ctoa1.doppler)

            cwtoa1 = correct_toa(ss, cwtoa, params)
            @test cwtoa1.corrected_dminfo == cwtoa.corrected_dminfo

            @test @ballocated(correct_toa($ss, $cwtoa, $params)) == 0

            display(ss)
        end
    end
end