@testset "components" begin
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

    params = (
        PHOFF = dimensionless(1e-6),
        PEPOCH = time(53470.0 * day_to_s),
        F_ = frequency(100.0),
        F = (frequency(0.0), GQ(-1e-14, -2)),
        DMEPOCH = time(53470.0 * day_to_s),
        DM = (GQ(10.0, -1), GQ(1e-4, -2)),
        POSEPOCH = time(53470.0 * day_to_s),
        ELAT = dimensionless(1.2),
        ELONG = dimensionless(1.25),
        PX = GQ(3e-12, -1),
        PMELAT = GQ(-7e-16, -1),
        PMELONG = GQ(-5e-16, -1),
        EFAC = (dimensionless(1.1),),
        JUMP = (time(1e-6), time(1.2e-6)),
    )

    @testset "SolarSystem" begin
        ss = SolarSystem(true, true)
        @test ss.ecliptic_coordinates && ss.planet_shapiro

        ctoa1 = correct_toa(ss, ctoa, params)

        @test ctoa1.phase == ctoa.phase
        @test ctoa1.delay != ctoa.delay
        @test ctoa.doppler == 0 && ctoa1.doppler != 0
        @test !ctoa.barycentered && ctoa1.barycentered
        @test ctoa1.level == ctoa.level + 1

        ctoa2 = correct_toa(ss, ctoa1, params)
        @test (ctoa2.delay == ctoa1.delay) && (ctoa2.doppler == ctoa1.doppler)

        display(ss)
    end

    @testset "SolarWindDispersion" begin
        @test_throws AssertionError SolarWindDispersion(2)

        swd = SolarWindDispersion(0)
        @test dispersion_slope(swd, toa, params) == GQ(0.0, -1)
    end

    @testset "DispersionTaylor" begin
        dmt = DispersionTaylor()
        @test dispersion_slope(dmt, ctoa, params) == GQ(10.0, -1)
        @test delay(dmt, ctoa, params) ==
              dispersion_slope(dmt, ctoa, params) / ctoa.toa.observing_frequency^2
    end

    @testset "PhaseOffset" begin
        poff = PhaseOffset()
        @test phase(poff, ctoa, params) == dimensionless(-1e-6)
        @test phase(poff, ctzrtoa, params) == dimensionless(0.0)

        ctoa1 = correct_toa(poff, ctoa, params)
        @test ctoa1.delay == ctoa.delay
        @test ctoa1.phase ≈ ctoa.phase + phase(poff, ctoa, params)
        @test ctoa1.efac == ctoa.efac
        @test ctoa1.equad2 == ctoa.equad2
        @test ctoa1.spin_frequency == ctoa.spin_frequency == frequency(-1.0)

        ctzrtoa1 = correct_toa(poff, ctzrtoa, params)
        @test ctzrtoa1.delay == ctzrtoa.delay
        @test ctzrtoa1.phase == ctzrtoa.phase
        @test ctzrtoa1.efac == ctzrtoa.efac
        @test ctzrtoa1.equad2 == ctzrtoa.equad2
    end

    @testset "Spindown" begin
        spn = Spindown()
        @test phase(spn, ctoa, params) == dimensionless(0.0)
        @test spin_frequency(spn, ctoa, params) == frequency(100.0)

        ctoa1 = correct_toa(spn, ctoa, params)
        @test ctoa1.delay == ctoa.delay
        @test ctoa1.phase == ctoa.phase + phase(spn, ctoa, params)
        @test ctoa1.doppler == ctoa.doppler
        @test ctoa.spin_frequency == frequency(-1.0) &&
              ctoa1.spin_frequency > frequency(0.0)

    end

    @testset "PhaseJump" begin
        jump_mask = Matrix([1 0 0; 0 1 0])
        pjmp = PhaseJump(jump_mask)

        @test phase(pjmp, ctzrtoa, params) == 0
        @test phase(pjmp, ctoa, params) ≈ params.JUMP[1] * params.F_

        display(pjmp)
    end

    # @testset "Troposphere" begin
    #     tropo = Troposphere()
    #     @test delay(tropo, toa, params) == time(0.0)
    # end

    @testset "MeasurementNoise" begin
        wn = MeasurementNoise(UInt[1], UInt[0])
        @test efac(wn, ctoa, params) == params.EFAC[1]
        @test equad2(wn, ctoa, params) == time(0.0)^2
        ctoa1 = correct_toa(wn, ctoa, params)
        @test ctoa1.delay == ctoa.delay
        @test ctoa1.phase == ctoa.phase
        @test ctoa1.doppler == ctoa.doppler
        @test ctoa1.spin_frequency == ctoa.spin_frequency
        @test scaled_toa_error_sqr(ctoa1) > ctoa1.toa.error^2
        display(wn)
    end
end
