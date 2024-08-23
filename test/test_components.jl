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
        DM = (GQ(4e16, -1), GQ(1e11, -2)),
        FDJUMPDM = (GQ(1e11, -1), GQ(-1e11, -1)),
        POSEPOCH = time(53470.0 * day_to_s),
        ELAT = dimensionless(1.2),
        ELONG = dimensionless(1.25),
        PX = GQ(3e-12, -1),
        PMELAT = GQ(-7e-16, -1),
        PMELONG = GQ(-5e-16, -1),
        EFAC = (dimensionless(1.1),),
        JUMP = (time(1e-6), time(1.2e-6)),
        WXEPOCH = time(53470.0 * day_to_s),
        WXFREQ_ = (frequency(1e-9), frequency(2e-9), frequency(3e-9)),
        WXSIN_ = (time(1.2e-6), time(5.1e-7), time(2.5e-7)),
        WXCOS_ = (time(-1.3e-6), time(5.2e-7), time(2.6e-7)),
        DMWXEPOCH = time(53470.0 * day_to_s),
        DMWXFREQ_ = (frequency(1e-9), frequency(2e-9), frequency(3e-9)),
        DMWXSIN_ = (GQ(1.2e+4, -1), GQ(5.1e+3, -1), GQ(2.5e+2, -1)),
        DMWXCOS_ = (GQ(-1.3e+4, -1), GQ(5.2e+3, -1), GQ(2.6e+3, -1)),
        CMWXEPOCH = time(53470.0 * day_to_s),
        CMWXFREQ_ = (frequency(1e-9), frequency(2e-9), frequency(3e-9)),
        CMWXSIN_ = (GQ(1.2e+4, 1), GQ(5.1e+3, 1), GQ(2.5e+2, 1)),
        CMWXCOS_ = (GQ(-1.3e+4, 1), GQ(5.2e+3, 1), GQ(2.6e+3, 1)),
        NE_SW = GQ(1.6e8, -2),
        TNCHROMIDX = dimensionless(2.0),
        CMEPOCH = time(53470.0 * day_to_s),
        CM = (GQ(4e4, 1), GQ(1e-1, 0)),
        TASC = time(53470.0 * day_to_s),
        PB = time(8e4),
        PBDOT = dimensionless(1e-10),
        XPBDOT = dimensionless(0.0),
        FB = (frequency(1.25e-5), GQ(-1.5625e-20, -2)),
        A1 = distance(5.0),
        A1DOT = dimensionless(0.0),
        EPS1 = dimensionless(1e-5),
        EPS2 = dimensionless(-2e-5),
        EPS1DOT = frequency(0.0),
        EPS2DOT = frequency(0.0),
        M2 = mass(5e-9),
        SINI = dimensionless(0.5),
        ECC = dimensionless(0.5),
        EDOT = frequency(0.0),
        OM = dimensionless(0.1),
        OMDOT = frequency(0.0),
        T0 = time(53470.0 * day_to_s),
        DR = dimensionless(0.0),
        DTH = dimensionless(0.0),
        GAMMA = time(0.0),
        FD = (time(0.1), time(0.2)),
    )

    @testset "SolarSystem" begin
        ss = SolarSystem(true, true)
        @test ss.ecliptic_coordinates && ss.planet_shapiro

        ctoa1 = correct_toa(ss, ctoa, params)

        @test @ballocated(correct_toa($ss, $ctoa, $params)) == 0

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
        swd = SolarWindDispersion()
        @test_throws AssertionError dispersion_slope(swd, ctoa, params)
        ss = SolarSystem(true, true)
        ctoa1 = correct_toa(ss, ctoa, params)
        @test dispersion_slope(swd, ctoa1, params) != GQ(0.0, -1)

        @test @ballocated(delay($swd, $ctoa1, $params)) == 0
    end

    @testset "DispersionTaylor" begin
        dmt = DispersionTaylor()
        @test dispersion_slope(dmt, ctoa, params) == params.DM[1]
        @test delay(dmt, ctoa, params) ==
              dispersion_slope(dmt, ctoa, params) / ctoa.toa.observing_frequency^2

        @test @ballocated(delay($dmt, $ctoa, $params)) == 0
    end

    @testset "DispersionOffset" begin
        jump_mask = BitMatrix([1 0 0; 0 1 0])
        dmoff = DispersionOffset(jump_mask)

        @test dispersion_slope(dmoff, ctzrtoa, params) == GQ(0.0, -1)
        @test dispersion_slope(dmoff, ctoa, params) == -params.FDJUMPDM[1]

        @test @ballocated(delay($dmoff, $ctoa, $params)) == 0

        display(dmoff)

        jump_mask_ex = [1, 2, 0]
        dmoff_ex = ExclusiveDispersionOffset(jump_mask_ex)

        @test dispersion_slope(dmoff_ex, ctzrtoa, params) == GQ(0.0, -1)
        @test dispersion_slope(dmoff_ex, ctoa, params) == -params.FDJUMPDM[1]

        toa1 = TOA(
            time(Double64(53470.0 * day_to_s)),
            time(1e-6),
            frequency(2.5e9),
            dimensionless(Double64(0.0)),
            false,
            ephem,
            3,
        )
        ctoa1 = CorrectedTOA(toa1)
        @test dispersion_slope(dmoff_ex, ctoa1, params) ≈ GQ(0.0, -1)

        display(dmoff_ex)
    end

    @testset "ChromaticTaylor" begin
        cmt = ChromaticTaylor()
        dmt = DispersionTaylor()
        @test chromatic_slope(cmt, ctoa, params) == params.CM[1]

        # When TNCHROMIDX == 2, DM and CM that have equal numerical values
        # in the par file must give equal delays. In our units, this corresponds
        # to `value(DM) / value(CM) == 1e12`. See the `params` tuple above.
        @test delay(dmt, ctoa, params) == delay(cmt, ctoa, params)

        @test @ballocated(delay($cmt, $ctoa, $params)) == 0
    end

    @testset "FrequencyDependent" begin
        fd = FrequencyDependent()
        @test isfinite(delay(fd, ctoa, params))
        @test @ballocated(delay($fd, $ctoa, $params)) == 0
    end

    @testset "WaveX" begin
        wx = WaveX()
        @test isfinite(delay(wx, ctoa, params))
        @test @ballocated(delay($wx, $ctoa, $params)) == 0
    end

    @testset "DMWaveX" begin
        dmwx = DMWaveX()
        @test isfinite(delay(dmwx, ctoa, params))
        @test @ballocated(delay($dmwx, $ctoa, $params)) == 0
    end

    @testset "CMWaveX" begin
        cmwx = CMWaveX()
        @test isfinite(delay(cmwx, ctoa, params))
        @test @ballocated(delay($cmwx, $ctoa, $params)) == 0
    end

    @testset "BinaryELL1" begin
        ell1 = BinaryELL1(true)
        @test isfinite(delay(ell1, ctoa, params))
        display(ell1)

        ell1 = BinaryELL1(false)
        @test isfinite(delay(ell1, ctoa, params))
        display(ell1)

        @test @ballocated(delay($ell1, $ctoa, $params)) == 0

        params1 = (
            TASC = time(53470.0 * day_to_s),
            PB = time(8e4),
            PBDOT = dimensionless(1e-10),
            A1 = distance(5.0),
            A1DOT = dimensionless(0.0),
            EPS1 = dimensionless(1e-5),
            EPS2 = dimensionless(-2e-5),
            EPS1DOT = frequency(0.0),
            EPS2DOT = frequency(0.0),
        )
        @test isfinite(delay(ell1, ctoa, params1))
    end

    @testset "BinaryDD" begin
        @testset "mikkola" begin
            kepler = (u, e) -> u - e * sin(u)
            us = [-π / 4, 0.0, π / 4, 3 * π / 4, 4 * π / 3, 7 * π / 3]
            es = [0.0, 0.5, 0.3]
            for e in es
                for u in us
                    @test u ≈ Vela.mikkola(kepler(u, e), e)
                end
            end
        end

        toa1 = TOA(
            time(Double64(53471.0 * day_to_s)),
            time(1e-6),
            frequency(2.5e9),
            dimensionless(Double64(0.0)),
            false,
            ephem,
            1,
        )
        ctoa1 = CorrectedTOA(toa1)

        dd = BinaryDD(true)
        @test isfinite(delay(dd, ctoa1, params))
        display(dd)

        dd = BinaryDD(false)
        @test isfinite(delay(dd, ctoa1, params))
        display(dd)

        @test @ballocated(delay($dd, $ctoa1, $params)) == 0

        # params1 = (
        #     TASC = time(53470.0 * day_to_s),
        #     PB = time(8e4),
        #     PBDOT = dimensionless(1e-10),
        #     A1 = distance(5.0),
        #     A1DOT = dimensionless(0.0),
        #     EPS1 = dimensionless(1e-5),
        #     EPS2 = dimensionless(-2e-5),
        #     EPS1DOT = frequency(0.0),
        #     EPS2DOT = frequency(0.0),
        # )
        # @test isfinite(delay(ell1, ctoa, params1))
    end

    @testset "PhaseOffset" begin
        poff = PhaseOffset()
        @test phase(poff, ctoa, params) == dimensionless(-1e-6)
        @test phase(poff, ctzrtoa, params) == dimensionless(0.0)

        @test @ballocated(phase($poff, $ctoa, $params)) == 0

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

        @ballocated(correct_toa($spn, $ctoa, $params)) == 0
    end

    @testset "PhaseJump" begin
        jump_mask = BitMatrix([1 0 0; 0 1 0])
        pjmp = PhaseJump(jump_mask)

        @test phase(pjmp, ctzrtoa, params) == 0
        @test phase(pjmp, ctoa, params) ≈ params.JUMP[1] * params.F_

        @ballocated(phase($pjmp, $ctzrtoa, $params)) == 0

        display(pjmp)

        jump_mask_ex = [1, 2, 0]
        pjmp_ex = ExclusivePhaseJump(jump_mask_ex)

        @test phase(pjmp_ex, ctzrtoa, params) == 0
        @test phase(pjmp_ex, ctoa, params) ≈ params.JUMP[1] * params.F_

        toa1 = TOA(
            time(Double64(53470.0 * day_to_s)),
            time(1e-6),
            frequency(2.5e9),
            dimensionless(Double64(0.0)),
            false,
            ephem,
            3,
        )
        ctoa1 = CorrectedTOA(toa1)
        @test phase(pjmp_ex, ctoa1, params) ≈ dimensionless(0.0)

        display(pjmp_ex)
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

        @ballocated(correct_toa($wn, $ctoa, $params)) == 0

        @test ctoa1.delay == ctoa.delay
        @test ctoa1.phase == ctoa.phase
        @test ctoa1.doppler == ctoa.doppler
        @test ctoa1.spin_frequency == ctoa.spin_frequency
        @test scaled_toa_error_sqr(ctoa1) > ctoa1.toa.error^2
        display(wn)
    end
end
