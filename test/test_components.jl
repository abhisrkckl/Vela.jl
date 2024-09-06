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

    dminfo = DMInfo(GQ{-1}(1e16), GQ{-1}(1e11))
    wtoa = WidebandTOA(toa, dminfo)
    cwtoa = CorrectedWidebandTOA(wtoa)

    params = (
        PHOFF = dimensionless(1e-6),
        PEPOCH = time(53470.0 * day_to_s),
        F_ = frequency(100.0),
        F = (frequency(0.0), GQ{-2}(-1e-14)),
        DMEPOCH = time(53470.0 * day_to_s),
        DM = (GQ{-1}(4e16), GQ{-2}(1e11)),
        FDJUMPDM = (GQ{-1}(1e11), GQ{-1}(-1e11)),
        POSEPOCH = time(53470.0 * day_to_s),
        ELAT = dimensionless(1.2),
        ELONG = dimensionless(1.25),
        PX = GQ{-1}(3e-12),
        PMELAT = GQ{-1}(-7e-16),
        PMELONG = GQ{-1}(-5e-16),
        EFAC = (dimensionless(1.1),),
        JUMP = (time(1e-6), time(1.2e-6)),
        NE_SW = GQ{-2}(1.6e8),
        CMEPOCH = time(53470.0 * day_to_s),
        CM = (GQ{1}(4e4), GQ{0}(1e-1)),
        TASC = time(53470.0 * day_to_s),
        PB = time(8e4),
        PBDOT = dimensionless(1e-10),
        XPBDOT = dimensionless(0.0),
        FB = (frequency(1.25e-5), GQ{-2}(-1.5625e-20)),
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
        H3 = time(1e-9),
        STIGMA = dimensionless(0.02),
        SHAPMAX = dimensionless(5.0),
        FD = (time(0.1), time(0.2)),
    )

    include("test_solarsystem.jl")

    include("test_solarwind.jl")

    include("test_dispersion.jl")

    include("test_chromatic.jl")

    include("test_fd.jl")

    include("test_wavex.jl")

    @testset "BinaryELL1" begin
        for use_fbx in [true, false]
            ell1 = BinaryELL1(use_fbx)
            display(ell1)
            ctoa_1 = correct_toa(ell1, ctoa, params)
            @test isfinite(ctoa_1.delay) && isfinite(ctoa_1.doppler)
            @test @ballocated(correct_toa($ell1, $ctoa, $params)) == 0
        end
    end

    @testset "BinaryELL1H" begin
        for use_fbx in [true, false]
            ell1h = BinaryELL1H(use_fbx)
            display(ell1h)
            ctoa_1 = correct_toa(ell1h, ctoa, params)
            @test isfinite(ctoa_1.delay) && isfinite(ctoa_1.doppler)
            @test @ballocated(correct_toa($ell1h, $ctoa, $params)) == 0
        end
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

        for use_fbx in [true, false]
            dd = BinaryDD(use_fbx)
            ctoa_1 = correct_toa(dd, ctoa1, params)
            @test isfinite(ctoa_1.delay) && isfinite(ctoa_1.doppler)
            @test @ballocated(correct_toa($dd, $ctoa1, $params)) == 0
            display(dd)
        end
    end

    @testset "BinaryDDH" begin
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

        for use_fbx in [true, false]
            ddh = BinaryDDH(use_fbx)
            ctoa_1 = correct_toa(ddh, ctoa1, params)
            @test isfinite(ctoa_1.delay) && isfinite(ctoa_1.doppler)
            @test @ballocated(correct_toa($ddh, $ctoa1, $params)) == 0
            display(ddh)
        end
    end

    @testset "BinaryDDS" begin
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

        for use_fbx in [true, false]
            dds = BinaryDDS(use_fbx)
            ctoa_1 = correct_toa(dds, ctoa1, params)
            @test isfinite(ctoa_1.delay) && isfinite(ctoa_1.doppler)
            @test @ballocated(correct_toa($dds, $ctoa1, $params)) == 0
            display(dds)
        end
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
        @test ctoa1.spin_frequency == ctoa.spin_frequency == frequency(0.0)

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
        @test ctoa.spin_frequency == frequency(0.0) && ctoa1.spin_frequency > frequency(0.0)

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
end
