@testset "BinaryELL1" begin
    toa = default_toa()
    ctoa = TOACorrection()

    tzrtoa = default_tzrtoa()
    ctzrtoa = TOACorrection()

    params = (
        TASC = time((53470.0 - epoch_mjd) * day_to_s),
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
    )

    for use_fbx in [true, false]
        ell1 = BinaryELL1(use_fbx)
        display(ell1)
        ctoa_1 = correct_toa(ell1, toa, ctoa, params)
        @test isfinite(ctoa_1.delay) && isfinite(ctoa_1.doppler)
        @test @ballocated(correct_toa($ell1, $toa, $ctoa, $params)) == 0
    end
end

@testset "BinaryELL1H" begin
    toa = default_toa()
    ctoa = TOACorrection()

    tzrtoa = default_tzrtoa()
    ctzrtoa = TOACorrection()

    params = (
        TASC = time((53470.0 - epoch_mjd) * day_to_s),
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
        H3 = time(1e-9),
        STIGMA = dimensionless(0.02),
    )

    for use_fbx in [true, false]
        ell1h = BinaryELL1H(use_fbx)
        display(ell1h)
        ctoa_1 = correct_toa(ell1h, toa, ctoa, params)
        @test isfinite(ctoa_1.delay) && isfinite(ctoa_1.doppler)
        @test @ballocated(correct_toa($ell1h, $toa, $ctoa, $params)) == 0
    end
end

@testset "BinaryELL1k" begin
    toa = default_toa()
    ctoa = TOACorrection()

    tzrtoa = default_tzrtoa()
    ctzrtoa = TOACorrection()

    params = (
        TASC = time((53470.0 - epoch_mjd) * day_to_s),
        PB = time(8e4),
        PBDOT = dimensionless(1e-10),
        XPBDOT = dimensionless(0.0),
        FB = (frequency(1.25e-5), GQ{-2}(-1.5625e-20)),
        A1 = distance(5.0),
        A1DOT = dimensionless(0.0),
        EPS1 = dimensionless(1e-5),
        EPS2 = dimensionless(-2e-5),
        OMDOT = frequency(4e-9),
        LNEDOT = frequency(0.0),
        M2 = mass(5e-9),
        SINI = dimensionless(0.5),
    )

    for use_fbx in [true, false]
        ell1k = BinaryELL1k(use_fbx)
        display(ell1k)
        ctoa_1 = correct_toa(ell1k, toa, ctoa, params)
        @test isfinite(ctoa_1.delay) && isfinite(ctoa_1.doppler)
        @test @ballocated(correct_toa($ell1k, $toa, $ctoa, $params)) == 0
    end
end
