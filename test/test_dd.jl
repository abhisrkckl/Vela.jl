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
        time(Double64((53471.0 - epoch_mjd) * day_to_s)),
        time(1e-6),
        frequency(2.5e9),
        dimensionless(Double64(0.0)),
        default_ephem(),
        1,
    )
    ctoa1 = TOACorrection()

    params = (
        T0 = time((53470.0 - epoch_mjd) * day_to_s),
        PB = time(8e4),
        PBDOT = dimensionless(1e-10),
        XPBDOT = dimensionless(0.0),
        FB = (frequency(1.25e-5), GQ{-2}(-1.5625e-20)),
        A1 = distance(5.0),
        A1DOT = dimensionless(0.0),
        ECC = dimensionless(0.5),
        EDOT = frequency(0.0),
        OM = dimensionless(0.1),
        OMDOT = frequency(0.0),
        DR = dimensionless(0.0),
        DTH = dimensionless(0.0),
        GAMMA = time(0.0),
        M2 = mass(5e-9),
        SINI = dimensionless(0.5),
    )

    for use_fbx in [true, false]
        dd = BinaryDD(use_fbx)
        ctoa_1 = correct_toa(dd, toa1, ctoa1, params)
        @test isfinite(ctoa_1.delay) && isfinite(ctoa_1.doppler)
        @test @ballocated(correct_toa($dd, $toa1, $ctoa1, $params)) == 0
        display(dd)
    end
end

@testset "BinaryDDH" begin
    toa1 = TOA(
        time(Double64((53471.0 - epoch_mjd) * day_to_s)),
        time(1e-6),
        frequency(2.5e9),
        dimensionless(Double64(0.0)),
        default_ephem(),
        1,
    )
    ctoa1 = TOACorrection()

    params = (
        T0 = time((53470.0 - epoch_mjd) * day_to_s),
        PB = time(8e4),
        PBDOT = dimensionless(1e-10),
        XPBDOT = dimensionless(0.0),
        FB = (frequency(1.25e-5), GQ{-2}(-1.5625e-20)),
        A1 = distance(5.0),
        A1DOT = dimensionless(0.0),
        ECC = dimensionless(0.5),
        EDOT = frequency(0.0),
        OM = dimensionless(0.1),
        OMDOT = frequency(0.0),
        DR = dimensionless(0.0),
        DTH = dimensionless(0.0),
        GAMMA = time(0.0),
        H3 = time(1e-9),
        STIGMA = dimensionless(0.02),
    )

    for use_fbx in [true, false]
        ddh = BinaryDDH(use_fbx)
        ctoa_1 = correct_toa(ddh, toa1, ctoa1, params)
        @test isfinite(ctoa_1.delay) && isfinite(ctoa_1.doppler)
        @test @ballocated(correct_toa($ddh, $toa1, $ctoa1, $params)) == 0
        display(ddh)
    end
end

@testset "BinaryDDS" begin
    toa1 = TOA(
        time(Double64((53471.0 - epoch_mjd) * day_to_s)),
        time(1e-6),
        frequency(2.5e9),
        dimensionless(Double64(0.0)),
        default_ephem(),
        1,
    )
    ctoa1 = TOACorrection()

    params = (
        T0 = time((53470.0 - epoch_mjd) * day_to_s),
        PB = time(8e4),
        PBDOT = dimensionless(1e-10),
        XPBDOT = dimensionless(0.0),
        FB = (frequency(1.25e-5), GQ{-2}(-1.5625e-20)),
        A1 = distance(5.0),
        A1DOT = dimensionless(0.0),
        ECC = dimensionless(0.5),
        EDOT = frequency(0.0),
        OM = dimensionless(0.1),
        OMDOT = frequency(0.0),
        DR = dimensionless(0.0),
        DTH = dimensionless(0.0),
        GAMMA = time(0.0),
        M2 = mass(5e-9),
        SHAPMAX = dimensionless(5.0),
    )

    for use_fbx in [true, false]
        dds = BinaryDDS(use_fbx)
        ctoa_1 = correct_toa(dds, toa1, ctoa1, params)
        @test isfinite(ctoa_1.delay) && isfinite(ctoa_1.doppler)
        @test @ballocated(correct_toa($dds, $toa1, $ctoa1, $params)) == 0
        display(dds)
    end
end

@testset "BinaryDDK" begin
    toa1 = TOA(
        time(Double64((53471.0 - epoch_mjd) * day_to_s)),
        time(1e-6),
        frequency(2.5e9),
        dimensionless(Double64(0.0)),
        default_ephem(),
        1,
    )
    ctoa1 = TOACorrection()

    params_binary = (
        T0 = time((53470.0 - epoch_mjd) * day_to_s),
        PB = time(8e4),
        PBDOT = dimensionless(1e-10),
        XPBDOT = dimensionless(0.0),
        FB = (frequency(1.25e-5), GQ{-2}(-1.5625e-20)),
        A1 = distance(5.0),
        A1DOT = dimensionless(0.0),
        ECC = dimensionless(0.5),
        EDOT = frequency(0.0),
        OM = dimensionless(0.1),
        OMDOT = frequency(0.0),
        DR = dimensionless(0.0),
        DTH = dimensionless(0.0),
        GAMMA = time(0.0),
        M2 = mass(5e-9),
        KIN = dimensionless(0.5),
        KOM = dimensionless(0.15),
    )

    params_ecl = (
        POSEPOCH = time((53470.0 - epoch_mjd) * day_to_s),
        ELAT = dimensionless(1.2),
        ELONG = dimensionless(1.25),
        PX = GQ{-1}(3e-12),
        PMELAT = GQ{-1}(-7e-16),
        PMELONG = GQ{-1}(-5e-16),
    )

    params_eql = (
        POSEPOCH = time((53470.0 - epoch_mjd) * day_to_s),
        RAJ = dimensionless(1.2),
        DECJ = dimensionless(1.25),
        PX = GQ{-1}(3e-12),
        PMRA = GQ{-1}(-7e-16),
        PMDEC = GQ{-1}(-5e-16),
    )

    for ecl in [true, false]
        params = merge(params_binary, (ecl ? params_ecl : params_eql))
        for use_fbx in [true, false]
            ss = SolarSystem(ecl, false)
            ctoa2 = correct_toa(ss, toa1, ctoa1, params)

            ddk = BinaryDDK(use_fbx, ecl)
            ctoa3 = correct_toa(ddk, toa1, ctoa2, params)
            @test isfinite(ctoa3.delay) && isfinite(ctoa3.doppler)
            @test @ballocated(correct_toa($ddk, $toa1, $ctoa2, $params)) == 0
            display(ddk)
        end
    end
end
