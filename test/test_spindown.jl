@testset "Spindown" begin
    toa = default_toa()
    ctoa = CorrectedTOA(toa)

    tzrtoa = default_tzrtoa()
    ctzrtoa = CorrectedTOA(tzrtoa)

    params = (
        PEPOCH = time(53470.0 * day_to_s),
        F_ = frequency(100.0),
        F = (frequency(0.0), GQ{-2}(-1e-14)),
    )

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
