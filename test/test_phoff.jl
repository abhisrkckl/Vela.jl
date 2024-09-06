@testset "PhaseOffset" begin
    toa = default_toa()
    ctoa = CorrectedTOA(toa)

    tzrtoa = default_tzrtoa()
    ctzrtoa = CorrectedTOA(tzrtoa)

    params = (PHOFF = dimensionless(1e-6),)

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
