@testset "PhaseJump" begin
    toa = default_toa()
    ctoa = CorrectedTOA(toa)

    tzrtoa = default_tzrtoa()
    ctzrtoa = CorrectedTOA(tzrtoa)

    params = (
        PEPOCH = time(53470.0 * day_to_s),
        F_ = frequency(100.0),
        F = (frequency(0.0), GQ{-2}(-1e-14)),
        JUMP = (time(1e-6), time(1.2e-6)),
    )

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
        default_ephem(),
        3,
    )
    ctoa1 = CorrectedTOA(toa1)
    @test phase(pjmp_ex, ctoa1, params) ≈ dimensionless(0.0)

    display(pjmp_ex)
end
