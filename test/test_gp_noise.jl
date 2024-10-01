@testset "PowerlawRedNoiseGP" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (
        TNREDAMP = dimensionless(-13.5),
        TNREDGAM = dimensionless(3.5),
        PLREDFREQ = frequency(1e-8),
        PLREDEPOCH = time(day_to_s * (54000.0 - epoch_mjd)),
        PLREDSIN_ = map(dimensionless, (1.1, 0.12, -0.23)),
        PLREDCOS_ = map(dimensionless, (1.5, -1.11, -0.8)),
    )

    rn = PowerlawRedNoiseGP(3)
    @test length(rn.ln_js) == 3

    @test isfinite(delay(rn, toa, ctoa, params))
end
