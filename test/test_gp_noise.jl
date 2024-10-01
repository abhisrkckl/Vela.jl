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

@testset "PowerlawDispersionNoiseGP" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (
        TNDMAMP = dimensionless(-13.5),
        TNDMGAM = dimensionless(3.5),
        PLDMFREQ = frequency(1e-8),
        PLDMEPOCH = time(day_to_s * (54000.0 - epoch_mjd)),
        PLDMSIN_ = map(dimensionless, (1.1, 0.12, -0.23)),
        PLDMCOS_ = map(dimensionless, (1.5, -1.11, -0.8)),
    )

    dmn = PowerlawDispersionNoiseGP(3)
    @test length(dmn.ln_js) == 3

    @test isfinite(delay(dmn, toa, ctoa, params))
end

@testset "PowerlawChromaticNoiseGP" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (
        TNCHROMIDX = dimensionless(4.0),
        TNCHROMAMP = dimensionless(-13.5),
        TNCHROMGAM = dimensionless(3.5),
        PLCHROMFREQ = frequency(1e-8),
        PLCHROMEPOCH = time(day_to_s * (54000.0 - epoch_mjd)),
        PLCHROMSIN_ = map(dimensionless, (1.1, 0.12, -0.23)),
        PLCHROMCOS_ = map(dimensionless, (1.5, -1.11, -0.8)),
    )

    cmn = PowerlawChromaticNoiseGP(3)
    @test length(cmn.ln_js) == 3

    @test isfinite(delay(cmn, toa, ctoa, params))
end
