@testset "PowerlawRedNoiseGP" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (
        TNREDAMP = dimensionless(-13.5),
        TNREDGAM = dimensionless(3.5),
        PLREDFREQ = frequency(1e-8),
        PLREDEPOCH = time(day_to_s * (54000.0 - epoch_mjd)),
        PLREDSIN_ = map(dimensionless, (1.3, 2.0, 1.1, 0.12, -0.23)),
        PLREDCOS_ = map(dimensionless, (-1.6, 2.3, 1.5, -1.11, -0.8)),
    )

    rn = PowerlawRedNoiseGP(3, 2, 2.0)
    @test length(rn.ln_js) == 5

    @test isfinite(delay(rn, toa, ctoa, params))

    @test all(calc_noise_weights_inv(rn, params) .> 0)
    @test calc_noise_weights_inv(rn, params)[1:5] ==
          calc_noise_weights_inv(rn, params)[6:10]

    w1 = calc_noise_weights_inv(rn, params)[1:5]
    w2 = [
        value(
            1 / Vela.powerlaw(
                10.0^params.TNREDAMP,
                params.TNREDGAM,
                exp(rn.ln_js[j]) * params.PLREDFREQ,
                params.PLREDFREQ,
            ),
        ) for j = 1:5
    ]
    @test all(isapprox.(w1 ./ w2, 1.0, atol = 1e-6))
end

@testset "PowerlawDispersionNoiseGP" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (
        TNDMAMP = dimensionless(-13.5),
        TNDMGAM = dimensionless(3.5),
        PLDMFREQ = frequency(1e-8),
        PLDMEPOCH = time(day_to_s * (54000.0 - epoch_mjd)),
        PLDMSIN_ = map(dimensionless, (3.0, 2.0, 1.1, 0.12, -0.23)),
        PLDMCOS_ = map(dimensionless, (3.0, 2.0, 1.5, -1.11, -0.8)),
    )

    dmn = PowerlawDispersionNoiseGP(3, 2, 2.0)
    @test length(dmn.ln_js) == 5

    @test isfinite(delay(dmn, toa, ctoa, params))

    @test all(calc_noise_weights_inv(dmn, params) .> 0)
    @test calc_noise_weights_inv(dmn, params)[1:5] ==
          calc_noise_weights_inv(dmn, params)[6:10]

    @test calc_noise_weights_inv(dmn, params)[1:5] ≈ [
        value(
            1 / Vela.powerlaw(
                10.0^params.TNDMAMP,
                params.TNDMGAM,
                exp(dmn.ln_js[j]) * params.PLDMFREQ,
                params.PLDMFREQ,
            ),
        ) for j = 1:5
    ]
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
        PLCHROMSIN_ = map(dimensionless, (1.3, 2.0, 1.1, 0.12, -0.23)),
        PLCHROMCOS_ = map(dimensionless, (1.3, -2.0, 1.5, -1.11, -0.8)),
    )

    cmn = PowerlawChromaticNoiseGP(3, 2, 2.0)
    @test length(cmn.ln_js) == 5

    @test isfinite(delay(cmn, toa, ctoa, params))

    @test all(calc_noise_weights_inv(cmn, params) .> 0)
    @test calc_noise_weights_inv(cmn, params)[1:5] ==
          calc_noise_weights_inv(cmn, params)[6:10]

    @test calc_noise_weights_inv(cmn, params)[1:5] ≈ [
        value(
            1 / Vela.powerlaw(
                10.0^params.TNCHROMAMP,
                params.TNCHROMGAM,
                exp(cmn.ln_js[j]) * params.PLCHROMFREQ,
                params.PLCHROMFREQ,
            ),
        ) for j = 1:5
    ]
end

@testset "PowerlawSolarWindNoiseGP" begin
    toa = default_toa()
    ctoa = TOACorrection()

    params = (
        POSEPOCH = time((53470.0 - epoch_mjd) * day_to_s),
        ELAT = dimensionless(1.2),
        ELONG = dimensionless(1.25),
        PX = GQ{-1}(3e-12),
        PMELAT = GQ{-1}(-7e-16),
        PMELONG = GQ{-1}(-5e-16),
        TNSWAMP = dimensionless(-13.5),
        TNSWGAM = dimensionless(3.5),
        PLSWFREQ = frequency(1e-8),
        PLSWEPOCH = time(day_to_s * (54000.0 - epoch_mjd)),
        PLSWSIN_ = map(dimensionless, (3.0, 2.0, 1.1, 0.12, -0.23)),
        PLSWCOS_ = map(dimensionless, (3.0, 2.0, 1.5, -1.11, -0.8)),
    )

    ss = SolarSystem(true, true)
    ctoa1 = correct_toa(ss, toa, ctoa, params)

    swn = PowerlawSolarWindNoiseGP(3, 2, 2.0)
    @test length(swn.ln_js) == 5

    @test isfinite(delay(swn, toa, ctoa1, params))

    @test all(calc_noise_weights_inv(swn, params) .> 0)
    @test calc_noise_weights_inv(swn, params)[1:5] ==
          calc_noise_weights_inv(swn, params)[6:10]

    @test calc_noise_weights_inv(swn, params)[1:5] ≈ [
        value(
            1 / Vela.powerlaw(
                10.0^params.TNSWAMP,
                params.TNSWGAM,
                exp(swn.ln_js[j]) * params.PLSWFREQ,
                params.PLSWFREQ,
            ),
        ) for j = 1:5
    ]

    @test Vela.get_gp_npars(swn) == 10
    @test length(get_marginalized_param_names(swn)) == 10

    display(swn)
end
