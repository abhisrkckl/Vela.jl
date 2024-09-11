@testset "TimingModel" begin
    @testset "WhiteNoiseKernel" begin
        pulsar_name = "J1234+6789"
        ephem_name = "DE440"
        clock = "TT(BIPM2021)"
        units = "TDB"
        components = (Spindown(), PhaseOffset())
        kernel = WhiteNoiseKernel()


        pepoch = Parameter(:PEPOCH, time(56000.0 * day_to_s), true, "day", float(day_to_s))
        f0 = Parameter(:F0, frequency(100.0), false, "Hz", 1.0)
        f1 = Parameter(:F1, GQ{-2}(-1e-14), false, "Hz^2", 1.0)
        f_mpar = MultiParameter(:F, [f0, f1])
        phoff = Parameter(:PHOFF, dimensionless(0.0), false, "", 1.0)
        ph = ParamHandler([pepoch, phoff], [f_mpar])

        toaval = time(parse(Double64, "4610197611.8464445127"))
        freq = frequency(1.4e9)
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
        tzrtoa = make_tzr_toa(toaval, freq, true, ephem)

        phoff_prior = SimplePrior{:PHOFF}(Uniform(-0.5, 0.5))
        f0_prior = SimplePriorMulti{:F,UInt(1)}(truncated(Normal(100.0, 1e-7); lower = 0.0))
        f1_prior = SimplePriorMulti{:F,UInt(2)}(Uniform(-1.01e-14, -0.9e-14))
        priors = (phoff_prior, f0_prior, f1_prior)

        m = TimingModel(
            pulsar_name,
            ephem_name,
            clock,
            units,
            components,
            kernel,
            ph,
            tzrtoa,
            priors,
        )

        @test get_free_param_names(m) == ["PHOFF", "F0", "F1"]
    end
end
