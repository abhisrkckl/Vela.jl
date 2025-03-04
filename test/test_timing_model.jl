@testset "TimingModel" begin
    pulsar_name = "J1234+6789"
    ephem_name = "DE440"
    clock = "TT(BIPM2021)"
    units = "TDB"
    components = (Spindown(), PhaseOffset())

    toaval = time(parse(Double64, "4610197611.8464445127"))
    freq = frequency(1.4e9)
    ephem = default_ephem()
    tzrtoa = make_tzr_toa(toaval, freq, ephem)

    pepoch = Parameter(
        :PEPOCH,
        time((56000.0 - epoch_mjd) * day_to_s),
        true,
        "day",
        float(day_to_s),
    )
    f0 = Parameter(:F0, frequency(100.0), false, "Hz", 1.0)
    f1 = Parameter(:F1, GQ{-2}(-1e-14), false, "Hz^2", 1.0)
    f_mpar = MultiParameter(:F, [f0, f1])
    phoff = Parameter(:PHOFF, dimensionless(0.0), false, "", 1.0)
    ph = ParamHandler([pepoch, phoff], [f_mpar])

    phoff_prior = SimplePrior{:PHOFF}(Uniform(-0.5, 0.5))
    f0_prior = SimplePriorMulti{:F,UInt(1)}(truncated(Normal(100.0, 1e-7); lower = 0.0))
    f1_prior = SimplePriorMulti{:F,UInt(2)}(Uniform(-1.01e-14, -0.9e-14))

    @testset "white noise kernel" begin
        wn_kernel = WhiteNoiseKernel()

        priors = (phoff_prior, f0_prior, f1_prior)
        ph = ParamHandler([pepoch, phoff], [f_mpar])

        m1 = TimingModel(
            pulsar_name,
            ephem_name,
            clock,
            units,
            time(epoch_mjd * day_to_s),
            components,
            wn_kernel,
            ph,
            tzrtoa,
            priors,
        )
        @test get_free_param_names(m1) == ["PHOFF", "F0", "F1"]
    end

    @testset "ecorr_kernel" begin
        ecorr_groups = [
            EcorrGroup(1, 8, 0),
            EcorrGroup(9, 16, 1),
            EcorrGroup(17, 24, 1),
            EcorrGroup(25, 32, 2),
            EcorrGroup(33, 44, 2),
        ]
        ecorr_kernel = EcorrKernel(ecorr_groups)

        display(ecorr_kernel)

        ecorr1 = Parameter(:ECORR1, time(1e-6), false, "us", 1e-6)
        ecorr2 = Parameter(:ECORR2, time(1e-6), false, "us", 1e-6)
        ecorr_mpar = MultiParameter(:ECORR, [ecorr1, ecorr2])
        ph = ParamHandler([pepoch, phoff], [f_mpar, ecorr_mpar])

        ecorr1_prior = SimplePriorMulti{:ECORR,UInt(1)}(Uniform(1e-9, 1e-5))
        ecorr2_prior = SimplePriorMulti{:ECORR,UInt(1)}(Uniform(1e-9, 1e-5))
        priors = (phoff_prior, f0_prior, f1_prior, ecorr1_prior, ecorr2_prior)

        m2 = TimingModel(
            pulsar_name,
            ephem_name,
            clock,
            units,
            time(epoch_mjd * day_to_s),
            components,
            ecorr_kernel,
            ph,
            tzrtoa,
            priors,
        )

        @test get_free_param_names(m2) == ["PHOFF", "F0", "F1", "ECORR1", "ECORR2"]
        @test length(get_free_param_units(m2)) == length(get_free_param_names(m2))
        @test length(get_free_param_prefixes(m2)) == length(get_free_param_names(m2))
    end
end
