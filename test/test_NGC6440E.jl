@testset "NGC6440E" begin
    model, toas = Vela.load_pulsar_data("datafiles/NGC6440E.jlso")

    pepoch_mjd = value(model.epoch) / day_to_s

    @testset "model info" begin
        @test model.pulsar_name == "J1748-2021E"
        @test model.ephem == "DE440"
        @test model.clock == "TT(BIPM2021)"
        @test model.units == "TDB"
    end

    @testset "file save" begin
        Vela.save_pulsar_data("__test.jlso", model, toas)
        @test isfile("__test.jlso")
    end

    @testset "copy" begin
        m2 = TimingModel(
            model.pulsar_name,
            model.ephem,
            model.clock,
            model.units,
            model.epoch,
            model.components,
            model.kernel,
            model.param_handler,
            model.tzr_toa,
            model.priors,
        )
    end

    @testset "read_toas" begin
        @test !any([is_tzr(toa) for toa in toas])
        @test length(toas) == 62
        @test all([
            frequency(1e9) < toa.observing_frequency < frequency(2.5e9) for toa in toas
        ])
        @test all([
            time((53470.0 - pepoch_mjd) * day_to_s) <
            toa.value <
            time((54200.0 - pepoch_mjd) * day_to_s) for toa in toas
        ])
        @test all([modf(toa.pulse_number.x)[1] == 0 for toa in toas])
        @test all([toa.error > time(0.0) for toa in toas])
    end

    @testset "tzr_toa" begin
        tzrtoa = model.tzr_toa
        @test is_tzr(tzrtoa)
        @test tzrtoa.error > time(0.0)
        @test tzrtoa.pulse_number == dimensionless(0.0)
        @test frequency(1e9) < tzrtoa.observing_frequency < frequency(2.5e9)
        @test time((53470.0 - pepoch_mjd) * day_to_s) <
              tzrtoa.value <
              time((54200.0 - pepoch_mjd) * day_to_s)
    end

    @testset "param_handler" begin
        param_handler = model.param_handler
        @test Set(get_free_param_names(param_handler)) ==
              Set(["F0", "F1", "PHOFF", "RAJ", "DECJ", "DM", "EFAC1", "EQUAD1"])
        @test get_num_timing_params(model) == 6
        @test length(param_handler.multi_params) + length(param_handler.single_params) ==
              length(param_handler._default_params_tuple)
        @test length(get_free_param_names(param_handler)) ==
              length(param_handler._free_indices)
        @test sizeof(param_handler._default_params_tuple) ==
              sizeof(GQ{0,Float64}) * length(param_handler._default_values)

        @test length(get_free_param_labels(model)) == length(get_free_param_names(model))

        @test length(
            read_param_values_to_vector(
                model.param_handler,
                model.param_handler._default_params_tuple,
            ),
        ) == length(param_handler._free_indices)

        @test length(
            read_param_values_to_vector(model, model.param_handler._default_params_tuple),
        ) == length(param_handler._free_indices)

        @test all(isfinite.(get_scale_factors(model))) &&
              !any(iszero.(get_scale_factors(model)))
    end

    @testset "components" begin
        components = model.components
        @test length(components) == 5

        @test isa(components[1], SolarSystem)
        @test !components[1].ecliptic_coordinates
        @test !components[1].planet_shapiro

        # @test isa(components[2], Troposphere)

        #,@test isa(components[3], SolarWindDispersion)
        # @test components[3].model == 0

        @test isa(components[2], DispersionTaylor)

        @test isa(components[3], Spindown)

        @test isa(components[4], PhaseOffset)

        # @test all([!isa(c, Troposphere) for c in components])
    end

    @testset "repr" begin
        @test startswith(string(toas[1]), "TOA")
        display(toas)
        display(toas[1])
        # @test startswith(string(model.tzr_toa), "TZRTOA")
        display(model.tzr_toa)
        @test startswith(string(model), "TimingModel")
        display(model)
    end

    @testset "form_residuals" begin
        params = model.param_handler._default_params_tuple
        res = form_residuals(model, toas, params)
        @test all(abs(r) < 3 * toa.error for (r, toa) in zip(res, toas))
    end

    @testset "_calc_resids_and_Ndiag" begin
        params = model.param_handler._default_params_tuple
        y, Ndiag = Vela._calc_resids_and_Ndiag(model, toas, params)
        @test length(y) == length(Ndiag)
        @test all(Ndiag .> 0)
    end

    @testset "calc_chi2" begin
        calc_chi2_s = get_chi2_serial_func(model, toas)
        calc_chi2_p = get_chi2_parallel_func(model, toas)

        calc_chi2 = get_chi2_func(model, toas)

        parv = read_param_values_to_vector(model.param_handler)
        # parnp = PyArray(parv)
        chi2_s = calc_chi2_s(parv)
        chi2_p = calc_chi2_p(parv)
        @test chi2_s / degrees_of_freedom(model, toas) < 1.2
        @test chi2_s ≈ chi2_p

        params = model.param_handler._default_params_tuple
        @test @ballocated(Vela.calc_chi2_serial($model, $toas, $params)) == 0
    end

    @testset "calc_lnlike" begin
        calc_lnlike_s = get_lnlike_serial_func(model, toas)
        calc_lnlike_p = get_lnlike_parallel_func(model, toas)

        _ = get_lnlike_func(model, toas)

        params = model.param_handler._default_params_tuple
        parv = read_param_values_to_vector(model)
        # parnp = PyArray(parv)
        @test calc_lnlike_s(parv) ≈ calc_lnlike_p(parv)

        @test @ballocated(Vela.calc_lnlike_serial($model, $toas, $params)) == 0

        parv1 = read_param_values_to_vector(model.param_handler, params)
        parv1[end-2] *= 2
        @test calc_lnlike(model, toas, parv1) < calc_lnlike(model, toas, params)
    end

    @testset "priors" begin
        calc_lnprior = get_lnprior_func(model)
        params = model.param_handler._default_params_tuple
        parv = read_param_values_to_vector(model.param_handler, params)
        @test isfinite(calc_lnprior(model.param_handler._default_params_tuple))
        @test calc_lnprior(params) == calc_lnprior(parv)

        prior_transform = get_prior_transform_func(model)
        halfs = fill(0.5, length(parv))
        @test all(isfinite.(prior_transform(halfs)))
        # @test all(prior_transform(halfs) .≈ parv)
    end

    @testset "posterior" begin
        params = model.param_handler._default_params_tuple
        parv = read_param_values_to_vector(model.param_handler, params)

        calc_lnpost_ = get_lnpost_func(model, toas)
        @test isfinite(calc_lnpost_(params))

        calc_lnpost_vec = get_lnpost_func(model, toas, true)
        paramss = transpose([parv parv parv])
        @test allequal(calc_lnpost_vec(paramss))
        @test calc_lnpost_vec(paramss)[1] ≈ calc_lnpost_(params)
    end

    @testset "pulsar" begin
        psr = Pulsar(model, toas)

        params = model.param_handler._default_params_tuple
        parv = read_param_values_to_vector(model.param_handler, params)

        @test calc_lnlike(psr, parv) ≈ calc_lnlike_serial(psr, parv)

        @test isfinite(calc_lnprior(psr, parv))

        cube = 0.5 .* ones(length(parv))
        @test isfinite(calc_lnprior(psr, prior_transform(psr, cube)))

        @test calc_lnpost(psr, parv) == calc_lnpost(psr, params)
        @test calc_lnpost(psr, parv) ≈ calc_lnpost_serial(psr, parv)

        paramss = transpose([parv parv parv])
        @test allequal(calc_lnpost_vectorized(psr, paramss))
    end
end
