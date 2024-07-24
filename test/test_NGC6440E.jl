@testset "NGC6440E" begin
    model, toas = read_pulsar("datafiles/NGC6440E.jlso")

    @testset "read_toas" begin
        @test !any([toa.tzr for toa in toas])
        @test length(toas) == 62
        @test all([
            frequency(1e9) < toa.observing_frequency < frequency(2.5e9) for toa in toas
        ])
        @test all([
            time(53470.0 * day_to_s) < toa.value < time(54200.0 * day_to_s) for toa in toas
        ])
        @test all([modf(toa.pulse_number.x)[1] == 0 for toa in toas])
        @test all([toa.error > time(0.0) for toa in toas])

    end

    @testset "read_tzr_toa" begin
        tzrtoa = model.tzr_toa
        @test tzrtoa.tzr
        @test tzrtoa.error > time(0.0)
        @test tzrtoa.pulse_number == dimensionless(0.0)
        @test frequency(1e9) < tzrtoa.observing_frequency < frequency(2.5e9)
        @test time(53470.0 * day_to_s) < tzrtoa.value < time(54200.0 * day_to_s)
    end

    @testset "read_param_handler" begin
        param_handler = model.param_handler
        @test Set(get_free_param_names(param_handler)) ==
              Set(["F0", "F1", "PHOFF", "RAJ", "DECJ", "DM"])
        @test length(param_handler.multi_params) + length(param_handler.single_params) ==
              length(param_handler._default_params_tuple)
        @test length(get_free_param_names(param_handler)) ==
              length(param_handler._free_indices)
        @test sizeof(param_handler._default_params_tuple) ==
              sizeof(GQ{Float64}) * length(param_handler._default_quantities)

        @test length(
            read_param_values_to_vector(
                model.param_handler,
                model.param_handler._default_params_tuple,
            ),
        ) == length(param_handler._free_indices)
    end

    @testset "read_components" begin
        components = model.components
        @test length(components) == 4

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

    @testset "form_residual" begin
        params = model.param_handler._default_params_tuple
        tzrphase = calc_tzr_phase(model, params)
        res = form_residual(model, toas[1], params, tzrphase)
        @test abs(res) < 3 * toas[1].error
    end

    @testset "calc_chi2" begin
        calc_chi2_s = get_chi2_serial_func(model, toas)
        calc_chi2_p = get_chi2_parallel_func(model, toas)
        params = model.param_handler._default_params_tuple
        parv = read_param_values_to_vector(model.param_handler, params)
        # parnp = PyArray(parv)
        chi2_s = calc_chi2_s(parv)
        chi2_p = calc_chi2_p(parv)
        @test chi2_s / length(toas) < 1.1
        @test chi2_s ≈ chi2_p
    end

    @testset "calc_lnlike" begin
        calc_lnlike_s = get_lnlike_serial_func(model, toas)
        calc_lnlike_p = get_lnlike_parallel_func(model, toas)
        params = model.param_handler._default_params_tuple
        parv = read_param_values_to_vector(model.param_handler, params)
        # parnp = PyArray(parv)
        @test calc_lnlike_s(parv) ≈ calc_lnlike_p(parv)

        @test @ballocated(Vela.calc_lnlike_serial($model, $toas, $params)) == 0

        parv1 = read_param_values_to_vector(model.param_handler, params)
        parv1[end] *= 2
        @test calc_lnlike(model, toas, parv1) < calc_lnlike(model, toas, params)
    end

    # @testset "plot_summary" begin
    #     plotfile = plot_pulsar_summary("datafiles/NGC6440E.hdf5")
    #     @test isfile(plotfile)
    # end
end
