using Vela
using Test
using GeometricUnits
using LinearAlgebra
using Quadmath
using JuliaFormatter
using HDF5
using BenchmarkTools

const day_to_s = 86400

@testset "Vela" verbose = true begin

    ssb_obs_pos = distance.((18.0354099, 450.01472245, 195.05827732))
    ssb_obs_vel = speed.((-9.96231954e-05, 3.31555854e-06, 1.12968547e-06))
    obs_sun_pos = distance.((-15.89483533, -450.17185232, -195.18212616))
    obs_jupiter_pos = distance.((-1610.1796849, -1706.87348483, -681.22381513))
    obs_saturn_pos = distance.((-2392.85651431, 3109.13626083, 1405.71912274))
    obs_venus_pos = distance.((140.85922773, -217.65571843, -74.64804201))
    obs_uranus_pos = distance.((9936.62957939, -3089.07377113, -1486.17339104))
    obs_neptune_pos = distance.((11518.60924426, -9405.0693235, -4126.91030657))
    obs_earth_pos = distance.((0.01199435, 0.01159591, -0.01316261))

    @testset "ephemeris" begin
        # Wrong dimensions
        @test_throws AssertionError SolarSystemEphemeris(
            ssb_obs_vel,
            ssb_obs_vel,
            obs_sun_pos,
            obs_jupiter_pos,
            obs_saturn_pos,
            obs_venus_pos,
            obs_uranus_pos,
            obs_neptune_pos,
            obs_earth_pos,
        )

        # Wrong dimensions
        @test_throws AssertionError SolarSystemEphemeris(
            ssb_obs_pos,
            ssb_obs_pos,
            obs_sun_pos,
            obs_jupiter_pos,
            obs_saturn_pos,
            obs_venus_pos,
            obs_uranus_pos,
            obs_neptune_pos,
            obs_earth_pos,
        )

        # Wrong dimensions
        @test_throws AssertionError SolarSystemEphemeris(
            ssb_obs_pos,
            ssb_obs_vel,
            ssb_obs_vel,
            obs_jupiter_pos,
            obs_saturn_pos,
            obs_venus_pos,
            obs_uranus_pos,
            obs_neptune_pos,
            obs_earth_pos,
        )

        # ssb_obs_vel is too large.
        @test_throws AssertionError SolarSystemEphemeris(
            ssb_obs_pos,
            1e6 .* ssb_obs_vel,
            obs_sun_pos,
            obs_jupiter_pos,
            obs_saturn_pos,
            obs_venus_pos,
            obs_uranus_pos,
            obs_neptune_pos,
            obs_earth_pos,
        )

        ephem_vecs = SolarSystemEphemeris(
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

        # ssb_obs_pos and obs_sun_pos should be less than the Aphelion distance
        @test sqrt(dot(ephem_vecs.ssb_obs_pos, ephem_vecs.ssb_obs_pos)) < distance(509.0)
        @test sqrt(dot(ephem_vecs.obs_sun_pos, ephem_vecs.obs_sun_pos)) < distance(509.0)

        # ssb_obs_pos and obs_sun_pos should have a small angle.
        @test acos(
            -dot(ephem_vecs.ssb_obs_pos, ephem_vecs.obs_sun_pos) / sqrt(
                dot(ephem_vecs.ssb_obs_pos, ephem_vecs.ssb_obs_pos) *
                dot(ephem_vecs.obs_sun_pos, ephem_vecs.obs_sun_pos),
            ),
        ) < 0.01

        # |ssb_obs_vel| should be less than the speed of light
        @test dot(ephem_vecs.ssb_obs_vel, ephem_vecs.ssb_obs_vel) < 1

    end

    @testset "toa" begin
        toaval = time(parse(Float128, "4610197611.8464445127"))
        toaerr = time(1e-6)
        freq = frequency(1.4e9)
        pulse_number = dimensionless(Float128(1000.0))
        barycentered = false

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

        # TOA value should be of type GQ{Float128}.
        @test_throws MethodError TOA(
            time(4610197611.8),
            toaerr,
            freq,
            pulse_number,
            barycentered,
            ephem,
        )

        # Wrong dimensions for TOA value.
        @test_throws AssertionError TOA(
            dimensionless(toaval.x),
            toaerr,
            freq,
            pulse_number,
            barycentered,
            ephem,
        )

        # Wrong dimensions for TOA error.
        @test_throws AssertionError TOA(
            toaval,
            dimensionless(1e-6),
            freq,
            pulse_number,
            barycentered,
            ephem,
        )

        # Wrong dimensions for TOA observing_frequency.
        @test_throws AssertionError TOA(
            toaval,
            toaerr,
            time(1.4e9),
            phase,
            barycentered,
            ephem,
        )

        # Wrong dimensions for TOA pulse_number.
        @test_throws AssertionError TOA(
            toaval,
            toaerr,
            freq,
            time(1000.0),
            barycentered,
            ephem,
        )

        toa1 = TOA(toaval, toaerr, freq, pulse_number, barycentered, ephem)
        @test !toa1.tzr

        ctoa1 = CorrectedTOA(toa1)
        @test ctoa1.level == 0

        dt = time(1.0)
        ctoa2 = correct_toa(ctoa1; delay = dt)
        @test ctoa2.delay == ctoa1.delay + dt
        @test corrected_toa_value(ctoa2) == corrected_toa_value(ctoa1) - dt
        @test ctoa2.phase == ctoa1.phase
        @test ctoa2.efac == ctoa1.efac
        @test ctoa2.equad2 == ctoa1.equad2
        @test ctoa2.doppler == ctoa1.doppler
        @test ctoa2.barycentered == ctoa1.barycentered
        @test ctoa2.level == 1

        dphi = dimensionless(0.3)
        ctoa3 = correct_toa(ctoa2; phase = dphi)
        @test ctoa3.delay == ctoa2.delay
        @test ctoa3.phase == ctoa2.phase + dphi
        @test phase_residual(ctoa3) == phase_residual(ctoa2) + dphi
        @test ctoa3.efac == ctoa2.efac
        @test ctoa3.equad2 == ctoa2.equad2
        @test ctoa3.doppler == ctoa2.doppler
        @test ctoa3.barycentered == ctoa2.barycentered
        @test ctoa3.level == 2

        efac = dimensionless(1.1)
        equad2 = time(1e-6)^2
        ctoa4 = correct_toa(ctoa3; efac = efac, equad2 = equad2)
        @test ctoa4.delay == ctoa3.delay
        @test ctoa4.phase == ctoa3.phase
        @test ctoa4.efac == ctoa3.efac * efac
        @test ctoa4.equad2 == ctoa3.equad2 + equad2
        @test scaled_toa_error_sqr(ctoa4) ≈ (scaled_toa_error_sqr(ctoa3) + equad2) * efac^2
        @test ctoa4.doppler == ctoa3.doppler
        @test ctoa4.barycentered == ctoa3.barycentered
        @test ctoa4.level == 3

        @testset "tzr_toa" begin
            tzrtoa = make_tzr_toa(toaval, freq, true, ephem)
            @test tzrtoa.tzr
            @test tzrtoa.barycentered
            @test tzrtoa.error == time(0.0)

            ctzrtoa = CorrectedTOA(tzrtoa)
            @test ctzrtoa.level == 0
        end
    end

    @testset "parameter & param handler" begin
        pepoch = Parameter(:PEPOCH, time(56000.0 * day_to_s), true, "day", day_to_s)
        @test pepoch.name == :PEPOCH
        @test_throws AssertionError read_param(pepoch, 56000.0 * day_to_s)

        f0 = Parameter(:F0, frequency(100.0), false, "Hz", 1.0)
        @test read_param(f0, 101.0) == frequency(101.0)

        f1 = Parameter(:F1, GQ(-1e-14, -2), false, "Hz^2", 1.0)
        @test read_param(f1, -1.1e-14) == GQ(-1.1e-14, -2)

        mparF = MultiParameter(:F, [f0, f1])
        @test length(mparF.parameters) == 2

        ph = ParamHandler([pepoch], [mparF])
        @test get_free_param_names(ph) == ["F0", "F1"]
        @test Set(keys(ph._default_params_tuple)) == Set([:PEPOCH, :F])

        params = read_params(ph, [100.01, -1.01e-14])
        @test Set(keys(params)) == Set([:PEPOCH, :F])
        @test params.F == (frequency(100.01), GQ(-1.01e-14, -2))
    end

    @testset "components" begin
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

        toa = TOA(
            time(Float128(53470.0 * day_to_s)),
            time(1e-6),
            frequency(2.5e9),
            dimensionless(Float128(0.0)),
            false,
            ephem,
        )
        ctoa = CorrectedTOA(toa)

        tzrtoa =
            make_tzr_toa(time(Float128(53475.0 * day_to_s)), frequency(2.5e9), false, ephem)
        ctzrtoa = CorrectedTOA(tzrtoa)

        params = (
            PHOFF = dimensionless(1e-6),
            PEPOCH = time(53470.0 * day_to_s),
            F = (frequency(100.0), GQ(-1e-14, -2)),
            DMEPOCH = time(53470.0 * day_to_s),
            DM = (GQ(10.0, -1), GQ(1e-4, -2)),
            POSEPOCH = time(53470.0 * day_to_s),
            ELAT = dimensionless(1.2),
            ELONG = dimensionless(1.25),
            PX = GQ(3e-12, -1),
            PMELAT = GQ(-7e-16, -1),
            PMELONG = GQ(-5e-16, -1),
        )

        @testset "SolarSystem" begin
            ss = SolarSystem(true, true)
            @test ss.ecliptic_coordinates && ss.planet_shapiro

            ctoa1 = correct_toa(ss, ctoa, params)

            @test ctoa1.phase == ctoa.phase
            @test ctoa1.delay != ctoa.delay
            @test ctoa.doppler == 0 && ctoa1.doppler != 0
            @test !ctoa.barycentered && ctoa1.barycentered
            @test ctoa1.level == ctoa.level + 1

            ctoa2 = correct_toa(ss, ctoa1, params)
            @test (ctoa2.delay == ctoa1.delay) && (ctoa2.doppler == ctoa1.doppler)

            @test isfinite(lnprior(ss, params))
            @test lnprior(
                ss,
                (;
                    ELAT = dimensionless(1.2),
                    ELONG = dimensionless(1.25),
                    PX = GQ(-3e-12, -1),
                    PMELAT = GQ(-7e-16, -1),
                    PMELONG = GQ(-5e-16, -1),
                ),
            ) == -Inf
            @test @ballocated(lnprior($ss, $params)) == 0
        end

        @testset "PhaseOffset" begin
            poff = PhaseOffset()
            @test phase(poff, ctoa, params) == dimensionless(-1e-6)
            @test phase(poff, ctzrtoa, params) == dimensionless(0.0)

            ctoa1 = correct_toa(poff, ctoa, params)
            @test ctoa1.delay == ctoa.delay
            @test ctoa1.phase ≈ ctoa.phase + phase(poff, ctoa, params)
            @test ctoa1.efac == ctoa.efac
            @test ctoa1.equad2 == ctoa.equad2
            @test ctoa1.spin_frequency == ctoa.spin_frequency == frequency(-1.0)

            ctzrtoa1 = correct_toa(poff, ctzrtoa, params)
            @test ctzrtoa1.delay == ctzrtoa.delay
            @test ctzrtoa1.phase == ctzrtoa.phase
            @test ctzrtoa1.efac == ctzrtoa.efac
            @test ctzrtoa1.equad2 == ctzrtoa.equad2

            @test lnprior(poff, (; PHOFF = dimensionless(0.1))) == 0.0
            @test lnprior(poff, (; PHOFF = dimensionless(0.6))) == -Inf
            @test @ballocated(lnprior($poff, (; PHOFF = dimensionless(0.1)))) == 0
        end

        @testset "Spindown" begin
            spn = Spindown()
            @test phase(spn, ctoa, params) == dimensionless(0.0)
            @test spin_frequency(spn, ctoa, params) == frequency(100.0)

            ctoa1 = correct_toa(spn, ctoa, params)
            @test ctoa1.delay == ctoa.delay
            @test ctoa1.phase == ctoa.phase + phase(spn, ctoa, params)
            @test ctoa1.doppler == ctoa.doppler
            @test ctoa.spin_frequency == frequency(-1.0) &&
                  ctoa1.spin_frequency > frequency(0.0)

            @test isfinite(lnprior(spn, (; F = (frequency(101.1), GQ(-1e-14, -2)))))
            @test lnprior(spn, (; F = (frequency(2001.1), GQ(-1e-14, -2)))) == -Inf
            @test lnprior(spn, (; F = (frequency(200.1), GQ(-1e-1, -2)))) == -Inf
            @test @ballocated(lnprior($spn, (; F = (frequency(101.1), GQ(-1e-14, -2))))) ==
                  0
        end

        # @testset "Troposphere" begin
        #     tropo = Troposphere()
        #     @test delay(tropo, toa, params) == time(0.0)
        # end

        @testset "DispersionTaylor" begin
            dmt = DispersionTaylor()
            @test dispersion_slope(dmt, ctoa, params) == GQ(10.0, -1)
            @test delay(dmt, ctoa, params) ==
                  dispersion_slope(dmt, ctoa, params) / ctoa.toa.observing_frequency^2

            @test isfinite(lnprior(dmt, (; DM = (GQ(6e+16, -1), GQ(-4e+11, -2)))))
            @test lnprior(dmt, (; DM = (GQ(-6e+16, -1), GQ(-4e+11, -2)))) == -Inf
            @test lnprior(dmt, (; DM = (GQ(6e+16, -1), GQ(-4e+20, -2)))) == -Inf
            @test @ballocated(lnprior($dmt, (; DM = (GQ(6e+16, -2), GQ(-4e+11, -2))))) == 0
        end

        @testset "SolarWindDispersion" begin
            @test_throws AssertionError SolarWindDispersion(2)

            swd = SolarWindDispersion(0)
            @test dispersion_slope(swd, toa, params) == GQ(0.0, -1)
        end
    end

    @testset "pure_rotator" begin
        model, toas = read_model_and_toas("pure_rotator.hdf5")

        @testset "model info" begin
            @test model.pulsar_name == "SIM0"
            @test model.ephem == "DE421"
            @test model.clock == "TT(BIPM2021)"
            @test model.units == "TDB"
        end

        @testset "read_toas" begin
            @test !any([toa.tzr for toa in toas])
            @test length(toas) == 100

            @test all([toa.observing_frequency == frequency(1.4e9) for toa in toas])
            @test all([
                time(53400.0 * day_to_s) < toa.value < time(56002.0 * day_to_s) for
                toa in toas
            ])
            @test all([modf(toa.pulse_number.x)[1] == 0 for toa in toas])
            @test all([toa.error > time(0.0) for toa in toas])
        end

        @testset "read_tzr_toa" begin
            tzrtoa = model.tzr_toa
            @test tzrtoa.tzr
            @test tzrtoa.error > time(0.0)
            @test modf(tzrtoa.pulse_number.x)[1] == 0
            @test frequency(1e9) < tzrtoa.observing_frequency < frequency(2.5e9)
            @test time(53400.0 * day_to_s) < tzrtoa.value < time(56002.0 * day_to_s)
        end

        @testset "read_param_handler" begin
            param_handler = model.param_handler
            @test Set(get_free_param_names(param_handler)) == Set(["F0", "F1", "PHOFF"])
            @test length(param_handler.multi_params) +
                  length(param_handler.single_params) ==
                  length(param_handler._default_params_tuple)
            @test length(get_free_param_names(param_handler)) ==
                  length(param_handler._free_indices)
            @test sizeof(param_handler._default_params_tuple) ==
                  sizeof(GQ{Float64}) * length(param_handler._default_quantities)

            params = model.param_handler._default_params_tuple
            parv = [params.PHOFF.x, params.F[1].x, params.F[2].x]
            @test read_param_values_to_vector(model.param_handler, params) == parv
        end

        @testset "read_components" begin
            components = model.components
            @test length(components) == 2
            @test isa(components[1], Spindown)
            @test isa(components[2], PhaseOffset)
        end

        @testset "repr" begin
            @test startswith(string(toas[1]), "TOA")
            display(toas)
            display(toas[1])
            @test startswith(string(model.tzr_toa), "TZRTOA")
            display(model.tzr_toa)
            @test startswith(string(model), "TimingModel")
            display(model)
        end

        params = model.param_handler._default_params_tuple
        parv = read_param_values_to_vector(model.param_handler, params)

        @testset "correct_toa" begin
            ctoa = correct_toa(model, toas[1], params)
            @test corrected_toa_value(ctoa) ≈ toas[1].value
            @test doppler_shifted_spin_frequency(ctoa) ≈ ctoa.spin_frequency
        end

        @testset "form_residual" begin
            resid = form_residual(model, toas[1], params, dimensionless(0.0))
            @test abs(resid) < time(1e-2)
        end

        @testset "calc_chi2" begin
            chi2 = calc_chi2(model, toas, params)
            nfree = length(get_free_param_names(model.param_handler))
            @test isa(chi2, Float64)
            @test chi2 / (length(toas) - nfree) < dimensionless(1.5)

            @test chi2 ≈ calc_chi2(model, toas, parv)
            @test chi2 ≈ Vela.calc_chi2_serial(model, toas, params)
            @test chi2 ≈ Vela.calc_chi2_serial(model, toas, parv)

            @test @ballocated(Vela.calc_chi2_serial($model, $toas, $params)) == 0
        end

        @testset "calc_lnlike" begin
            lnlike_func = get_lnlike_func(model, toas)
            lnlike_serial_func = Vela.get_lnlike_serial_func(model, toas)

            # lnlike = calc_lnlike(model, toas, params)
            lnlike = lnlike_func(params)
            @test isa(lnlike, Float64)
            @test isfinite(lnlike)

            @test lnlike ≈ lnlike_func(parv)
            @test lnlike ≈ lnlike_serial_func(params)
            @test lnlike ≈ lnlike_serial_func(parv)

            @test @ballocated(Vela.calc_lnlike_serial($model, $toas, $params)) == 0
        end

        @testset "lnprior" begin
            @test_broken isfinite(lnprior(model, model.param_handler._default_params_tuple))
        end
    end

    @testset "NGC6440E" begin
        model, toas = read_model_and_toas("NGC6440E.hdf5")

        @testset "read_toas" begin
            @test !any([toa.tzr for toa in toas])
            @test length(toas) == 62
            @test all([
                frequency(1e9) < toa.observing_frequency < frequency(2.5e9) for toa in toas
            ])
            @test all([
                time(53470.0 * day_to_s) < toa.value < time(54200.0 * day_to_s) for
                toa in toas
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
            @test length(param_handler.multi_params) +
                  length(param_handler.single_params) ==
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
            tzrphase = Vela.calc_tzr_phase(model, params)
            res = form_residual(model, toas[1], params, tzrphase)
            @test abs(res) < 3 * toas[1].error
        end

        @testset "calc_chi2" begin
            params = model.param_handler._default_params_tuple
            chi2 = calc_chi2(model, toas, params)
            @test chi2 / length(toas) < 1.1
            @test chi2 == Vela.calc_chi2_serial(model, toas, params)
        end

        @testset "calc_lnlike" begin
            params = model.param_handler._default_params_tuple
            @test calc_lnlike(model, toas, params) ==
                  Vela.calc_lnlike_serial(model, toas, params)
            @test @ballocated(Vela.calc_lnlike_serial($model, $toas, $params)) == 0

            parv1 = read_param_values_to_vector(model.param_handler, params)
            parv1[end] *= 2
            @test calc_lnlike(model, toas, parv1) < calc_lnlike(model, toas, params)
        end

        @testset "calc_lnprior" begin
            @test isfinite(calc_lnprior(model, model.param_handler._default_params_tuple))

            # This doesn't allocate when model is a const.
            @test_broken @ballocated(
                calc_lnprior($model, $model.param_handler._default_params_tuple)
            ) == 0
        end

        @testset "calc_lnpost" begin
            @test isfinite(
                calc_lnpost(model, toas, model.param_handler._default_params_tuple),
            )

            @test_broken @ballocated(
                calc_lnpost($model, $toas, $model.param_handler._default_params_tuple)
            ) == 0
        end

        @testset "plot_summary" begin
            plotfile = plot_pulsar_summary("NGC6440E.hdf5")
            @test isfile(plotfile)
        end
    end

    @testset "formatting" begin
        @test format(Vela)
    end
end
