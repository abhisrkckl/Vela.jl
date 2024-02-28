using Vela
using Test
using GeometricUnits
using LinearAlgebra
using Quadmath
using JuliaFormatter
using HDF5

const day_to_s = 86400

@testset "Vela" verbose = true begin
    ssb_obs_pos = distance.([18.0354099, 450.01472245, 195.05827732])
    ssb_obs_vel = speed.([-9.96231954e-05, 3.31555854e-06, 1.12968547e-06])
    obs_sun_pos = distance.([-15.89483533, -450.17185232, -195.18212616])
    obs_jupiter_pos = distance.([-1610.1796849, -1706.87348483, -681.22381513])
    obs_saturn_pos = distance.([-2392.85651431, 3109.13626083, 1405.71912274])
    obs_venus_pos = distance.([140.85922773, -217.65571843, -74.64804201])
    obs_uranus_pos = distance.([9936.62957939, -3089.07377113, -1486.17339104])
    obs_neptune_pos = distance.([11518.60924426, -9405.0693235, -4126.91030657])
    obs_earth_pos = distance.([0.01199435, 0.01159591, -0.01316261])

    @testset "toa" begin
        @test_throws AssertionError EphemerisVectors(
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
        @test_throws AssertionError EphemerisVectors(
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
        @test_throws AssertionError EphemerisVectors(
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
        @test_throws AssertionError EphemerisVectors(
            ssb_obs_pos,
            1e6 * ssb_obs_vel,
            obs_sun_pos,
            obs_jupiter_pos,
            obs_saturn_pos,
            obs_venus_pos,
            obs_uranus_pos,
            obs_neptune_pos,
            obs_earth_pos,
        )

        ephem_vecs = EphemerisVectors(
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

        # Aphelion distance
        @test sqrt(dot(ephem_vecs.ssb_obs_pos, ephem_vecs.ssb_obs_pos)) < distance(509.0)
        @test sqrt(dot(ephem_vecs.obs_sun_pos, ephem_vecs.obs_sun_pos)) < distance(509.0)

        # Angle
        @test acos(
            -dot(ephem_vecs.ssb_obs_pos, ephem_vecs.obs_sun_pos) / sqrt(
                dot(ephem_vecs.ssb_obs_pos, ephem_vecs.ssb_obs_pos) *
                dot(ephem_vecs.obs_sun_pos, ephem_vecs.obs_sun_pos),
            ),
        ) < 0.01

        # Speed of light
        @test dot(ephem_vecs.ssb_obs_vel, ephem_vecs.ssb_obs_vel) < 1

        toaval = time(parse(Float128, "4610197611.8464445127"))
        toaerr = time(1e-6)
        freq = frequency(1.4e9)
        phase = dimensionless(Float128(1000.0))

        @test_throws MethodError TOA(time(4610197611.8), toaerr, freq, phase, ephem_vecs)
        @test_throws AssertionError TOA(
            dimensionless(parse(Float128, "4610197611.8464445127")),
            toaerr,
            freq,
            phase,
            ephem_vecs,
        )
        @test_throws AssertionError TOA(
            toaval,
            dimensionless(1e-6),
            freq,
            phase,
            ephem_vecs,
        )
        @test_throws AssertionError TOA(toaval, toaerr, time(1.4e9), phase, ephem_vecs)
        @test_throws MethodError TOA(
            toaval,
            toaerr,
            freq,
            dimensionless(1000.0),
            ephem_vecs,
        )
        @test_throws AssertionError TOA(
            toaval,
            toaerr,
            freq,
            time(Float128(1000.0)),
            ephem_vecs,
        )

        # toa1 = TOA(toaval, toaerr, freq, phase)
        # @test is_barycentered(toa1)
        # @test !toa1.tzr
        # @test toa1.level == 0

        toa1 = TOA(toaval, toaerr, freq, phase, ephem_vecs)
        # @test !is_barycentered(toa1)
        @test !toa1.tzr
        @test toa1.level == 0

        dt = time(1.0)
        toa3 = deepcopy(toa1)
        correct_toa_delay!(toa3, dt)
        @test toa3.value == toa1.value - dt
        @test toa3.error == toa1.error
        @test toa3.observing_frequency == toa1.observing_frequency
        @test toa3.phase == toa1.phase
        # @test is_barycentered(toa3) == is_barycentered(toa1)
        @test toa3.tzr == toa1.tzr
        @test toa3.level == 1

        dphi = dimensionless(0.3)
        toa4 = deepcopy(toa3)
        correct_toa_phase!(toa4, dphi)
        @test toa4.value == toa3.value
        @test toa4.error == toa3.error
        @test toa4.observing_frequency == toa3.observing_frequency
        @test toa4.phase == toa3.phase + dphi
        # @test is_barycentered(toa4) == is_barycentered(toa3)
        @test toa4.tzr == toa3.tzr
        @test toa4.level == 2

        @testset "tzr_toa" begin
            tzrtoa = make_tzr_toa(toaval, freq, ephem_vecs)
            @test tzrtoa.tzr
            # @test !is_barycentered(tzrtoa)
            @test tzrtoa.error == time(0.0)
            @test tzrtoa.level == 0
        end
    end

    @testset "parameter & param handler" begin
        pepoch = Parameter("pepoch", time(56000.0 * day_to_s), true)
        @test pepoch.display_name == "PEPOCH"
        @test_throws AssertionError read_param(pepoch, 56000.0 * day_to_s)

        f0 = Parameter("F0", frequency(100.0), false)
        @test read_param(f0, 101.0) == frequency(101.0)

        f1 = Parameter("F1", GQ(-1e-14, -2), false)
        @test read_param(f1, -1.1e-14) == GQ(-1.1e-14, -2)

        mparT = MultiParameter("PEPOCH", [pepoch])
        @test length(mparT.parameters) == 1

        mparF = MultiParameter("F", [f0, f1])
        @test length(mparF.parameters) == 2

        ph = ParamHandler([mparT, mparF])
        @test get_free_param_names(ph) == ["F0", "F1"]
        @test keys(ph._default_params_dict) == Set(["PEPOCH", "F"])

        params_dict = read_params(ph, [100.01, -1.01e-14])
        @test keys(params_dict) == Set(["PEPOCH", "F"])
        @test params_dict["F"] == [frequency(100.01), GQ(-1.01e-14, -2)]
    end

    @testset "components" begin
        ephem_vecs = EphemerisVectors(
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
            ephem_vecs,
        )
        tzrtoa =
            make_tzr_toa(time(Float128(53475.0 * day_to_s)), frequency(2.5e9), ephem_vecs)

        params = Dict(
            "PHOFF" => [dimensionless(1e-6)],
            "PEPOCH" => [time(53470.0 * day_to_s)],
            "F" => [frequency(100.0), GQ(-1e-14, -2)],
            "DMEPOCH" => [time(53470.0 * day_to_s)],
            "DM" => [GQ(10.0, -1), GQ(1e-4, -2)],
        )

        @testset "PhaseOffset" begin
            poff = PhaseOffset()
            poff_params = read_params_from_dict(poff, params)
            @test phase(poff, toa, poff_params) == dimensionless(-1e-6)
            @test phase(poff, tzrtoa, poff_params) == dimensionless(0.0)

            ctoas = [deepcopy(toa)]
            correct_toas!(poff, ctoas, params)
            @test ctoas[1].value == toa.value
            @test ctoas[1].error == toa.error
            @test ctoas[1].observing_frequency == toa.observing_frequency
            @test toa.spin_frequency == frequency(-1.0) &&
                  ctoas[1].spin_frequency == frequency(-1.0)
            @test ctoas[1].phase ≈ toa.phase + phase(poff, toa, poff_params)

            ctzrtoa_ = [deepcopy(tzrtoa)]
            correct_toas!(poff, ctzrtoa_, params)
            ctzrtoa = ctzrtoa_[1]
            @test ctzrtoa.value == tzrtoa.value
            @test ctzrtoa.error == tzrtoa.error
            @test ctzrtoa.observing_frequency == tzrtoa.observing_frequency
            @test toa.spin_frequency == frequency(-1.0) &&
                  ctzrtoa.spin_frequency == frequency(-1.0)
            @test ctzrtoa.phase == tzrtoa.phase
        end

        @testset "Spindown" begin
            spn = Spindown()
            spn_params = read_params_from_dict(spn, params)
            @test phase(spn, toa, spn_params) == dimensionless(0.0)
            @test spin_frequency(spn, toa, spn_params) == frequency(100.0)

            ctoas = [deepcopy(toa)]
            correct_toas!(spn, ctoas, params)
            ctoa = ctoas[1]
            @test ctoa.value == toa.value
            @test ctoa.error == toa.error
            @test ctoa.observing_frequency == toa.observing_frequency
            @test toa.spin_frequency == frequency(-1.0) &&
                  ctoa.spin_frequency > frequency(0.0)
            @test ctoa.phase == toa.phase + phase(spn, toa, spn_params)
        end

        @testset "SolarSystem" begin
            ss = SolarSystem(true, true, true)
            @test ss.ecliptic_coordinates && ss.planet_shapiro && ss.proper_motion
            @test delay(ss, toa, params) == time(0.0)
        end

        @testset "Troposphere" begin
            tropo = Troposphere()
            @test delay(tropo, toa, params) == time(0.0)
        end

        @testset "DispersionTaylor" begin
            dmt = DispersionTaylor()
            @test dispersion_slope(dmt, toa, params) == GQ(10.0, -1)
            @test delay(dmt, toa, params) ==
                  dispersion_slope(dmt, toa, params) / toa.observing_frequency^2
        end

        @testset "SolarWindDispersion" begin
            @test_throws AssertionError SolarWindDispersion(2)

            swd = SolarWindDispersion(0)
            @test dispersion_slope(swd, toa, params) == GQ(0.0, -1)
        end
    end

    @testset "Test datasets" begin
        @testset "NGC6440E" begin
            model, toas = read_model_and_toas("NGC6440E.hdf5")

            @testset "read_toas" begin
                @test !any([toa.tzr for toa in toas])
                @test length(toas) == 62
                @test all([toa.level == 0 for toa in toas])
                @test all([
                    frequency(1e9) < toa.observing_frequency < frequency(2.5e9) for
                    toa in toas
                ])
                @test all([
                    time(53470.0 * day_to_s) < toa.value < time(54200.0 * day_to_s) for
                    toa in toas
                ])
                @test all([modf(toa.phase.x)[1] == 0 for toa in toas])
                @test all([toa.error > time(0.0) for toa in toas])

            end

            @testset "read_tzr_toa" begin
                tzrtoa = model.tzr_toa
                @test tzrtoa.tzr
                @test tzrtoa.level == 0
                @test tzrtoa.error > time(0.0)
                @test modf(tzrtoa.phase.x)[1] == 0
                @test frequency(1e9) < tzrtoa.observing_frequency < frequency(2.5e9)
                @test time(53470.0 * day_to_s) < tzrtoa.value < time(54200.0 * day_to_s)
            end

            @testset "read_param_handler" begin
                param_handler = model.param_handler
                @test length(param_handler.multi_params) ==
                      length(param_handler._default_params_dict)
                @test Set(get_free_param_names(param_handler)) ==
                      Set(["F0", "F1", "DECJ", "RAJ", "DM", "PHOFF"])
            end

            @testset "read_components" begin
                components = model.components
                @test length(components) == 6

                @test isa(components[1], SolarSystem)
                @test !components[1].ecliptic_coordinates
                @test !components[1].planet_shapiro
                @test !components[1].proper_motion

                @test isa(components[2], Troposphere)

                @test isa(components[3], SolarWindDispersion)
                @test components[3].model == 0

                @test isa(components[4], DispersionTaylor)

                @test isa(components[5], Spindown)

                @test isa(components[6], PhaseOffset)

                # @test all([!isa(c, Troposphere) for c in components])
            end
        end

        @testset "pure_rotator" begin
            model, toas = read_model_and_toas("pure_rotator.hdf5")

            @testset "read_toas" begin
                @test !any([toa.tzr for toa in toas])
                @test length(toas) == 100
                @test all([toa.level == 0 for toa in toas])
                @test all([toa.observing_frequency == frequency(1.4e9) for toa in toas])
                @test all([
                    time(53400.0 * day_to_s) < toa.value < time(56002.0 * day_to_s) for
                    toa in toas
                ])
                @test all([modf(toa.phase.x)[1] == 0 for toa in toas])
                @test all([toa.error > time(0.0) for toa in toas])
            end

            @testset "read_tzr_toa" begin
                tzrtoa = model.tzr_toa
                @test tzrtoa.tzr
                @test tzrtoa.level == 0
                @test tzrtoa.error > time(0.0)
                @test modf(tzrtoa.phase.x)[1] == 0
                @test frequency(1e9) < tzrtoa.observing_frequency < frequency(2.5e9)
                @test time(53400.0 * day_to_s) < tzrtoa.value < time(56002.0 * day_to_s)
            end

            @testset "read_param_handler" begin
                param_handler = model.param_handler
                @test length(param_handler.multi_params) ==
                      length(param_handler._default_params_dict)
                @test Set(get_free_param_names(param_handler)) == Set(["F0", "F1", "PHOFF"])
            end

            @testset "read_components" begin
                components = model.components
                @test length(components) == 2
                @test isa(components[1], Spindown)
                @test isa(components[2], PhaseOffset)
            end

            @testset "correct_toas!" begin
                ctoas = correct_toas(model, toas, model.param_handler._default_params_dict)
                @test all(ctoa.value == toa.value for (ctoa, toa) in zip(ctoas, toas))
            end

            @testset "form_residuals" begin
                resids =
                    form_residuals(model, toas, model.param_handler._default_params_dict)
                @test all(isfinite.(resids))
            end

            @testset "calc_chi2" begin
                chi2 = calc_chi2(model, toas, model.param_handler._default_params_dict)
                nfree = length(get_free_param_names(model.param_handler))
                @test chi2 / (length(toas) - nfree) < dimensionless(Float128(1.5))
            end

            @testset "calc_lnlike" begin
                lnlike = calc_lnlike(model, toas, model.param_handler._default_params_dict)
                @test isfinite(lnlike)
            end
        end
    end

    @testset "formatting" begin
        @test format(Vela)
    end
end
