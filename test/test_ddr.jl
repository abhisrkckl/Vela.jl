@testset "BinaryDDR" begin
    make_ddr_toa(Δt) = TOA(
        time(Double64((54000.0 - epoch_mjd) * day_to_s + Δt)),
        time(1e-6),
        frequency(2.5e9),
        dimensionless(Double64(0.0)),
        default_ephem(),
        1,
    )

    base_params = (
        TASC = time((54000.0 - epoch_mjd) * day_to_s),
        PB = time(86400.0),
        PBDOT = dimensionless(0.0),
        A1 = distance(5.0),
        A1DOT = dimensionless(0.0),
        EPS1 = dimensionless(0.02),
        EPS2 = dimensionless(-0.03),
        M2 = 0.8 * Vela.M_SUN,
        COSI = dimensionless(0.5),
    )

    @testset "regular Kepler helpers" begin
        for λv in (-7π / 3, -π / 4, 0.0, π / 4, 7π / 3)
            λ = dimensionless(λv)
            F, D, c_e, s_e = Vela.solve_F(λ, dimensionless(0.0), dimensionless(0.0))
            λ_red, _ = Vela.reduce_longitude(λ)
            @test value(F) ≈ value(λ_red) atol = 5e-15
            @test D ≈ dimensionless(1.0)
            @test c_e ≈ dimensionless(0.0)
            @test s_e ≈ dimensionless(0.0)
        end

        h = dimensionless(0.18)
        k = dimensionless(0.24)
        F0 = dimensionless(1.1)
        λ = F0 - k * sin(F0) + h * cos(F0)
        F, D, c_e, s_e = Vela.solve_F(λ, h, k)
        @test value(F) ≈ value(F0) atol = 2e-15
        @test value(F - k * sin(F) + h * cos(F)) ≈ value(λ) atol = 2e-15
        @test D ≈ 1 - k * cos(F) - h * sin(F)

        X, Y, _, _, _, _ = Vela.static_XY(F, h, k)
        @test value(X * X + Y * Y) ≈ value((1 - c_e) * (1 - c_e)) atol = 2e-15
        _, _, c_e_tasc, s_e_tasc = Vela.solve_F(dimensionless(0.0), h, k)
        @test value(Vela.q_nu_minus_M(c_e_tasc, s_e_tasc, h * h + k * k)) ≈
              value(Vela.q_at_tasc(h, k)) atol = 2e-15 rtol = 0

        for Δtv in (0.0, 21600.0, 86400.0 * 365.25), pv in (0.0, -1e-12, 1e-12)
            Δt = time(Δtv)
            PB = time(86400.0)
            p = dimensionless(pv)
            λ, λdot = Vela.mean_longitude(Δt, PB, p)
            phase = Δt / PB
            @test value(λ) ≈ value(2π * (phase - 0.5 * p * phase * phase)) atol = 1e-14
            @test value(λdot) ≈ value(2π * (1 / PB - p * Δt / (PB * PB))) atol = 1e-20
        end
    end

    @testset "delay against PINT" begin
        ddr_pk = BinaryDDR(false, false, true, false, false, false)
        expected_pk = (-0.19858979711316324, 4.8951043709516995, -0.2014760001122599)
        for (Δt, expected) in zip((0.0, 21600.0, 43200.0), expected_pk)
            ctoa = correct_toa(ddr_pk, make_ddr_toa(Δt), TOACorrection(), base_params)
            @test value(ctoa.delay) ≈ expected atol = 1e-12 rtol = 0
            @test isfinite(value(ctoa.doppler))
        end

        toa = make_ddr_toa(21600.0)
        correct_toa(ddr_pk, toa, TOACorrection(), base_params)
        @test @ballocated(correct_toa($ddr_pk, $toa, $(TOACorrection()), $base_params)) == 0

        ddr_pheno = BinaryDDR(false, false, false, false, false, false)
        pheno_params = merge(base_params, (GGAMMA = time(0.0), OMDOT = frequency(0.0)))
        expected_pheno = (-0.1984501848750591, 4.8953235689803245, -0.2015923196527405)
        for (Δt, expected) in zip((0.0, 21600.0, 43200.0), expected_pheno)
            ctoa = correct_toa(ddr_pheno, make_ddr_toa(Δt), TOACorrection(), pheno_params)
            @test value(ctoa.delay) ≈ expected atol = 1e-12 rtol = 0
        end

        ddr_fbx = BinaryDDR(true, false, true, false, false, false)
        fbx_params = Base.structdiff(base_params, NamedTuple{(:PB, :PBDOT)})
        fbx_params = merge(fbx_params, (FB = (frequency(1.0 / 86400.0),),))
        expected_fbx = (-0.19858979711316324, 4.8951043709516995, -0.2014760001122591)
        for (Δt, expected) in zip((0.0, 21600.0, 43200.0), expected_fbx)
            ctoa = correct_toa(ddr_fbx, make_ddr_toa(Δt), TOACorrection(), fbx_params)
            @test value(ctoa.delay) ≈ expected atol = 1e-12 rtol = 0
        end

        state = Vela.DDRState(ddr_pk, make_ddr_toa(0.0), TOACorrection(), base_params)
        mp, _ = Vela.pulsar_mass(
            state.n,
            base_params.A1,
            base_params.M2 / Vela.M_SUN,
            base_params.COSI,
        )
        @test value(mp) ≈ 0.7741081166 rtol = 1e-10
        @test value(state.g_gamma) ≈ 0.0071937506 rtol = 1e-8

        Ω = dimensionless(π / 6)
        px_rad = dimensionless(1e-3 * π / (180 * 3600))
        I, J = Vela.IJ_from_v(-px_rad, dimensionless(0.0), Ω)
        expected_geo = (-0.19858982234067224, 4.895104376332743, -0.20147597688237923)
        for (Δt, expected) in zip((0.0, 21600.0, 43200.0), expected_geo)
            st = Vela.DDRState(ddr_pk, make_ddr_toa(Δt), TOACorrection(), base_params)
            B_S = Vela._shapiro_B_S_squared_norm(st.c_e, st.X, st.Y, st.c, st.s, I, J, Ω)
            st_geo = Vela.DDRState(
                st.x,
                st.n,
                st.c,
                st.s,
                st.c_e,
                st.s_e,
                st.g_gamma,
                st.X,
                st.Y,
                st.dX,
                st.dY,
                st.d2X,
                st.d2Y,
                I,
                J,
                B_S,
                st.m2,
                true,
            )
            d = Vela.rømer_einstein_delay(ddr_pk, st_geo)
            dp = Vela.d_rømer_einstein_delay_d_F(ddr_pk, st_geo)
            dpp = Vela.d2_rømer_einstein_delay_d_F2(ddr_pk, st_geo)
            nhat = st_geo.n / (1 - st_geo.c_e)
            d_inv =
                d * (
                    1 - nhat * dp +
                    (nhat * dp) * (nhat * dp) +
                    0.5 * (nhat * nhat) * d * dpp -
                    0.5 * st_geo.s_e / (1 - st_geo.c_e) * (nhat * nhat) * d * dp
                )
            delay = d_inv + Vela.shapiro_delay(ddr_pk, st_geo)
            @test value(delay) ≈ expected atol = 1e-12 rtol = 0
        end

        conjunction_params = merge(
            pheno_params,
            (
                EPS1 = dimensionless(0.0),
                EPS2 = dimensionless(0.0),
                COSI = dimensionless(1e-10),
            ),
        )
        conjunction_state = Vela.DDRState(
            ddr_pheno,
            make_ddr_toa(21600.0),
            TOACorrection(),
            conjunction_params,
        )
        @test conjunction_state.valid
        @test value(conjunction_state.B_S) ≈ 5e-21 rtol = 1e-12
        conjunction_correction = correct_toa(
            ddr_pheno,
            make_ddr_toa(21600.0),
            TOACorrection(),
            conjunction_params,
        )
        @test isfinite(value(conjunction_correction.delay))
    end

    @testset "kinematic period derivative" begin
        mas_to_rad = 1e-3 * π / (180 * 3600)
        px = dimensionless(mas_to_rad) / Vela.AU
        pm = frequency(10 * mas_to_rad / (365.25 * day_to_s))
        n = 2π / base_params.PB
        mc = base_params.M2 / Vela.M_SUN
        mp, _ = Vela.pulsar_mass(n, base_params.A1, mc, base_params.COSI)
        E2 = base_params.EPS1 * base_params.EPS1 + base_params.EPS2 * base_params.EPS2
        p_gw = Vela.pbdot_gw(n, mp, mc, E2)
        p_shk = base_params.PB * pm * pm / px
        p_gal =
            base_params.PB *
            Vela._galactic_acceleration_los(1 / px, dimensionless(1.0), dimensionless(0.2))
        p = p_gw + p_shk + p_gal

        @test value(p_gw) ≈ -1.170160493104432e-14 rtol = 1e-10
        @test value(p_shk) ≈ 2.0988692470710383e-14 rtol = 1e-10
        @test value(p_gal) ≈ -4.722805781109057e-15 rtol = 1e-10
        @test value(p) ≈ 4.564281758557005e-15 rtol = 1e-10

        l = dimensionless(1.0)
        b = dimensionless(0.2)
        sb, cb = sincos(b)
        sl, cl = sincos(l)
        n_gal = (cb * cl, cb * sl, sb)
        rotation = Vela.ICRS_TO_GAL
        n_icrs = (
            rotation[1][1] * n_gal[1] +
            rotation[2][1] * n_gal[2] +
            rotation[3][1] * n_gal[3],
            rotation[1][2] * n_gal[1] +
            rotation[2][2] * n_gal[2] +
            rotation[3][2] * n_gal[3],
            rotation[1][3] * n_gal[1] +
            rotation[2][3] * n_gal[2] +
            rotation[3][3] * n_gal[3],
        )
        raj = atan(n_icrs[2], n_icrs[1])
        decj = atan(n_icrs[3], sqrt(n_icrs[1] * n_icrs[1] + n_icrs[2] * n_icrs[2]))
        compose_params = (
            PB = base_params.PB,
            XPBDOT = dimensionless(0.0),
            PX = px,
            TGEO = base_params.TASC,
            POSEPOCH = base_params.TASC,
            RAJ = raj,
            DECJ = decj,
            PMRA = pm,
            PMDEC = frequency(0.0),
        )
        composed, composed_shk, composed_gal, composed_gw = Vela._compose_p(
            BinaryDDR(false, false, true, true, false, true),
            n,
            mp,
            mc,
            E2,
            compose_params,
        )
        @test value(composed_gw) ≈ value(p_gw) rtol = 1e-13
        @test value(composed_shk) ≈ value(p_shk) rtol = 1e-13
        @test value(composed_gal) ≈ value(p_gal) rtol = 1e-13
        @test value(composed) ≈ value(p) rtol = 1e-13
    end

    @testset "geometry and invalid mass domain" begin
        mas_to_rad = 1e-3 * π / (180 * 3600)
        l_expected = dimensionless(1.0)
        b_expected = dimensionless(0.2)
        sb, cb = sincos(b_expected)
        sl, cl = sincos(l_expected)
        n_gal = (cb * cl, cb * sl, sb)
        rotation = Vela.ICRS_TO_GAL
        n_icrs = (
            rotation[1][1] * n_gal[1] +
            rotation[2][1] * n_gal[2] +
            rotation[3][1] * n_gal[3],
            rotation[1][2] * n_gal[1] +
            rotation[2][2] * n_gal[2] +
            rotation[3][2] * n_gal[3],
            rotation[1][3] * n_gal[1] +
            rotation[2][3] * n_gal[2] +
            rotation[3][3] * n_gal[3],
        )
        raj = atan(n_icrs[2], n_icrs[1])
        decj = atan(n_icrs[3], sqrt(n_icrs[1] * n_icrs[1] + n_icrs[2] * n_icrs[2]))
        geometry_params = merge(
            base_params,
            (
                TGEO = base_params.TASC,
                POSEPOCH = base_params.TASC,
                PX = dimensionless(mas_to_rad) / Vela.AU,
                KOM = dimensionless(π / 6),
                RAJ = raj,
                DECJ = decj,
                PMRA = frequency(0.0),
                PMDEC = frequency(0.0),
            ),
        )
        l, b = Vela._galactic_direction(
            BinaryDDR(false, false, true, false, false, true),
            geometry_params,
        )
        @test value(l) ≈ value(l_expected) atol = 2e-15
        @test value(b) ≈ value(b_expected) atol = 2e-15

        n_ecliptic = Vela._icrs_to_ecliptic(n_icrs)
        elong = atan(n_ecliptic[2], n_ecliptic[1])
        elat = atan(
            n_ecliptic[3],
            sqrt(n_ecliptic[1] * n_ecliptic[1] + n_ecliptic[2] * n_ecliptic[2]),
        )
        ecliptic_params = merge(
            Base.structdiff(geometry_params, NamedTuple{(:RAJ, :DECJ, :PMRA, :PMDEC)}),
            (ELONG = elong, ELAT = elat, PMELONG = frequency(0.0), PMELAT = frequency(0.0)),
        )
        l_ecl, b_ecl = Vela._galactic_direction(
            BinaryDDR(false, true, true, false, false, true),
            ecliptic_params,
        )
        @test value(l_ecl) ≈ value(l_expected) atol = 2e-15
        @test value(b_ecl) ≈ value(b_expected) atol = 2e-15

        ddr_geo = BinaryDDR(false, false, true, false, true, true)
        ctoa = correct_toa(ddr_geo, make_ddr_toa(0.0), TOACorrection(), geometry_params)
        @test isfinite(value(ctoa.delay))
        @test isfinite(value(ctoa.doppler))

        invalid_params = merge(base_params, (M2 = 0.2 * Vela.M_SUN,))
        invalid_ddr = BinaryDDR(false, false, true, false, false, false)
        invalid_ctoa =
            correct_toa(invalid_ddr, make_ddr_toa(0.0), TOACorrection(), invalid_params)
        @test isnan(value(invalid_ctoa.delay))
        @test isnan(value(invalid_ctoa.doppler))
        @test_throws AssertionError correct_toa_delay(
            TOACorrection();
            doppler = dimensionless(Inf),
        )
        @test_throws AssertionError correct_toa_phase(
            TOACorrection();
            delta_spin_frequency = frequency(Inf),
        )

        components = (invalid_ddr, Spindown())
        ph = ParamHandler(
            [Parameter(:COSI, dimensionless(0.5), false, "", 1.0, false)],
            MultiParameter[],
        )
        cosi_prior = SimplePrior{:COSI}(Vela.Uniform(-1.0, 1.0), Vela.USER_DEFINED_PRIOR)
        tzrtoa = make_tzr_toa(
            time(Double64((54000.0 - epoch_mjd) * day_to_s)),
            frequency(2.5e9),
            default_ephem(),
        )
        model = TimingModel(
            "J1234+5678",
            "DE440",
            "TT(BIPM)",
            "TDB",
            time(epoch_mjd * day_to_s),
            components,
            WhiteNoiseKernel(),
            ph,
            tzrtoa,
            (cosi_prior,),
        )
        likelihood_params = merge(
            invalid_params,
            (F = (frequency(100.0), GQ{-2}(0.0)), F_ = frequency(0.0)),
        )
        toas = [make_ddr_toa(0.0)]
        @test Vela.calc_chi2_serial(model, toas, likelihood_params) == Inf
        @test calc_lnlike_serial(model, toas, likelihood_params) == -Inf
        @test calc_lnpost_serial(model, toas, likelihood_params) == -Inf

        woodbury_model = TimingModel(
            "J1234+5678",
            "DE440",
            "TT(BIPM)",
            "TDB",
            time(epoch_mjd * day_to_s),
            components,
            WoodburyKernel(
                WhiteNoiseKernel(),
                (PowerlawRedNoiseGP(1, 0, 2.0),),
                zeros(1, 2),
            ),
            ph,
            tzrtoa,
            (cosi_prior,),
        )
        @test calc_lnlike_serial(woodbury_model, toas, likelihood_params) == -Inf
    end

    display(BinaryDDR(false, false, true, false, false, false))
end
