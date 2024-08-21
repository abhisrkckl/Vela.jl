orbital_phase(Δt::GQ, params::NamedTuple, use_fbx::Bool) =
    2π * (
        use_fbx ? taylor_horner_integral(Δt, params.FB, dimensionless(0.0)) :
        (Δt / params.PB - 0.5 * params.PBDOT * (Δt / params.PB)^2)
    )

mean_motion(Δt::GQ, params::NamedTuple, use_fbx::Bool) =
    2π * (use_fbx ? taylor_horner(Δt, params.FB) : 1 / (params.PB + params.PBDOT * Δt))
