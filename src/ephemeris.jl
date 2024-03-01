export EphemerisVectors

struct EphemerisVectors
    ssb_obs_pos::Vector{GQ{Float64}}
    ssb_obs_vel::Vector{GQ{Float64}}
    obs_sun_pos::Vector{GQ{Float64}}
    obs_jupiter_pos::Vector{GQ{Float64}}
    obs_saturn_pos::Vector{GQ{Float64}}
    obs_venus_pos::Vector{GQ{Float64}}
    obs_uranus_pos::Vector{GQ{Float64}}
    obs_neptune_pos::Vector{GQ{Float64}}
    obs_earth_pos::Vector{GQ{Float64}}

    function EphemerisVectors(
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
        @assert all([xi.d == 1 for xi in ssb_obs_pos]) "Dimension mismatch in ssb_obs_pos."
        @assert all([vi.d == 0 for vi in ssb_obs_vel]) "Dimension mismatch in ssb_obs_vel."
        @assert all([xi.d == 1 for xi in obs_sun_pos]) "Dimension mismatch in obs_sun_pos."
        @assert all([xi.d == 1 for xi in obs_jupiter_pos]) "Dimension mismatch in obs_jupiter_pos."
        @assert all([xi.d == 1 for xi in obs_saturn_pos]) "Dimension mismatch in obs_saturn_pos."
        @assert all([xi.d == 1 for xi in obs_venus_pos]) "Dimension mismatch in obs_venus_pos."
        @assert all([xi.d == 1 for xi in obs_uranus_pos]) "Dimension mismatch in obs_uranus_pos."
        @assert all([xi.d == 1 for xi in obs_neptune_pos]) "Dimension mismatch in obs_neptune_pos."
        @assert all([xi.d == 1 for xi in obs_earth_pos]) "Dimension mismatch in obs_earth_pos."
        @assert dot(ssb_obs_vel, ssb_obs_vel) < 1 "Magnitude of ssb_obs_vel should be less than 1."

        return new(
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
    end
end
