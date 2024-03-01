using GeometricUnits
using Quadmath
using HDF5
using JSON

export read_toas, read_tzr_toa, read_param_handler, read_components, read_model_and_toas

function read_ephem_vectors(toa_data::NamedTuple)::EphemerisVectors
    ssb_obs_pos =
        distance.([toa_data.ssb_obs_pos_x, toa_data.ssb_obs_pos_y, toa_data.ssb_obs_pos_z])
    ssb_obs_vel =
        speed.([toa_data.ssb_obs_vel_x, toa_data.ssb_obs_vel_y, toa_data.ssb_obs_vel_z])
    obs_sun_pos =
        distance.([toa_data.obs_sun_pos_x, toa_data.obs_sun_pos_y, toa_data.obs_sun_pos_z])
    obs_jupiter_pos =
        distance.([
            toa_data.obs_jupiter_pos_x,
            toa_data.obs_jupiter_pos_y,
            toa_data.obs_jupiter_pos_z,
        ])
    obs_saturn_pos =
        distance.([
            toa_data.obs_saturn_pos_x,
            toa_data.obs_saturn_pos_y,
            toa_data.obs_saturn_pos_z,
        ])
    obs_venus_pos =
        distance.([
            toa_data.obs_venus_pos_x,
            toa_data.obs_venus_pos_y,
            toa_data.obs_venus_pos_z,
        ])
    obs_uranus_pos =
        distance.([
            toa_data.obs_uranus_pos_x,
            toa_data.obs_uranus_pos_y,
            toa_data.obs_uranus_pos_z,
        ])
    obs_neptune_pos =
        distance.([
            toa_data.obs_neptune_pos_x,
            toa_data.obs_neptune_pos_y,
            toa_data.obs_neptune_pos_z,
        ])
    obs_earth_pos =
        distance.([
            toa_data.obs_earth_pos_x,
            toa_data.obs_earth_pos_y,
            toa_data.obs_earth_pos_z,
        ])
    return EphemerisVectors(
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

function read_toa(toa_data::NamedTuple, tzr = false)::TOA
    value = time(Float128(toa_data.tdb_int) + Float128(toa_data.tdb_frac))
    error = time(toa_data.error)
    freq = frequency(toa_data.frequency)
    phase = dimensionless(Float128(toa_data.phase_int) + Float128(toa_data.phase_frac))

    # ephem_vectors = read_ephem_vectors(toa_data)

    return TOA(value, error, freq, phase, frequency(-1.0), false, tzr, 0)
end

read_toas(f::HDF5.File)::Vector{TOA} = [read_toa(toa_data) for toa_data in read(f["TOAs"])]

read_tzr_toa(f::HDF5.File) = read_toa(read(f["TZRTOA"])[1], true)

function read_param_handler(f::HDF5.File)
    params_info = JSON.parse(read(f["Parameters"]))

    multi_params = Vector{MultiParameter}()
    for (mpname, param_dicts) in params_info
        params = [
            Parameter(
                pdict["display_name"],
                GQ(pdict["default_value"], pdict["dimension"]),
                pdict["frozen"],
            ) for pdict in param_dicts
        ]
        push!(multi_params, MultiParameter(mpname, params))
    end

    # params = [
    #     Parameter(
    #         pd["name"],
    #         GQ(pd["default_value"], pd["dimension"]),
    #         pd["frozen"],
    #     ) for pd in param_dicts
    # ]
    return ParamHandler(multi_params)
end

function read_components(f::HDF5.File)
    component_dicts = JSON.parse(read(f["Components"]))

    components = Vector{Component}()
    for cdict in component_dicts
        if cdict["name"] == "Spindown"
            push!(components, Spindown())
        elseif cdict["name"] == "PhaseOffset"
            push!(components, PhaseOffset())
        elseif cdict["name"] == "SolarSystem"
            push!(
                components,
                SolarSystem(
                    cdict["ecliptic_coordinates"],
                    cdict["proper_motion"],
                    cdict["planet_shapiro"],
                ),
            )
        elseif cdict["name"] == "Troposphere"
            push!(components, Troposphere())
        elseif cdict["name"] == "DispersionTaylor"
            push!(components, DispersionTaylor())
        elseif cdict["name"] == "SolarWindDispersion"
            swm = cdict["model"]
            @assert swm in [0, 1]
            push!(components, SolarWindDispersion(swm))
        end
    end

    return Tuple(components)
end

function read_model_and_toas(filename::String)
    h5open(filename) do f
        toas = read_toas(f)

        param_handler = read_param_handler(f)
        components = read_components(f)
        tzr_toa = read_tzr_toa(f)
        model = TimingModel(components, param_handler, tzr_toa)

        return model, toas
    end
end
