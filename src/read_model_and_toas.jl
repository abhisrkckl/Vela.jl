export read_toas, read_tzr_toa, read_param_handler, read_components, read_model_and_toas

function read_ephemeris(toa_data::NamedTuple)::SolarSystemEphemeris
    ssb_obs_pos =
        distance.((toa_data.ssb_obs_pos_x, toa_data.ssb_obs_pos_y, toa_data.ssb_obs_pos_z))
    ssb_obs_vel =
        speed.((toa_data.ssb_obs_vel_x, toa_data.ssb_obs_vel_y, toa_data.ssb_obs_vel_z))
    obs_sun_pos =
        distance.((toa_data.obs_sun_pos_x, toa_data.obs_sun_pos_y, toa_data.obs_sun_pos_z))
    obs_jupiter_pos =
        distance.((
            toa_data.obs_jupiter_pos_x,
            toa_data.obs_jupiter_pos_y,
            toa_data.obs_jupiter_pos_z,
        ))
    obs_saturn_pos =
        distance.((
            toa_data.obs_saturn_pos_x,
            toa_data.obs_saturn_pos_y,
            toa_data.obs_saturn_pos_z,
        ))
    obs_venus_pos =
        distance.((
            toa_data.obs_venus_pos_x,
            toa_data.obs_venus_pos_y,
            toa_data.obs_venus_pos_z,
        ))
    obs_uranus_pos =
        distance.((
            toa_data.obs_uranus_pos_x,
            toa_data.obs_uranus_pos_y,
            toa_data.obs_uranus_pos_z,
        ))
    obs_neptune_pos =
        distance.((
            toa_data.obs_neptune_pos_x,
            toa_data.obs_neptune_pos_y,
            toa_data.obs_neptune_pos_z,
        ))
    obs_earth_pos =
        distance.((
            toa_data.obs_earth_pos_x,
            toa_data.obs_earth_pos_y,
            toa_data.obs_earth_pos_z,
        ))
    return SolarSystemEphemeris(
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

function read_toa(toa_data::NamedTuple, index, tzr = false)::TOA
    value = time(Double64(toa_data.tdb_int, toa_data.tdb_frac))
    error = time(toa_data.error)
    freq = frequency(toa_data.frequency)
    pulse_number = -dimensionless(Double64(toa_data.phase_int, toa_data.phase_frac))
    ephem = read_ephemeris(toa_data)

    return TOA(value, error, freq, pulse_number, false, tzr, ephem, index)
end

read_toas(f::HDF5.File)::Vector{TOA} =
    [read_toa(toa_data, index) for (index, toa_data) in enumerate(read(f["TOAs"]))]

read_tzr_toa(f::HDF5.File) = read_toa(read(f["TZRTOA"])[1], 0, true)

_read_single_parameter(spdict::Dict) = Parameter(
    Symbol(spdict["name"]),
    GQ(spdict["default_value"], spdict["dimension"]),
    spdict["frozen"],
    spdict["original_units"],
    spdict["unit_conversion_factor"],
)

_read_multi_parameter(mpdict::Dict) =
    MultiParameter(Symbol(mpdict["name"]), map(_read_single_parameter, mpdict["elements"]))

function read_param_handler(f::HDF5.File)
    single_params_info, multi_params_info = JSON.parse(read(f["Parameters"]))

    single_params = map(_read_single_parameter, single_params_info)
    multi_params = map(_read_multi_parameter, multi_params_info)

    return ParamHandler(single_params, multi_params)
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
                SolarSystem(cdict["ecliptic_coordinates"], cdict["planet_shapiro"]),
            )
            # elseif cdict["name"] == "Troposphere"
            # #     push!(components, Troposphere())
        elseif cdict["name"] == "DispersionTaylor"
            push!(components, DispersionTaylor())
            # elseif cdict["name"] == "SolarWindDispersion"
            #     swm = cdict["model"]
            #     @assert swm in [0, 1]
            #     push!(components, SolarWindDispersion(swm))
        end
    end

    return Tuple(components)
end

read_info(f::HDF5.File) = JSON.parse(read(f["Info"]))

function read_model_and_toas(filename::String)
    HDF5.h5open(filename) do f
        toas = read_toas(f)

        param_handler = read_param_handler(f)
        components = read_components(f)
        tzr_toa = read_tzr_toa(f)
        info = read_info(f)

        model = TimingModel(
            info["pulsar_name"],
            info["ephem"],
            info["clock"],
            info["units"],
            components,
            param_handler,
            tzr_toa,
        )

        return model, toas
    end
end
