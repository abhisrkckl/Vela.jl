using GeometricUnits
using Quadmath
using HDF5
using JSON

export read_toas, read_tzr_toa, read_param_handler, read_components

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

    ephem_vectors = read_ephem_vectors(toa_data)

    return TOA(value, error, freq, phase, ephem_vectors, 0, tzr)
end

read_toas(f::HDF5.File)::Vector{TOA} = [read_toa(toa_data) for toa_data in read(f["TOAs"])]

read_tzr_toa(f::HDF5.File) = read_toa(read(f["TZRTOA"])[1], true)

function read_param_handler(f::HDF5.File)
    param_dicts = JSON.parse(read(f["Parameters"]))
    params = [
        Parameter(
            pd["name"],
            GQ(pd["default_value"], pd["dimension"]),
            GQ(-Inf, pd["dimension"]),
            GQ(Inf, pd["dimension"]),
            pd["frozen"],
        ) for pd in param_dicts
    ]
    return ParamHandler(params)
end

function read_components(f::HDF5.File)
    component_dicts = JSON.parse(read(f["Components"]))

    components = Vector{Component}()
    for cdict in component_dicts
        if cdict["name"] == "Spindown"
            push!(components, Spindown(cdict["number_of_terms"]))
        elseif cdict["name"] == "PhaseOffset"
            push!(components, PhaseOffset())
        elseif cdict["name"] == "SolarSystem"
            push!(components, SolarSystem(cdict["ecliptic_coordinates"], cdict["planet_shapiro"]))
        end
    end

    return components
end

# function _read_model(pint_model, pint_toas)
#     tzr_toa = _read_tzr_toa(pint_model, pint_toas)
#     params = _read_params(pint_model)

#     return TimingModel(ParamHandler(params), tzr_toa)
# end

# function read_model_and_toas(parfile::String, timfile::String)
#     setup_log = pyimport("pint.logging" => "setup")
#     setup_log(level = "WARNING")

#     get_model_and_toas = pyimport("pint.models" => "get_model_and_toas")
#     pint_model, pint_toas = get_model_and_toas(parfile, timfile, add_tzr_to_model = true)
#     pint_toas.compute_pulse_numbers(pint_model)

#     if !pyin("PhaseOffset", pint_model.components)
#         PhaseOffset = pyimport("pint.models" => "PhaseOffset")
#         pint_model.add_component(PhaseOffset())
#         pint_model.PHOFF.frozen = pybuiltins.False
#     end

#     toas = _read_toas(pint_toas)
#     model = _read_model(pint_model, pint_toas)

#     return model, toas
# end
