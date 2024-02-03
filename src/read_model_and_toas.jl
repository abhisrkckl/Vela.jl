using GeometricUnits
using Quadmath
using HDF5

export read_toas

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

# function _parse_astropy_quantity(astropy_quantity, scale_factor = 1)
#     u = pyimport("astropy" => "units")

#     scaled_quantity = (astropy_quantity * scale_factor).si

#     value = pyconvert(Float64, scaled_quantity.value)

#     if isempty(scaled_quantity.unit.bases) ||
#        pyconvert(Bool, scaled_quantity.unit.bases == pylist([u.rad]))
#         return dimensionless(value)
#     elseif pyconvert(Bool, scaled_quantity.unit.bases == pylist([u.s]))
#         d = pyconvert(Int, scaled_quantity.unit.powers[0])
#         return GQ(value, d)
#     elseif pyconvert(Bool, pyset(scaled_quantity.unit.bases) == pyset([u.s, u.rad]))
#         d = pyconvert(
#             Int,
#             scaled_quantity.unit.powers[scaled_quantity.unit.bases.index(u.s)],
#         )
#         return GQ(value, d)
#     else
#         throw(
#             DomainError(
#                 scaled_quantity,
#                 "The scaled quantity has an unsupported unit. Check the scale_factor.",
#             ),
#         )
#     end
# end

# function _get_scale_factor(pint_param)
#     dmconst = pyimport("pint" => "DMconst")
#     u = pyimport("astropy" => "units")
#     c = pyimport("astropy" => "constants")

#     scale_factors = Dict("DM" => dmconst, "NE_SW" => c.c * dmconst, "PX" => c.c / u.au)

#     param_name = pyconvert(String, pint_param.name)
#     prefix =
#         pyhasattr(pint_param, "prefix") ? pyconvert(String, pint_param.prefix) : Nothing

#     if param_name in keys(scale_factors)
#         return scale_factors[param_name]
#     elseif !isnothing(prefix) && prefix in keys(scale_factors)
#         return scale_factors[prefix]
#     else
#         return 1
#     end
# end

# function _read_params(pint_model)
#     FloatParameter, MaskParameter, PrefixParameter, AngleParameter, MJDParameter = pyimport(
#         "pint.models.parameter" => (
#             "floatParameter",
#             "maskParameter",
#             "prefixParameter",
#             "AngleParameter",
#             "MJDParameter",
#         ),
#     )

#     ignore_params =
#         ["START", "FINISH", "RM", "CHI2", "CHI2R", "TRES", "DMRES", "TZRMJD", "TZRFRQ"]

#     params = Vector{Parameter}([])
#     for param_name in pint_model.params
#         if pyin(param_name, ignore_params) ||
#            !pyisinstance(
#                pint_model[param_name],
#                (
#                    FloatParameter,
#                    MJDParameter,
#                    AngleParameter,
#                    MaskParameter,
#                    PrefixParameter,
#                ),
#            ) ||
#            pyis(pint_model[param_name].quantity, pybuiltins.None)
#             continue
#         end

#         pint_param = pint_model[param_name]

#         frozen = pyconvert(Bool, pint_param.frozen)

#         scale_factor = _get_scale_factor(pint_param)

#         default_quantity =
#             pyisinstance(pint_param, MJDParameter) ?
#             time(f64_convert(pint_param.value * day_to_s)) :
#             _parse_astropy_quantity(pint_param.quantity, scale_factor)

#         param = Parameter(
#             pyconvert(String, param_name),
#             default_quantity,
#             quantity_like(default_quantity, -Inf),
#             quantity_like(default_quantity, Inf),
#             frozen,
#         )

#         push!(params, param)
#     end

#     return params
# end

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
