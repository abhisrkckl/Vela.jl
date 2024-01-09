using GeometricUnits
using Quadmath
using PythonCall

export _read_toas, _read_tzr_toa

parse_quantity(quantity_func, convert_func, values) =
    [quantity_func(convert_func(val)) for val in values]

f128_convert(val) = parse(Float128, string(val))
f64_convert(val) = pyconvert(Float64, val)
vec3_convert(vec) = pyconvert(Vector{Float64}, vec)

const day_to_s = 86400

function _read_toas(pint_toas)::Vector{TOA}
    vals = parse_quantity(time, f128_convert, pint_toas.table["tdbld"].value * day_to_s)
    errs = parse_quantity(time, f64_convert, pint_toas.table["error"].quantity.si.value)
    freqs =
        parse_quantity(frequency, f64_convert, pint_toas.table["freq"].quantity.si.value)
    ssb_obs_poss = parse_quantity(
        (xs) -> distance.(xs),
        vec3_convert,
        pint_toas.table["ssb_obs_pos"].quantity.to_value("lightsecond"),
    )
    ssb_obs_vels = parse_quantity(
        (xs) -> speed.(xs),
        vec3_convert,
        pint_toas.table["ssb_obs_vel"].quantity.to_value("lightsecond/second"),
    )
    obs_sun_poss = parse_quantity(
        (xs) -> distance.(xs),
        vec3_convert,
        pint_toas.table["obs_sun_pos"].quantity.to_value("lightsecond"),
    )
    phases = parse_quantity(
        dimensionless,
        f128_convert,
        (
            pint_toas.table["delta_pulse_number"].value -
            pint_toas.table["pulse_number"].value
        ),
    )

    return [
        TOA(
            val,
            err,
            freq,
            phase,
            (
                iszero(ssb_obs_pos) ? nothing :
                EphemerisVectors(ssb_obs_pos, ssb_obs_vel, obs_sun_pos)
            ),
        ) for (val, err, freq, phase, ssb_obs_pos, ssb_obs_vel, obs_sun_pos) in
        zip(vals, errs, freqs, phases, ssb_obs_poss, ssb_obs_vels, obs_sun_poss)
    ]
end

function _read_tzr_toa(pint_model, pint_toas)::TOA
    pint_tzr_toa = pint_model.get_TZR_toa(pint_toas)

    tzr_tdb = time(f128_convert(pint_tzr_toa.table["tdbld"].value[0] * day_to_s))
    tzr_frq = frequency(f64_convert(pint_tzr_toa.table["freq"].quantity.si.value[0]))
    tzr_ssb_obs_pos =
        distance.(
            vec3_convert(
                pint_tzr_toa.table["ssb_obs_pos"].quantity.to_value("lightsecond")[0],
            )
        )
    tzr_ssb_obs_vel =
        speed.(
            vec3_convert(
                pint_tzr_toa.table["ssb_obs_vel"].quantity.to_value("lightsecond/second")[0],
            )
        )
    tzr_obs_sun_pos =
        distance.(
            vec3_convert(
                pint_tzr_toa.table["obs_sun_pos"].quantity.to_value("lightsecond")[0],
            )
        )
    tzr_ephem = EphemerisVectors(tzr_ssb_obs_pos, tzr_ssb_obs_vel, tzr_obs_sun_pos)

    return make_tzr_toa(tzr_tdb, tzr_frq, tzr_ephem)
end

# function _read_model(pint_model, pint_toas)
#     tzr_toa = _read_tzr_toa(pint_model, pint_toas)
#     params = _read_params(pint_model)
    
#     return TimingModel(ParamHandler(params), tzr_toa)
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
#     for param_name in param_names
#         if pyin(param_name, ignore_params)
#             or ! pyisinstance(
#                 pint_model[param_name],
#                 (
#                     FloatParameter,
#                     MJDParameter,
#                     AngleParameter,
#                     MaskParameter,
#                     PrefixParameter,
#                 ),
#             )
#             continue
#         end

#         frozen = pyconvert(Bool, pint_model[param_name].frozen)

#         default_quantity =
#             pyisinstance(pint_model[param_name], MJDParameter) ?
#             time(f64_convert(pint_model[param_name].value * day_to_s)) :
#             parse_astropy_quantity(pint_model[param_name].quantity)

#         param = Parameter(param_name, default_quantity, angle(-Inf), angle(Inf), frozen)

#         push!(params, param)
#     end

#     return params
# end

# function _parse_astropy_quantity(astropy_quantity, scale_factor=1)
#     u = pyimport("astropy" => "units")

#     scaled_quantity = (astropy_quantity * scale_factor).si

#     if isempty(scaled_quantity.unit.bases)
#         return dimensionless(scaled_quantity.value)
#     elseif scaled_quantity.unit.bases == pylist([u.rad])
#         return dimensionless(scaled_quantity.value)
#     elseif scaled_quantity.unit.bases == pylist([u.s])
#         return GQ(scaled_quantity.value, scaled_quantity.unit.powers[0])
#     elseif pyset(scaled_quantity.unit.bases) == pyset([u.s, u.rad])
#         return GQ(scaled_quantity.value, scaled_quantity.unit.powers[scaled_quantity.unit.bases.index(u.s)])
#     else
#         throw(DomainError(scaled_quantity, "The scaled quantity has an unsupported unit. Check the scale_factor."))
#     end
# end

# function read_model_and_toas(parfile::String, timfile::String)
#     get_model_and_toas = pyimport("pint.models" => "get_model_and_toas")
#     setup_log = pyimport("pint.logging" => "setup")

#     setup_log(level = "WARNING")

#     pint_model, pint_toas = get_model_and_toas(parfile, timfile, add_tzr_to_model = true)
#     pint_toas.compute_pulse_numbers(pint_model)

#     toas = _read_toas(pint_toas)
#     model = _read_model(pint_model, pint_toas)

#     return model, toas
# end
