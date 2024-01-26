from pint.models import get_model_and_toas, TimingModel
from pint.models.parameter import floatParameter, MJDParameter, AngleParameter, maskParameter, prefixParameter
from pint.toa import TOAs
from pint import DMconst
import numpy as np
from astropy import units as u
from astropy import constants as c
from astropy.table import Table
import json

day_to_s = 86400


def write_tim(toas: TOAs, file_prefix: str):
    assert toas.planets
    assert toas.get_pulse_numbers() is not None

    tdbs = toas.table["tdbld"].value * day_to_s
    tdbs_frac, tdbs_int = np.modf(tdbs)
    tdbs_frac = tdbs_frac.astype(float)
    tdbs_int = tdbs_int.astype(float)

    phases = -toas.table["pulse_number"].value + toas.table["delta_pulse_number"].value
    phases_frac, phases_int = np.modf(phases)
    phases_frac = phases_frac.astype(float)
    phases_int = phases_int.astype(float)

    errs = toas.get_errors().si.value
    freqs = toas.get_freqs().si.value

    ssb_obs_poss = toas.table["ssb_obs_pos"].quantity.to_value("lightsecond")
    ssb_obs_vels = toas.table["ssb_obs_vel"].quantity.to_value("lightsecond/s")
    obs_sun_poss = toas.table["obs_sun_pos"].quantity.to_value("lightsecond")

    obs_jupiter_poss = toas.table["obs_jupiter_pos"].quantity.to_value("lightsecond")
    obs_saturn_poss = toas.table["obs_saturn_pos"].quantity.to_value("lightsecond")
    obs_venus_poss = toas.table["obs_venus_pos"].quantity.to_value("lightsecond")
    obs_uranus_poss = toas.table["obs_uranus_pos"].quantity.to_value("lightsecond")
    obs_neptune_poss = toas.table["obs_neptune_pos"].quantity.to_value("lightsecond")
    obs_earth_poss = toas.table["obs_earth_pos"].quantity.to_value("lightsecond")

    data = np.vstack(
        (
            tdbs_int,
            tdbs_frac,
            phases_int,
            phases_frac,
            errs,
            freqs,
            ssb_obs_poss.T,
            ssb_obs_vels.T,
            obs_sun_poss.T,
            obs_jupiter_poss.T,
            obs_saturn_poss.T,
            obs_venus_poss.T,
            obs_uranus_poss.T,
            obs_neptune_poss.T,
            obs_earth_poss.T,
        )
    ).T

    colnames = [
        "tdb_int",
        "tdb_frac",
        "phase_int",
        "phase_frac",
        "error",
        "frequency",
        "ssb_obs_pos_x",
        "ssb_obs_pos_y",
        "ssb_obs_p os_z",
        "ssb_obs_vel_x",
        "ssb_obs_vel_y",
        "ssb_obs_vel_z",
        "obs_sun_pos_x",
        "obs_sun_pos_y",
        "obs_sun_pos_z",
        "obs_jupiter_pos_x",
        "obs_jupiter_pos_y",
        "obs_jupiter_pos_z",
        "obs_saturn_pos_x",
        "obs_saturn_pos_y",
        "obs_saturn_pos_z",
        "obs_venus_pos_x",
        "obs_venus_pos_y",
        "obs_venus_pos_z",
        "obs_uranus_pos_x",
        "obs_uranus_pos_y",
        "obs_uranus_pos_z",
        "obs_neptune_pos_x",
        "obs_neptune_pos_y",
        "obs_neptune_pos_z",
        "obs_earth_pos_x",
        "obs_earth_pos_y",
        "obs_earth_pos_z",
    ]

    table = Table(data=data, names=colnames)
    table.write(f"{file_prefix}.tim.hdf5", format="hdf5", path="TOAs")


def _get_scale_factor(param):
    scale_factors = {
        "DM" : DMconst, 
        "NE_SW" : c.c * DMconst, 
        "PX" : c.c / u.au,
    }

    if param.name in scale_factors:
        return scale_factors[param.name]
    elif hasattr(param, "prefix") and param.prefix in scale_factors:
        return scale_factors[param.prefix]
    else:
        return 1

def _parse_quantity(quantity, scale_factor=1):
    scaled_quantity = (quantity * scale_factor).si

    value = scaled_quantity.value

    if len(scaled_quantity.unit.bases) == 0 or scaled_quantity.unit.bases == [u.rad]:
        return value, 0
    elif scaled_quantity.unit.bases == [u.s]:
        d = scaled_quantity.unit.powers[0]
        return value, d
    elif set(scaled_quantity.unit.bases) == {u.s, u.rad}:
        d = scaled_quantity.unit.powers[scaled_quantity.unit.bases.index(u.s)]
        return value, d
    else:
        raise ValueError(
                "The scaled quantity has an unsupported unit. Check the scale_factor.",
            )

def write_params(model:TimingModel, file_prefix:str):
    ignore_params = ["START", "FINISH", "RM", "CHI2", "CHI2R", "TRES", "DMRES", "TZRMJD", "TZRFRQ", "SWM"]

    params = []
    for param_name in model.params:
        if param_name in ignore_params or not isinstance(
               model[param_name],
               (
                   floatParameter,
                   MJDParameter,
                   AngleParameter,
                   maskParameter,
                   prefixParameter,
               ),
           ) or model[param_name].quantity is None:
            continue

        param = model[param_name]

        frozen = param.frozen

        scale_factor = _get_scale_factor(param)

        value, dim = (param.value * day_to_s, 1) if isinstance(param, MJDParameter) else _parse_quantity(param.quantity, scale_factor)

        params.append({
            "name": param_name,
            "default_value": float(value),
            "dimension": dim,
            "frozen": frozen
        })
    
    with open(f"{file_prefix}.params.json", "w") as f:
        json.dump(params, f, indent=True)
    # print(param_name, value, int(frozen))
        
def write_model(model:TimingModel, file_prefix:str):
    component_names = list(model.components.keys())

    components = []
    for component_name in component_names:
        if component_name in ["AbsPhase", "SolarSystemShapiro"]:
            continue
        elif component_name == "Spindown":
            components.append({
                "name": "Spindown",
                "number_of_terms": len(model.get_spin_terms())
            })
        elif component_name.startswith("Astrometry"):
            components.append({
                "name": "SolarSystem",
                "coordinates": component_name[len("Astrometry"):].upper(),
                "planet_shapiro": model.PLANET_SHAPIRO.value
            })
        elif component_name == "SolarWindDispersion":
            components.append({
                "name": "SolarWind",
                "model": int(model.SWM.value)
            })
        elif component_name == "DispersionDM":
            components.append({
                "name": "DispersionTaylor",
                "number_of_terms": len(model.get_DM_terms())
            })
        else:
            components.append({
                "name": component_name
            })
    
    with open(f"{file_prefix}.model.json", "w") as f:
        json.dump(components, f, indent=True)