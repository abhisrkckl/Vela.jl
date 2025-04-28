# Parameters

Each `TimingModel` has a number of model parameters, some of which are free and some frozen.
Some of these parameters are single parameters, like the sky coordinates `RAJ` and `DECJ`, and are 
represented by the `Parameter` type. Some come as families of similar parameters, such as the
pulsar frequency `F0` and its derivatives `F1`, `F2`, etc. These are represented as 
`MultiParameter`s. These contain information about parameter names, their original units (used
by `PINT`), the scaling factors for converting to and from the `Vela.jl` units, default values,
whether they are free/frozen, etc. 
```@docs
Parameter
MultiParameter
```

The `correct_toa` methods expect the parameters to be passed in as a `NamedTuple` containing 
`GQ`s for single `Parameter`s and `NTuple`s of `GQ`s for `MultiParameter`s. In general, the 
samplers the parameter values input from a sampler will be some type of ordered collection of
`Float64`s like a `list`, `ndarray`, `Vector`, etc. This is converted into a `NamedTuple` of 
appropriate structure by the `read_params` function. The opposite can be achieved using the
`read_param_values_to_vector` function. The information needed to do these transformations 
is stored in the `ParamHandler` type.
```@docs
ParamHandler
read_params
read_param_values_to_vector
```

We also provide utility functions for getting ordered lists of parameter names, prefixes, units, etc. 
(in the `PINT` convention).
```@docs
get_free_param_names
get_free_param_prefixes
get_free_param_units
get_free_param_labels
get_scale_factors
```

A list of parameters and their units can be found [here](https://nanograv-pint.readthedocs.io/en/latest/timingmodels.html#supported-parameters).