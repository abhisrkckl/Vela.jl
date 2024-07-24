# [Unreleased]
## Added
- `load_pulsar_data()` and `save_pulsar_data()` functions
- `basis_dot` function
## Changed
- Reorganized source files into subdirectories
## Fixed
## Removed
- `read_model_and_toas()` function. Data is now read from `JLSO` files created using `pint2vela`
- `plot_summary()` function. (This is better done in Python.)
- Tests using `PyArray` (This speeds up the test suite)
- Support for Julia 1.9

# [0.0.2] - 2024-07-24
## Added
- `CHANGELOG` file
- Environment variables for safe Python interoperability in the `README` file
- `index` field in `TOA`
- `MeasurementNoise` component (`EFAC`s and `EQUAD`s)
- `get_scale_factors()` function.
- Assertion in `read_params` to make sure that the input has the correct number of values.
- `PhaseJump` component (`JUMP`s)
## Changed
- Split `F0` into two `Float64` variables (`F_` & `F0`) to preserve precision
- Use `DoubleFloats` instead of `Quadmath` to represent TOA values (it's faster)
- Rearrange code and tests into multiple files
- Rearrange test data files
- Use `@spawn` and `fetch` instead of atomic operations for parallel chi2 and likelihood.
- Move chi2 functions into a separate file `chi2.jl`
- Move the higher order functions in `pyinter.jl` to `chi2.jl` and `likelihood.jl`
- Updated `README` to use `pint2vela`
## Fixed
## Removed

# [0.0.1] - 2024-07-10
## Added
- `TimingModel` to represent the timing & noise model
- Hierarchy of `Component` types
- `SolarSystem` component (solar system delays)
- `DispersionTaylor` component (interstellar dispersion as a Taylor series)
- `Spindown` component (pulsar spindown as a Taylor series)
- `PhaseOffset` component (overall phase offset between physical TOAs and the TZR TOA)
- `TimingModel` to represent the timing & noise model
- `TOA` type to represent narrowband TOAs
- `CorrectedTOA` type to represent accumulated corrections to a `TOA`.
- `SolarSystemEphemeris` type to store solar system ephemerides
- `ParamHandler` class and its friends to convert parameter vectors to named tuples
- `correct_toa()` function
- Parallel and serial versions of the `chi2` and `lnlike` functions
- `read_model_and_toas()` to read data from `HDF5` files (created using `pint2vela.py`)
- `pure_rotator` and `NGC6440E` examples
- GitHub Actions for CI Tests and CodeCov upload
- Basic `README` file
- MIT Licence
## Changed
## Fixed
## Removed
