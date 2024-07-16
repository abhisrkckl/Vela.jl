# Unreleased
## Added
- `MeasurementNoise` component (EFACs and EQUADs)
- CHANGELOG file
- Environment variables for safe Python interoperability in the README file
## Changed
- Split `F0` into two `Float64` variables (`F_` & `F0`) to preserve precision
- Use `DoubleFloats` instead of `Quadmath` to represent TOA values (it's faster)
- Rearrange code and tests into multiple files
- Rearrange test data files
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
- `TOA` type to represent narrowband TOAs
- `CorrectedTOA` type to represent accumulated corrections to a `TOA`.
- `SolarSystemEphemeris` type to store solar system ephemerides
- `ParamHandler` class and its friends to convert parameter vectors to named tuples
- `correct_toa` function
- Parallel and serial versions of the `chi2` and `lnlike` functions
- `read_model_and_toas` to read data from `HDF5` files (created using `pint2vela.py`)
- `pure_rotator` and `NGC6440E` examples
- GitHub Actions for CI Tests and CodeCov upload
- Basic README file
- MIT Licence
## Changed
## Fixed
## Removed