# [0.0.1] - 2024-07-10
- First tagged version
## Added
- `TimingModel` with components `SolarSystem`, `Dispersion`, `Spindown`, and `PhaseOffset`
- `TOA` type to represent narrowband TOAs
- `SolarSystemEphemeris` type to store solar system ephemerides
- `ParamHandler` class and its friends to convert parameter vectors to named tuples
- Parallel and serial versions of the `chi2` and `lnlike` functions
- `read_model_and_toas` to read data from `HDF5` files (created using `pint2vela.py`)
- `pure_rotator` and `NGC6440E` examples
- GitHub Actions for CI Tests and CodeCov upload
- Basic README file
- MIT Licence
## Changed
## Fixed
## Removed