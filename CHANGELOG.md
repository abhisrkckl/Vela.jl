# [Unreleased]
## Added
- `KinematicDelayComponent` as the abstract base class for `DelayComponent`s that also contribute a Doppler factor.
- Doppler factor in `BinaryELL1`
- `BinaryDD` and `BinaryDDH` models
- Doppler factor in `BinaryDD`
- `get_free_param_labels()` function
- Examples - `J0613-0200.sim`, `J1856-3754.sim`, `J1802-2124.sim`, `J0955-6150.sim`, `J1208-5936.sim`, `sim6`, `sim_dd`
- Test Python formatting using `black`
## Changed
- Exposed `cheat_prior_scale` and `custom_prior_dists` options in `read_model_and_toas()`
- Made changes according to the `GeometricUnits` API changes (`GQ` now represents dimensions as a type parameter)
- Use the github version of `PINT` for testing
## Fixed
- Tests now handle par files without `PHOFF` properly.
- Proper motion computation
- Default of `CorrectedTOA.spin_frequency`
- Shapiro delay expression for `BinaryDDBase`
- True anomaly computation in `DDState`
## Removed

# [0.0.3] - 2024-08-22
## Added
- `load_pulsar_data()` and `save_pulsar_data()` functions
- `par_tim_to_jlso()` function and `par_tim-to-jlso` script in `pint2vela`
- `Prior` as the abstract base class for prior distributions
- `SimplePriorBase`, `SimplePrior` and `SimplePriorMulti` to represent priors that can be factorized parameter-wise.
- `distr()`, `lnprior()` and `prior_transform()` functions.
- `get_lnprior_func()` and `get_prior_transform_func()` functions.
- `get_default_priors()` function in `pint2vela`
- More methods for `get_free_param_names()`, `read_param_values_to_vector()`, and `get_scale_factors()` for convenience
- Added CI tests for `pint2vela`
- An alternative implementation of `PhaseJump` for mutually exclusive JUMPs
- `basis_dot` function
- Simple solar wind model (Edwards+ 2006) (`SolarWind`)
- Variable-index chromatic delay as a Taylor expansion (`ChromaticTaylor`)
- Auto-generation of HTML documentation using `Documenter`
- `docs-CI` tests
- Examples - `sim1`, `sim_jump`, `sim_jump_ex`, `sim_fdjumpdm`, `sim_sw`, `sim_cm`, `sim_fd`, `J0613-0200.InPTA.NB`, `J1857+0943.InPTA.NB`
- Tests corresponding to the example datasets
- System-dependent DM offsets (`DispersionOffset`)
- `compare_residuals.py` script in `examples`.
- Codecov upload for `pint2vela`
- Fourier series representation of achromatic red noise (`WaveX`), DM noise (`DMWaveX`), and chromatic noise (`ChromaticCM`)
- `get_lnpost_func` function
- Frequency-dependent profile variability corrections (`FrequencyDependent`)
- Memory allocation tests for all components in `test_components.jl`
- `mean_anomaly` and `mean_motion` functions
- `BinaryELL1` model
- Memory allocation tests in the `pint2vela` test suite.
## Changed
- Reorganized source files into subdirectories
- Replaced `par` and `tim` files for testing with `JLSO` files
- Moved `pint2vela.py` to separate repo, added it as a submodule.
- Rearranged `pint2vela` code into multiple files.
- Merged `pint2vela` into the main repo.
- Replaced `pint2vela` test datasets with symlinks.
- Updated `README.md`
- Moved `setup.py` from `Vela.jl/pint2vela` to `Vela.jl`
## Fixed
- `show` method for `MeasurementNoise`
- Copy the `toas` inside `get_lnlike_parallel_func`, `get_lnlike_serial_func`, `get_chi2_serial_func`, and `get_chi2_parallel_func` to avoid repeated allocations. 
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
