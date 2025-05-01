# Unreleased
## Added
- Support for wideband TOAs in `WoodburyKernel`
- Symlink for `pyvela/examples` in the root directory
- Instructions for updating `Vela.jl` and `pyvela`
- arXiv link in README
- CITATION file
- In `SPNTA.from_pint()`, compute planetary ephemerides if they are absent in the input `TOAs` object.
- Support for log-spaced frequencies in Fourier GP components
- Display default values in `$ pyvela --help`
- Check that `Vela` and `pyvela` have the same version.
- `PriorSourceType` enum, `source_type` attribute in `Prior` types
- `SPNTA.full_prior_dict()` and `SPNTA.save_results()` methods
- `unscale_prior_args()` function
- Save the full prior information along with the results
- Plot priors in `pyvela-plot`
## Changed
- Residual plot in `pyvela-plot` script
- Made GP noise marginalization the default in `SPNTA`
- Move `info_dict()`, `save_new_parfile()`, `save_resids()` from `pyvela_script` to `SPNTA`
## Fixed
- Julia example script `run_example.jl`
- In `_gls_lnlike_serial`, return  -inf if `Σinv` is not positive definite.
## Removed
- Saving maximum-posterior par file in `pyvela` script

# 0.0.9
## Added
- Store `toas_pint` in `SPNTA` (This is `None` when a `JLSO` file is loaded.)
- `get_num_timing_params()` function and `SPNTA.ntmdim` attribute
- Marginalized GP noise using Woodbury lemma
  - `WoodburyKernel` type with inner `WhiteNoiseKernel`
  - `is_gp_noise()`, `get_gp_npars()`, and `calc_noise_weights_inv()` functions
  - corresponding likelihood implementation
  - `marginalize_gp_noise` option in `SPNTA` class and `pyvela` script
  - Plot whitened residuals in `pyvela-script`
  - `has_marginalized_gp_noise` property and `get_marginalized_gp_noise_realization()` and `whitened_time_residuals()` methods in `SPNTA`
- `VelaFitter` class
## Changed
- Made some properties of `SPNTA` cached.
- Store a deepcopy of `model_pint` instead of the mutated original in `SPNTA`. The mutated version is stored as `model_pint_modified`.
- The parameter order now has timing parameters first and then noise parameters.
- Refactored parameter attribute functions.
- `get_kernel()` now supports Woorbury kernel with ECORRs.
- Big documentation update
- Made `TOA` creation slightly faster
## Fixed
- Priors in examples
- Creation of ECORR mask in `get_kernel` in `pyvela` (This was mistakenly being created from EFACs)
## Removed
- `pure_rotator.jlso` test data file
- Stub for `Troposphere` component

# 0.0.8
## Added
- Use `HDFBackend` in `pyvela` script
- `scaled_dm_uncertainties()` and `dm_residuals()` methods in `SPNTA` class
- `pyvela-compare` script
- Save maximum-posterior and median par files in `pyvela` script
- `get_free_param_prefixes` function and `SPNTA.param_prefixes` property
- Save free parameter prefixes, default values, and post-fit residuals in the `pyvela` script.
- `-T` option in `pyvela` script to save the true parameters for simulation studies.
- More help messages in `pyvela` script
- `pyvela-plot` script
- `check` option in `SPNTA` constructor
- Force rewrite (`-f`) option in `pyvela` script
- Take JLSO files as input in the `pyvela` script (`-J` option)
- Informative error messages in `assert` statements
## Changed
- Renamed `SPNTA.maxlike_params` -> `SPNTA.default_params`
- Throw an error if the output directory exists in `pyvela` script
- Save only the basename of input files in the summary file in in `pyvela` script
- Changed `par_tim-to-jlso.py` into an installable script `pyvela-jlso`.
## Fixed
- Correctly avoid likelihood computation when the parameter is outside prior range.
- Added `DoubleFloats`, `Distributions`, and `GeometricUnits` to installation instructions
- Improved error message for missing prior info.
## Removed
- `compare_residuals.py` example script
- Separate `README` for `pyvela`

# 0.0.7
## Added
- `time_residuals`, `scaled_toa_uncertainties`, and `model_dm`, `from_pint` methods in `SPNTA`
- NE_SW derivatives in `SolarWind`
- `Pulsar` class
- `pyvela` script
- `get_free_param_units` function
- "Known issues" page in documentation
- DMX: Piecewise-constant model for DM (`DispersionPiecewise`) 
- CMX: Piecewise-constant model for CM (`ChromaticPiecewise`) 
- `BinaryELL1k` model
## Changed
- Reduced the number of `pyvela` tests
## Fixed
- DMX and CMX cannot be used alongside other stochastic DM variation models 
## Removed
- `test_alloc` and `test_likelihood` from `test_data_files.py`
- Unnecessary test & example data files

# 0.0.6
## Added
- `pyvela.SPNTA` class
- `get_unit_conversion_factor` function
- `prior_scaling` methods and `scale_prior_args` function
## Changed
- Reorganized `pyvela` code
- Rerun failures in `pyvela` CI tests.
- Example scripts now use `pyvela.SPNTA`
- Updated documentation to use `pyvela.SPNTA` and installation instructions
- Made the repository public and deployed the documentation website.
- Migrated Python package settings from `setup.py` to `pyproject.toml`
- Migrated from `Codecov` to `Coveralls` for code coverage
- Prior JSON files now accept parameters in "normal" units.
- `SPNTA` constructor now accepts "raw" prior dictionaries similar to the prior files rather than those containing actual `Distribution` objects.
- Renamed `parse_custom_prior_file` to `process_custom_priors`
- Renamed `read_model_and_toas` to `convert_model_and_toas`
- Updated lowest versions for Python dependencies
## Fixed
- Memory allocation in GP components
## Removed
- Some unnecessary test datasets
- Installation instructions from README file

# [0.0.5]
## Added
- Check for `ECORR` exclusivity
- Vectorize option in `get_lnpost_func`
- GP noise models - `PowerlawRedNoiseGP`, `PowerlawDispersionNoiseGP`, `PowerlawChromaticNoiseGP`
- Examples - `sim3.gp`, `sim4.gp`, `sim6.gp`
- Detailed documentation
- Abstract base classes `RedNoiseBase`, `DispersionNoiseBase`, `ChromaticNoiseBase`, `FrequencyDependentBase`
## Changed
- Renamed `pint2vela` -> `pyvela`
- Infer `is_tzr(toa)` from `toa.index`
- Infer `is_barycentered(toa)` from `toa.ephem.ssb_obs_pos`
- Renamed `CorrectedTOA`, `CorrectedWidebandTOA`, and `CorrectedDMInfo` respectively to `TOACorrection`, `WidebandTOACorrection`, and `DMInfoCorrection`
- Removed the `TOA` from `TOACorrection` (it's passed separately to functions). This reduces copy overhead. 
- Avoid unnecessary repeated computations in `Spindown` 
- Specialized methods of `taylor_horner` and `taylor_horner_integral` for faster execution
- Split `correct_toa` into more specialized methods `correct_toa_delay`, `correct_toa_phase`, and `correct_toa_error`
- Multiplication instead of power in dispersion delay
- Subtract `PEPOCH` from all TOA and `MJDParameter` values; add `epoch` attribute to `TimingModel`
- Moved `basis_dot` to `jump.jl`
- Better priors for `KOM`, `SINI`, `KIN`, `STIGMA`, `SHAPMAX`, `EFAC`, etc.
## Fixed
- Unnecessary repetition of `sin` and `cos` in ecliptic coordinate conversion
- Bug in `_ecorr_lnlike_group`
## Removed
- `level` attribute from `TOACorrection`
- `obs_earth_pos` from `SolarSystemEphemeris`

# [0.0.4]
## Added
- Doppler factor in `BinaryELL1`
- `BinaryDD`, `BinaryDDH`, `BinaryDDS`, `BinaryELL1H`, and `BinaryDDK` models
- Doppler factor in `BinaryDD`
- `get_free_param_labels()` function
- Examples - `J0613-0200.sim`, `J1856-3754.sim`, `J1802-2124.sim`, `J0955-6150.sim`, `J1208-5936.sim`, `J2302+4442.sim`, `J1227-6208.sim`, `sim6`, `sim_dd`, `sim_ddk`, `sim_sw.wb`, `sim_dmjump`, `sim_dmwn`, `sim2`, `sim_glitch`
- Test Python formatting using `black`
- Use `BinaryDD` for par files with the BT model.
- Basic wideband timing implementation
  - `TOABase` as the base type for `TOA` and `WidebandTOA`
  - `WidebandTOA` as the composition of `TOA` and `DMInfo`. The latter contains DM measurement and error.
  - `CorrectedTOABase` as the base type for `CorrectedTOA` and `CorrectedWidebandTOA`
  - `CorrectedWidebandTOA` as the composition of `CorrectedTOA` and `CorrectedDMInfo`. The latter contains model DM and DM error.
  - Methods of `correct_toa()`, `form_residual()`, `calc_chi2()`, `calc_lnlike()`, and `calc_lnpost()` that act on `WidebandTOA`s and relatives.
  - `load_pulsar_data()` and `save_pulsar_data()` now work with wideband TOAs
  - `pint2vela` can now read wideband TOAs.
- `degrees_of_freedom()` and `reduced_chi2()` functions
- Wideband DM offsets (`DMJUMP`s)
- `DispersionMeasurementNoise` component (`DMEFAC`s and `DMEQUAD`s) for wideband data
- `ECORR` implementation
  - `Kernel` as the abstract base class for all likelihood kernels; added the `kernel` member in `TimingModel`.
  - `WhiteNoiseKernel` marks white noise-only likelihood computation.
  - `EcorrKernel` marks likelihood computation with only white noise and ECORR.
  - Methods of `calc_chi2` and `calc_lnlike` specialized for `WhiteNoiseKernel` and `EcorrKernel`
  - `ecorr_sort` function in `pint2vela`
  - `get_kernel` function in `pint2vela` that constructs `Kernel`s based on the `PINT` `TimingModel` object
- `Glitch` component
- `FrequencyDependentJump` component (`FDJUMP`s)
- Documentation - Getting Started and Explanation pages
## Changed
- Exposed `cheat_prior_scale` and `custom_prior_dists` options in `read_model_and_toas()`
- Made changes according to the `GeometricUnits` API changes (`GQ` now represents dimensions as a type parameter)
- Use the github version of `PINT` for testing
- Split `test_components.jl` into multiple files.
- Simplify the type hierarchy of `Component`s. Now all `Component`s are `TOA`-uncorrelated by definition.
## Fixed
- Tests now handle par files without `PHOFF` properly.
- Proper motion computation
- Default of `CorrectedTOA.spin_frequency`
- Shapiro delay expression for `BinaryDDBase`
- True anomaly computation in `DDState`
- Scale factor computation in `get_default_prior()`
- Bug in `correct_toa()` (`ssb_psr_pos` was being set incorrectly)
- `black` failure now shows up as a CI failure.
- Handling of `prefixParameters` with `Time` quantities
- Don't apply `MeasurementNoise` to TZR TOAs
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
