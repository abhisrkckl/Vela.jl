# Explanation

This page provides somewhat technical descriptions about the theory, internal workings, and design 
choices of `Vela.jl`. If you are a beginner please go to Tutorial page first.

## Quantities

It turns out that we can write the entire pulsar timing formula (e.g., see Equations 1-2 of 
[Susobhanan+ 2024](http://doi.org/10.3847/1538-4357/ad59f7)) can be expressed such that all
quantities appearing therein have dimensions of the form \[T^n\]. In practice, this is achieved
by absorbing certain constants appearing in the timing formula into measurable quantities.

For example, `DM -> DMconst*DM`; `M2 -> G*M2/c^3`; `PX -> PX*c/AU`; etc.

This means that we can represent all quantities in pulsar timing in the form `x*s^p` after some scaling,
where `x` is the value of the scaled quantity in SI units, `s` is second, and `p` is an
integer. This is implemented in the [GeometricUnits.jl](https://github.com/abhisrkckl/GeometricUnits.jl/)
package as the `GQ{p,F<:AbstractFloat}` type (with `p ∈ Integer`). This package overloads all the
arithmetic and comparison operators as well as elementary mathematical functions for the `GQ` type.
i.e., `GQ` types can be used just like `Number`s in most places through the magic of mutiple dispatch. 
(`GQ`s do not behave identically to `Number`s in some contexts, so `GQ` is not a subtype of `Number`.)
It also defines iterators and such for the `GQ` type so that we can use it with packages like
`LinearAlgebra.jl`

Note that the dimensionality `p` is a type parameter, which means that the dimensional correctness will
be enforced by the Julia language at "compile time", and it will refuse to execute dimensionally 
incorrect expressions. This provides strong assurances for code correctness. Further, since `p` is 
a type parameter, there is no run time penalty for ensuring dimensional correctness.

## Precision

Pulsar timing is one of the most precise techniques in science. This means that we often deal with
quantities of immense dynamic range which cannot be represented using double-precision floating
point (`Float64`) numbers. The quantities where more than double precision is necessary include the
pulse time of arrivals (TOAs), pulse phases, and pulsar rotational frequencies. 

In `tempo2` and `PINT`, these quantities are represented using the `long double` type available in 
C/C++ (`PINT` uses `numpy.longdouble`, which uses C `long double` internally). Unfortunately,
`sizeof(long double)` is not fixed by the C and C++ standards, and in some hardware it can be 
the same as the `double` type (e.g., the Apple Mx machines). In such cases, `tempo2` falls back to 
the `__float128` type which is available as a compiler extension in `gcc`, whereas `PINT` 
does not work at all.

To avoid this harware dependency, `Vela.jl` represents these quantities using the `Double64` type
available in `DoubleFloats.jl`. This package implements the double-double arithmetic 
([Decker 1971](https://doi.org/10.1007/BF01397083)) which treats an extended-precision number 
as a sum of two `Float64`s. Further, `Double64` is faster than the software-implemented `Float128` 
type from the `Quadmath.jl` package which uses `__float128` under the hood.

Please note that although the core of the `Vela.jl` package should be hardware-independent, its
full functionality won't be available in machines where `PINT` won't work, because it relies on 
`PINT` to do certain one-time computations such as clock corrections.

## Pulse Times of Arrival (TOAs) 

The TOA is the fundamental measurable quantity in *conventional* pulsar timing. They are measured
my folding the time series observations using a known pulsar ephemeris, and then matching the 
resulting integrated pulse profile against a template. This procedure produces a TOA measurement
(usually expressed as an MJD measured against the observatory clock) along with an uncertainty
which is assumed to be Gaussian. The TOA measurements and uncertainties, along with some metadata 
(such as observing frequency, observatory name, observatory reciever and backend information, etc)
is saved in a text file known as a `tim` file. Often, multiple TOAs are measured from the same 
observation by splitting it into multiple frequency sub-bands. This is known as narrowband timing.

Alternatively, we can measure a TOA and a dispersion measure simultaneously from an observation 
using the *wideband* timing method ([Pennucci 2019](http://doi.org/10.3847/1538-4357/aaf6ef)). 
This is useful for dealing with frequency-dependent profile evolution more accurately. It also 
helps reduce the number of TOAs significantly, thus reducing the computational cost of analyzing 
them.
  
In `Vela.jl`, all TOAs are represented by the abstract type `TOABase`. Narrowband TOAs are
represented using the `TOA` type and wideband TOAs are represented by the `WidebandTOA` type.
The latter is a composition of a `TOA` object and a `DMInfo` object that represents a wideband
DM measurement.

A `TOA` object contains the clock-corrected TOA value in the TDB timescale, uncertainty, and 
observing frequency, along with the solar system ephemeris evaluated at that instance.
The clock corrections from the observatory timescale to TBD and the solar system ephemeris
computations are precomputed using `PINT`, which in turn uses `astropy` underneath. The `TOA` objects
do not contain any other metadata unlike the `PINT` `TOAs` objects or the `tempo2` `observation` 
objects. Everything that depend on the metadata are precomputed and stored elsewhere (see below).
This helps with computational efficiency by avoiding repeated string operations.

A performance-critical assumption made in `Vela.jl` is that the TOAs are immutable. This assumption 
is not possible in general pulsar timing packages such as `tempo2` and `PINT` since they allow interactive
removal and flagging of TOAs. `Vela.jl` does not have this use case since it is only meant to do
pulsar timing & noise analysis on TOAs that have already gone through data combination & timing stages.

## The Timing & Noise Model and Timing Residuals

The pulsar timing & noise model (a.k.a. pulsar ephemeris) is a generative mathematical model for the
TOA measurements and uncertainties (and the DM measurements & uncertainties in the wideband paradigm).
The timing residuals are the differences between measured TOAs and the TOAs predicted by the timing model.

The physical and instrumental effects that affect one `TOA` at a time are represented as `Component`s, and
the effects that are correlated across multiple `TOA`s are modeled as a `Kernel`.

There are three types of `Component`s. `DelayComponent`s correct the measured TOA by subtracting certain 
delays. `PhaseComponent`s compute the various contributions to the rotational phase of the pulsar using
delay-corrected TOAs. `WhiteNoiseComponent`s modify the measurement uncertainties in various ways. 

A `TOA` with the various corrections mentioned above is represented as a `CorrectedTOA` (or as a 
`CorrectedWidebandTOA` in the case of a `WidebandTOA`). Each `Component` provides a `correct_toa` method 
which acts on a `CorrectedTOA` / `CorrectedWidebandTOA` using model parameters (represented as 
`NamedTuple{Symbol,GQ}`s) and produces a new object of the same type.

It should be noted that the action of different `Components` do not commute. Therefore, they must be 
applied in the correct order to get sensible results.

The order followed by `Vela.jl` is roughly as follows:

    1. Delay corrections
        a. Solar system effects
        b. Interstellar medium effects
        c. Pulsar binary effects
    2. Phase computation
        a. Pulsar rotational effects
    3. Uncertainty corrections
        a. Measurement noise corrections

These operations finally result in the pulse phase being computed, which is generally not an integer
(implemented in the `phase_residual` method). The phase residual is the pulse phase minus the expected 
pulse number (usually taken to be the integer closest to the pulse phase). The time residual is the 
phase residual divided by the instantaneous topocentric frequency (implemented in the `doppler_shifted_spin_frequency`
method). 

The DM residuals can be similarly computed for `WidebandTOA`s.

A type hierarchy of all available `Component`s is shown below:

```
Component
├─ DelayComponent
│  ├─ ChromaticComponent
│  │  ├─ CMWaveX
│  │  └─ ChromaticTaylor
│  ├─ FrequencyDependent
│  ├─ SolarSystem
│  ├─ Troposphere
│  ├─ BinaryComponent
│  │  ├─ BinaryDDBase
│  │  │  ├─ BinaryDD
│  │  │  ├─ BinaryDDH
│  │  │  ├─ BinaryDDK
│  │  │  └─ BinaryDDS
│  │  └─ BinaryELL1Base
│  │     ├─ BinaryELL1
│  │     └─ BinaryELL1H
│  ├─ DispersionComponent
│  │  ├─ DMWaveX
│  │  ├─ DispersionTaylor
│  │  ├─ SolarWind
│  │  │  └─ SolarWindDispersion
│  │  ├─ DispersionJumpBase
│  │  │  ├─ DispersionJump
│  │  │  └─ ExclusiveDispersionJump
│  │  └─ DispersionOffsetBase
│  │     ├─ DispersionOffset
│  │     └─ ExclusiveDispersionOffset
│  └─ WaveX
├─ PhaseComponent
│  ├─ Glitch
│  ├─ PhaseOffset
│  ├─ Spindown
│  └─ PhaseJumpBase
│     ├─ ExclusivePhaseJump
│     └─ PhaseJump
└─ WhiteNoiseComponent
   ├─ DispersionMeasurementNoise
   └─ MeasurementNoise
```

## The Likelihood Function

The pulsar timing log-likelihood function is given by 

``\ln L = -\frac{1}{2}r^T C^{-1} r - \frac{1}{2}\ln \det C``

where `C` is the covariance matrix that is represented by a `Kernel`, incorporating the measurement 
uncertainties and the various correlated noise processes.

The `calc_lnlike` and `calc_lnlike_serial` functions compute this log-likelihood function. The difference
between them is that the former parallelizes the computation using threads whereas the latter executes
serially. 

The `get_lnlike_func` function return a callable that can be passed on to sampling packages.
It chooses parallel or serial execution based on the number of available threads.

## Correlated Noise and `Kernel`s

The covariance matrix `C` appearing in the log-likelihood expression is represented by a 
`Kernel`. The simplest `Kernel` is a `WhiteNoiseKernel`, which represents the absence  
any correlated noise processes. `C` is diagonal in this case, and its inversion is trivial.

Correlated noise processes in pulsar timing are of two types: time-uncorrelated and 
time-correlated. `ECORR` is a time-uncorrelated noise process that is correlated between 
TOAs obtained from the same observation, but uncorrelated otherwise. The covariance matrix
with only white noise and `ECORR` is block-diagonal, and is represented by the `EcorrKernel`.

It turns out that the log-likelihood can be evaluated at linear time for both the `WhiteNoiseKernel`
and the `EcorrKernel`.

(Time-correlated noise processes are currently work in progress...)

## `Parameter`s

Each `TimingModel` has a number of model parameters, some of which are free and some frozen.
Some of these parameters are single parameters, like the phase offset `PHOFF`, and is 
represented by the `Parameter` type. Some come as families of similar parameters, such as the
pulsar frequency `F0` and its derivatives `F1`, `F2`, etc. These are represented as 
`MultiParameter`s. These contain information about parameter names, their original units (used
by `PINT`), the scaling factors for converting to and from the `Vela.jl` units, default values,
 whether they are free/frozen, etc. 

The `correct_toa` methods expect the parameters to be passed in as a `NamedTuple` containing 
`GQ`s for single `Parameter`s and `NTuple`s of `GQ`s for `MultiParameter`s. 

The `ParamHandler` type contains a collection of all the parameters of a model, and is used
to convert a `Vector{Float64}` of free parameter values provided by the sampler to a `NamedTuple` 
the `correct_toa` methods can understand. Each `TimingModel` has a `ParamHandler` attached to it.

## Prior and Posterior distributions

The prior distributions of each free parameter is represented as a `Prior` object. These use
the `Distribution`s defined by `Distributions.jl` under the hood. Each `Prior` has the 
`lnprior` and `prior_transform` methods which compute the log-prior distribution and the 
[prior transform function](http://kylebarbary.com/nestle/prior.html) for that parameter.
The former is necessary for MCMC samplers and the latter for nested samplers.
Please note that these functions act on `Float64`s rather than `GQ`s because the samplers 
only provide `Float64`s.

A `TimingModel` also has the `lnprior` and `prior_transform` methods; they evaluate the log-prior
and the prior transform over all free parameters.

The log-posterior is the sum of the log-likelihood and the log-prior up to an additive 
constant, and can be evaluated using the `lnpost` function. The `get_lnpost_func` function
returns a callable that can be passed on to samplers. 

Note that the expensive log-likelihood is evaluated only if the log-prior is finite. 

