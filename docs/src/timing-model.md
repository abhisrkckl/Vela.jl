# The Timing & Noise Model

The pulsar timing & noise model (a.k.a. pulsar ephemeris) is a generative mathematical model for the
TOA measurements and uncertainties (and the DM measurements & uncertainties in the wideband paradigm).
The timing residuals are the differences between measured TOAs and the TOAs predicted by the timing model.
See [Hobbs+ 2006](https://doi.org/10.1111/j.1365-2966.2006.10302.x) and 
[Edwards+ 2006](https://doi.org/10.1111/j.1365-2966.2006.10870.x) for details.

The physical and instrumental effects that affect one `TOA` at a time are represented as `Component`s, and
the effects that are correlated across multiple `TOA`s are modeled as a `Kernel`.

```@docs
Component
```

Each component has a `correct_toa` method which produces a TOA correction.
```@docs
correct_toa
```

`Kernel`s will be discussed in its own section.

## TOA corrections

Such a correction can be one or more of the following:
    1. A delay that modifies the TOA value
    2. A phase correction that modifies the TOA phase
    3. A Doppler correction that modifies the observing frequency and the pulsar spin frequency
    4. An EFAC or EQUAD that modifies the TOA uncertainty
    5. A DM correction that modifies the wideband DM measurement
    6. A DMEFAC or DMEQUAD that modifies the wideband DM uncertainty

The accumlated TOA corrections are represented by the `TOACorrection` and `WidebandTOACorrection`
types, which are derived from `TOACorrectionBase`.
```@docs
TOACorrectionBase
TOACorrection
WidebandTOACorrection
```

The following methods extract some of the intermediate-corrected quantities of interest.
```@docs
corrected_toa_value
doppler_corrected_observing_frequency
doppler_shifted_spin_frequency
scaled_toa_error_sqr
scaled_dm_error_sqr
phase_residual
dm_residual
```

There are three types of `Component`s as shown below based on what type of corrections they
provide.
```@eval
using InteractiveUtils
using AbstractTrees
using Vela
using Markdown

AbstractTrees.children(d::DataType) = subtypes(d)
Markdown.MD(Markdown.Code(repr_tree(Component, maxdepth=1)))
```

## Delay components

`DelayComponent`s correct the measured TOA by subtracting certain delays.

```@docs
DelayComponent
```

`DelayComponent` has further subtypes which represent different physical processes producing delays.
```@eval
using InteractiveUtils
using AbstractTrees
using Vela
using Markdown

AbstractTrees.children(d::DataType) = subtypes(d)
Markdown.MD(Markdown.Code(repr_tree(DelayComponent, maxdepth=1)))
```

### Solar system delays
```@docs
SolarSystem
```

In addition to a delay, `SolarSystem` also producess a Doppler correction which applies to the
observing frequency and the pulsar spin frequency. 

This component barycenters the TOA. It will skip TOAs that are already barycentered, e.g., TOAs
measured using space-based telescopes. The `is_barycentered` function checks whether a TOA has been
barycentered.
```@docs
is_barycentered
```

### Dispersion delays
`DispersionComponent`s represent the dispersion of the radio waves due to the free electrons present
along the line of sight to the pulsar. This may include the ionized interstellar medium as well as 
solar wind.

```@docs
DispersionComponent
```

`DispersionComponent`s provide a delay and a DM correction (for wideband TOAs), which are related by
the equation ``\Delta_{\text{DM}} = K * \text{DM} / \nu^2``. The dispersion correction is computed via 
the `dispersion_slope` function.
```@docs
dispersion_slope
```

The `DispersionComponent`s available in `Vela.jl` are
```@eval
using InteractiveUtils
using AbstractTrees
using Vela
using Markdown

AbstractTrees.children(d::DataType) = subtypes(d)
Markdown.MD(Markdown.Code(repr_tree(DispersionComponent)))
```

The most basic `DispersionComponent` is `DispersionTaylor`.
```@docs
DispersionTaylor
```

`DMWaveX` and `PowerlawDispersionNoiseGP` provide two representations of the
dispersion noise / stochastic DM variations.
```@docs
DMWaveX
PowerlawDispersionNoiseGP
```

`SolarWind` is a simple model for solar wind dispersion.
```@docs
SolarWindDispersion
```

`DispersionJump` and `DispersionOffset` represent two types of system-dependent dispersion 
offsets.
```@docs
DispersionJump
DispersionOffset
```

### Chromatic delays
Chromatic delays are similar to dispersion delays, but have a different powerlaw dependence 
on the observing frequency. Such delays can occur due to interstellar scattering or 
frequency-dependent dispersion. A chromatic delay is given by ``\Delta_{\text{CM}} = K * \text{CM} / \nu^\alpha``
where CM is called the chromatic measure, and ``\alpha`` is called the chromatic index.
The effect of chromatic effects on wideband TOAs is not well-understood.

```@eval
using InteractiveUtils
using AbstractTrees
using Vela
using Markdown

AbstractTrees.children(d::DataType) = subtypes(d)
Markdown.MD(Markdown.Code(repr_tree(ChromaticComponent)))
```

```@docs
ChromaticTaylor
CMWaveX
PowerlawChromaticNoiseGP
```

### Binary delays
Similar to solar system delays, the binary motion of the pulsar also introduces various 
delays to the TOAs, including Rømer delay, Shapiro delay, and Einstein delay. 
```@docs
BinaryComponent
```

Different binary models are used based on the properties of the binary orbit.
```@eval
using InteractiveUtils
using AbstractTrees
using Vela
using Markdown

AbstractTrees.children(d::DataType) = subtypes(d)
Markdown.MD(Markdown.Code(repr_tree(BinaryComponent)))
```

`Vela.jl` has two families of binary models. The `ELL1` family is used for nearly circular binaries
and the `DD` family is used for eccentric orbits. The different models are characterized by their treatment
of Shapiro delay, Kopeikin corrections, etc.
```@docs
BinaryDD
BinaryDDH
BinaryDDK
BinaryDDS
BinaryELL1
BinaryELL1H
```

## Phase components
A `PhaseComponent` contributes to the phase computation from a delay-corrected TOA.
```@docs
PhaseComponent
```

The hierarchy of `PhaseComponent`s is given below.
```@eval
using InteractiveUtils
using AbstractTrees
using Vela
using Markdown

AbstractTrees.children(d::DataType) = subtypes(d)
Markdown.MD(Markdown.Code(repr_tree(PhaseComponent)))
```

```@docs
Spindown
Glitch
PhaseOffset
PhaseJump
```

## White noise components
These components modify the TOA or DM uncertainty in some manner.

```@eval
using InteractiveUtils
using AbstractTrees
using Vela
using Markdown

AbstractTrees.children(d::DataType) = subtypes(d)
Markdown.MD(Markdown.Code(repr_tree(WhiteNoiseComponent)))
```

```@docs
MeasurementNoise
DispersionMeasurementNoise
```

## Order of components
It should be noted that the action of different `Components` do not commute in general. Therefore, 
they must be applied in the correct order to get sensible results. The order followed by `Vela.jl` 
is roughly as follows:

    1. Delay corrections
        a. Solar system effects
        b. Interstellar medium effects
        c. Pulsar binary effects
    2. Phase computation
        a. Pulsar rotational effects
    3. Uncertainty corrections
        a. Measurement noise corrections

## The TZR TOA
The pulse phases are measured with respect to a fictitious fiducial TOA called the TZR TOA.
This is represented using the `TOA` class, but is distinguished from physical TOAs using the 
`tzr` attribute. The `make_tzr_toa` function creates a TZR TOA and the `is_tzr` function checks 
whether a TOA is a TZR TOA.
```@docs
make_tzr_toa
is_tzr
```

## The `TimingModel` type
A timing & noise model is represented by the `TimingModel` type. 
```@docs
TimingModel
```

It has the following contents:
    1. Pulsar name (PSR), solar system ephemeris name (EPHEM), name of the TT timescale realization (CLOCK), etc.
    2. An ordered collection of `Components`
    3. A `Kernel`
    4. A `ParamHandler` containing information about model parameters.
    5. An ordered collection of `Prior`s for each free parameter

Some of these are explained in the following sections.