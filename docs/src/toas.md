# Pulse Times of Arrival (TOAs) 

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

```@docs
TOABase
TOA
DMInfo
WidebandTOA
```

A `TOA` object contains the clock-corrected TOA value in the TDB timescale, uncertainty, and 
observing frequency, along with the solar system ephemeris evaluated at that instance.
The clock corrections from the observatory timescale to TDB and the solar system ephemeris
computations are precomputed using `PINT`, which in turn uses `astropy` underneath. The solar
system ephemerides are represented using the `SolarSystemEphemeris` type.

```@docs
SolarSystemEphemeris
```

The `TOA` objects do not contain any other metadata unlike the `PINT` `TOAs` objects or 
the `tempo2` `observation` objects. Everything that depend on the metadata are precomputed and 
stored elsewhere. This helps with computational efficiency by avoiding repeated string operations.

A performance-critical assumption made in `Vela.jl` is that the TOAs are immutable. This assumption 
is not possible in general pulsar timing packages such as `tempo2` and `PINT` since they allow interactive
removal and flagging of TOAs. `Vela.jl` does not have this use case since it is only meant to do
pulsar timing & noise analysis on TOAs that have already gone through data combination & timing stages.