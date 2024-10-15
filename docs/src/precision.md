# Precision

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