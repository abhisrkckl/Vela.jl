# Quantities

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

```@docs
GQ
```

Note that the dimensionality `p` is a type parameter, which means that the dimensional correctness will
be enforced by the Julia language at "compile time", and it will refuse to execute dimensionally 
incorrect expressions. This provides strong assurances for code correctness. Further, since `p` is 
a type parameter, there is no run time penalty for ensuring dimensional correctness.