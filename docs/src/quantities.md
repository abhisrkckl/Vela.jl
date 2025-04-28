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

`GeometricUnits.jl` also implements the following operations.

1. Unary `+` and `-` operators for all `GQ` types
2. Binary `+` and `-` operators for `GQ` types with like dimensions
3. `*`, `/` operators for all `GQ` types
4. `^` operator for various cases where the output is a valid `GQ`
5. `sqrt` `cbrt`, `root` functions for cases where the output is a valid `GQ`
6. `==`, `!=`, `<`, `<=`, `>`, `>=`, `≈` operators for `GQ` types with like dimensions
7. Trigonometric functions (`sin`, `cos`, `sincos`, `tan`, `csc`, `sec`, `cot`) for dimensionless inputs
8. Inverse trigonometric functions (`asin`, `acos`, `atan`, `acsc`, `asec`, `acot`) for dimensionless inputs 
9. `atan` function for a pair of  `GQ`s types with like dimensions
10. `exp`, `exp2`, `exp10`, `log`, `log2`, `log10` functions for dimensionless inputs
11. `abs`, `sign`, `floor`, `ceil` functions for all `GQ` types
12. `isfinite` and `isnan` functions for all `GQ` types
13. `zero` and `oneunit` functions for all `GQ` types
14. `taylor_horner` and `taylor_horner_integral` functions