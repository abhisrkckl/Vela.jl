# The Likelihood Function and Kernels

The pulsar timing log-likelihood function is given by 

``\ln L = -\frac{1}{2}r^T C^{-1} r - \frac{1}{2}\ln \det C``

where `C` is the covariance matrix that is represented by a `Kernel`, incorporating the measurement 
uncertainties and the various correlated noise processes.

The `calc_lnlike` and `calc_lnlike_serial` functions compute this log-likelihood function. The difference
between them is that the former parallelizes the computation using threads over `TOA`s whereas the latter 
executes serially. 
```@docs
calc_lnlike
calc_lnlike_serial
```

The `get_lnlike_func` function return a callable that can be passed on to sampling packages.
It chooses parallel or serial execution based on the number of available threads. It also has a
`vectorize` option that evaluates the likelihood function over multiple sets of parameters parallely.
This is supported by some samplers like `emcee`, and is more efficient than parallelizing over `TOA`s.
There is also a `get_lnlike_serial_func` function that always returns the serial version of the callable.
```@docs
get_lnlike_func
get_lnlike_serial_func
```

## Kernels
The matrix operations appearing in the likelihood function expression are evaluated with
the help of `Kernel` objects.
```@docs
Kernel
```

Two types of `Kernel`s are currently supported.
```@eval
using InteractiveUtils
using AbstractTrees
using Vela
using Markdown

AbstractTrees.children(d::DataType) = subtypes(d)
Markdown.MD(Markdown.Code(repr_tree(Kernel)))
```

```@docs
WhiteNoiseKernel
EcorrKernel
```

More complex `Kernel`s supporting more general forms of ``C`` (including time-correlated noise) will be 
implemented in the future. Currently such processes are treated as `Components`, e.g., `PLRedNoiseGP`.