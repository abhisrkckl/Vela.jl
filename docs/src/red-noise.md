# Red noise models

A time-correlated stochastic process whose power spectral density is a decreasing function of
the conjugate frequency (not to be confused with observing frequency) is known as red noise.
Depending on the pulsar, different types of red noise processes will be present in the TOAs, including
spin noise, dispersion noise, and the stochastic gravitational wave background.

Such a process is usually modeled as truncated Fourier series.

``\Delta(t) = (\frac{\nu_{\text{ref}}}{\nu})^\alpha \sum_{j=1}^N \{ a_j\cos(2\pi j f_1(t-t_0)) + b_j\sin(2\pi j f_1(t-t_0)) \}\,.``

where ``\nu`` is the observing frequency, ``\nu_{\text{ref}}`` is a reference frequency, ``\alpha`` is the
chromatic index, ``N`` is the number of harmonics, ``f_1`` is the fundamental frequency, ``t_0`` is a fiducial 
time, and ``a_j`` and ``b_j`` are Fourier coefficients.

Two types of red noise models are currently available in `Vela.jl`. 

The `WaveX` family of models treat the Fourier coefficients as unconstrained free parameters 
(e.g., `WXSIN_` and `WXCOS_`).

The `PowerlawRedNoiseGP` family of models treat these coefficients as Gaussian random variables whose
variances, interpreted as power spectral densities, follow a power law spectrum.

``\left\langle a_j \right\rangle = \left\langle b_j \right\rangle = 0``

``\left\langle a_j a_k \right\rangle = \left\langle b_j b_k \right\rangle = \sigma_j^2 \delta_{jk}``

``\left\langle a_j b_k \right\rangle = 0``

``\sigma_j = P(f_j) = \frac{A^2}{12\pi^2 f_{\text{yr}}^3} f_1 \left(\frac{f_{\text{yr}}}{f}\right)^\gamma``

where ``A`` is the powerlaw amplitude ``\gamma`` is the powerlaw index, and ``f_{\text{yr}}=1 \text{yr}^{-1}``.
The prior parameters ``A`` and ``\gamma`` are also treated as free parameters and sampled over. However, in this
case, the geometry of the parameter space exhibits [Neil's funnel](https://crackedbassoon.com/writing/funneling)-like
geometry, which is hard for MCMC samplers to deal with. To avoid this, we use ``\bar{a}_j=a_j/\sigma_j`` and 
``\bar{b}_j=b_j/\sigma_j`` as free parameters where 

``\left\langle\bar{a}_j^2\right\rangle = \left\langle\bar{b}_j^2\right\rangle = 1\,.``

In `PowerlawRedNoiseGP`, this is represented by the parameters `PLREDSIN_` and `PLREDCOS_`. The powerlaw parameters
are `PLREDAMP` and `PNREDGAM`. Please note that the `PINT`-format `par` files do not support `PLREDSIN_` and `PLREDCOS_`.
Hence, while they are included in the posterior samples, they will not be included in the output par file.  

Similar representations also exist for dispersion noise (``\alpha=2``) and chromatic noise (see [Dispersion delays](@ref)
and [Chromatic delays](@ref)).

Other types of spectral models such as free spectrum, t-process, broken powerlaw, running power law, etc are not
yet implemented.