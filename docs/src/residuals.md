# Residuals

The timing residual corresponding to a TOA can be computed using the `form_residual` function. This function
also returns the DM residual in the case of wideband data.
```@docs
form_residual
```

Under the hood, it calls the `correct_toa` method for each component, and computes the phase residual,
which is the computed phase modulo 1 (this is done using the `phase_residual` function). The timing
residual is the phase residual divided by the topocentric spin frequency (computed using 
`doppler_shifted_spin_frequency`).

The `form_residuals` function computes the residuals for a collection of TOAs.
```@docs
form_residuals
```
