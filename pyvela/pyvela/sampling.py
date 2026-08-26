import warnings
import numpy as np
import emcee
from .spnta import SPNTA


def get_start_samples(spnta: SPNTA, s: float, nwalkers: int) -> np.ndarray:
    """Get starting samples for the MCMC. nwalkers is the number of samples
    to be returned."""

    p0_ = spnta.draw_from_prior(nwalkers)
    p0 = (
        ((1 - s) * spnta.default_params + s * p0_)
        if np.isfinite(spnta.lnpost(spnta.default_params))
        else p0_
    )
    return p0


def run_emcee_until_converged(
    sampler,
    initial_state,
    max_steps=100_000,
    initial_chunk=500,
    min_chunk=100,
    max_chunk=5000,
    tau_multiplier=50,
    tau_rtol=0.05,
    discard_for_tau=400,
    min_checks=2,
    progress=True,
):
    """
    Run an emcee chain until convergence, using an adaptive chunk size.

    Convergence requires:
        1. N > tau_multiplier * max(tau)
        2. The estimated autocorrelation time is stable to tau_rtol.
        3. The stability criterion has been satisfied for min_checks consecutive convergence checks.

    Parameters
    ----------
    log_prob_fn : callable
        Log-posterior function.

    initial_state : array-like
        Initial walker positions, shape (nwalkers, ndim).
        Ignored when resuming from a non-empty backend.

    max_steps : int
        Maximum total number of MCMC steps.

    initial_chunk : int
        Number of steps in the first chunk.

    min_chunk, max_chunk : int
        Bounds on the adaptive chunk size.

    tau_multiplier : float
        Require chain length N > tau_multiplier * max(tau).

    tau_rtol : float
        Maximum relative change in tau between convergence checks.

    min_checks : int
        Number of consecutive successful stability checks required.

    backend : emcee.backends.Backend, optional
        If supplied, allows the chain to be resumed.

    progress : bool
        Show emcee progress bar.

    Returns
    -------
    sampler : emcee.EnsembleSampler
        The sampler containing the complete chain.

    converged : bool
        Whether the convergence criterion was reached.

    tau : ndarray or None
        Final estimated autocorrelation times.
    """

    initial_state = np.asarray(initial_state)

    start = initial_state
    total_steps = 0

    tau_old = None
    stable_checks = 0
    tau = None

    chunk = initial_chunk

    while total_steps < max_steps:

        # Don't exceed max_steps.
        chunk = min(chunk, max_steps - total_steps)

        sampler.run_mcmc(
            start,
            chunk,
            progress=progress,
        )

        # Subsequent calls continue from the current sampler state.
        start = sampler.get_chain()[-1, :, :]

        total_steps = sampler.iteration

        # Need enough samples before tau estimation is meaningful.
        if total_steps < 50:
            continue

        try:
            tau = sampler.get_autocorr_time(tol=0, discard=discard_for_tau)
        except emcee.autocorr.AutocorrError:
            print(f"N = {total_steps}: " "chain too short for reliable tau estimate")

            # Increase chunk size gradually.
            chunk = min(
                max(chunk * 2, min_chunk),
                max_chunk,
            )
            continue

        tau_max = np.max(tau)

        # ------------------------------------------------------------
        # Criterion 1: chain is sufficiently long
        # ------------------------------------------------------------

        long_enough = total_steps > tau_multiplier * tau_max

        # ------------------------------------------------------------
        # Criterion 2: autocorrelation time has stabilized
        # ------------------------------------------------------------

        if tau_old is None:
            tau_stable = False
            stable_checks = 0
        else:
            relative_change = np.abs(tau - tau_old) / tau_old

            tau_stable = np.all(relative_change < tau_rtol)

            if tau_stable:
                stable_checks += 1
            else:
                stable_checks = 0

        # ------------------------------------------------------------
        # Diagnostics
        # ------------------------------------------------------------

        print(
            f"N = {total_steps:6d} | "
            f"tau_max = {tau_max:8.2f} | "
            f"N/tau = {total_steps / tau_max:7.1f} | "
            f"long = {long_enough} | "
            f"stable = {tau_stable} | "
            # f"stable_checks = {stable_checks}"
        )

        # ------------------------------------------------------------
        # Convergence
        # ------------------------------------------------------------

        if long_enough and stable_checks >= min_checks:
            print("Converged.")
            return True, tau

        tau_old = tau.copy()

        # ------------------------------------------------------------
        # Adaptive chunk size
        #
        # We want the next chunk to be a fraction of tau.
        # ------------------------------------------------------------

        target_chunk = int(0.5 * tau_max) if np.isfinite(tau_max) else initial_chunk

        chunk = np.clip(
            target_chunk,
            min_chunk,
            max_chunk,
        )

    warnings.warn("Maximum number of steps reached without convergence.")

    return False, tau
