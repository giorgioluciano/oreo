# oreo 1.1 (2026-09-28)

## Bug fixes

* `rpp_fft()`: the third time derivative of the Fourier-reconstructed strain and
  strain-rate coordinates used `A cos - B sin` instead of the analytical
  `A sin - B cos`. This affected only `Gp_t_dot`, `Gpp_t_dot` and `G_speed`;
  `Gp_t`, `Gpp_t`, `G_star_t`, `delta_t` and the T/N/B vectors were not affected.
  The same expression is present in the MATLAB SPPplus v2 code from which oreo
  was translated.
* `rpp_fft()` and `Rpp_num()`: `disp_stress` (and therefore `eq_strain_est`)
  mixed the normalized and raw strain-rate conventions (an extra `1/omega`
  factor). It is now consistent in both conventions.

## Changes

* New argument `norm_rate = TRUE` in `rpp_fft()` and `Rpp_num()`. The trajectory
  now uses the published normalized coordinate `[gamma, gamma_dot/omega, sigma]`,
  so `Gpp_t` is returned in Pa. In oreo 1.0 `Gpp_t` was `G''_t/omega`.
  Use `norm_rate = FALSE` to reproduce the oreo 1.0 convention.
* New self-test: `source(system.file("tests", "oreo_selftest.R", package = "oreo"))`
  checks the implementation against exact analytical results (54 checks,
  including phase-shifted strain and multi-cycle input).
* `Rpp_plot_v3.r` restored to the CRAN 1.0 version (the GitHub copy contained a
  syntax error that prevented installation from source).

## Acknowledgement

Thanks to Guo Huang (Northeast Agricultural University) for the detailed report
and independent verification.
