### LSST-Y1

| DESC-CCL (reference settings) vs CoCoA | 3x2pt | shear | $\gamma_t$ | $w(\theta)$ | points left out |
|---|---|---|---|---|---|
| fiducial | 8.757 | 1.3e-03 | 7.290 | 6.259 | 0 |
| $\Omega_m = 0.25$ | 6.559 | 4.7e-04 | 4.869 | 5.477 | 0 |
| $\Omega_m = 0.35$ | 10.502 | 3.2e-03 | 9.442 | 6.459 | 0 |
| $n_s = 0.92$ | 9.070 | 1.3e-03 | 7.563 | 6.547 | 0 |
| $n_s = 1.01$ | 8.454 | 1.4e-03 | 7.057 | 5.972 | 0 |

| Fiducial, one choice at a time | 3x2pt | shear | $\gamma_t$ | $w(\theta)$ | points left out |
|---|---|---|---|---|---|
| CoCoA default vs CoCoA high accuracy | 4.0e-03 | 4.9e-04 | 3.8e-04 | 3.7e-03 | 0 |
| DESC-CCL default FKEM/$\ell$ sampling vs reference | 0.345 | 1.3e-06 | 0.027 | 0.395 | 0 |
| both in Limber: DESC-CCL vs CoCoA | 0.158 | 1.3e-03 | 1.1e-03 | 0.157 | 0 |
| non-Limber effect in CoCoA | 159.553 | 0 | 1.979 | 151.434 | 0 |
| non-Limber effect in DESC-CCL | 149.109 | 0 | 11.212 | 134.644 | 0 |
| non-Limber effect in DESC-CCL, separable $P_{\rm lin}$ (diagnostic) | 144.235 | 0 | 2.063 | 136.712 | 0 |
| DESC-CCL, separable $P_{\rm lin}$ (diagnostic), vs CoCoA | 0.208 | 1.3e-03 | 0.064 | 0.165 | 0 |
| RSD in $\gamma_t$: effect in DESC-CCL | 0.358 | 0 | 0.358 | 0 | 0 |
| DESC-CCL $C_{gs}$ in Limber vs CoCoA | 9.228 | 1.3e-03 | 2.017 | 6.259 | 0 |
| DESC-CCL no RSD vs CoCoA | 115.402 | 1.3e-03 | 7.290 | 117.499 | 0 |
| DESC-CCL RSD also in $\gamma_t$ vs CoCoA | 8.703 | 1.3e-03 | 6.690 | 6.259 | 0 |
| DESC-CCL flat-sky FFTLog, bin centers vs CoCoA | 12.684 | 1.818 | 7.710 | 10.871 | 0 |
| DESC-CCL full-sky, bin centers vs CoCoA | 10.311 | 2.003 | 3.586 | 8.656 | 0 |
| DESC-CCL the CCL-benchmark modeling vs CoCoA | 153.590 | 5.692 | 14.694 | 166.482 | 0 |

| DESC-CCL run | pair | type | call | error |
|---|---|---|---|---|
| fid_ccl_numerics | l4-s0 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| fid_flat | l4-s0 | NG | correlation | Error CCL_ERROR_SPLINE: ccl_correlation.c: ccl_tracer_corr_fftlog(): failed to create spline |
| fid_norsd | l4-s0 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| fid_points | l4-s0 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| fid_ref | l4-s0 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| fid_rsd_gs | l4-s0 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| fid_separable | l4-s0 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| ns_hi_ref | l4-s0 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| ns_lo_ref | l4-s0 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| omm_hi_ref | l4-s0 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| omm_lo_ref | l4-s0 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |

### Roman-Real

| DESC-CCL (reference settings) vs CoCoA | 3x2pt | shear | $\gamma_t$ | $w(\theta)$ | points left out |
|---|---|---|---|---|---|
| fiducial | 0.228 | 0.058 | 0.269 | 0.096 | 13 |
| $\Omega_m = 0.25$ | 0.166 | 0.023 | 0.193 | 0.087 | 13 |
| $\Omega_m = 0.35$ | 0.291 | 0.123 | 0.360 | 0.098 | 13 |
| $n_s = 0.92$ | 0.236 | 0.055 | 0.280 | 0.099 | 13 |
| $n_s = 1.01$ | 0.215 | 0.055 | 0.259 | 0.091 | 13 |

| Fiducial, one choice at a time | 3x2pt | shear | $\gamma_t$ | $w(\theta)$ | points left out |
|---|---|---|---|---|---|
| CoCoA default vs CoCoA high accuracy | 7.3e-03 | 1.9e-03 | 3.4e-03 | 6.4e-03 | 0 |
| DESC-CCL default FKEM/$\ell$ sampling vs reference | 4.7e-04 | 2.5e-06 | 1.7e-04 | 3.4e-04 | 13 |
| both in Limber: DESC-CCL vs CoCoA | 0.031 | 0.058 | 0.060 | 0.011 | 0 |
| non-Limber effect in CoCoA | 9.845 | 0 | 0.488 | 8.706 | 0 |
| non-Limber effect in DESC-CCL | 9.293 | 0 | 0.791 | 7.910 | 13 |
| non-Limber effect in DESC-CCL, separable $P_{\rm lin}$ (diagnostic) | 9.263 | 0 | 0.469 | 8.189 | 13 |
| DESC-CCL, separable $P_{\rm lin}$ (diagnostic), vs CoCoA | 0.062 | 0.058 | 0.102 | 0.010 | 13 |
| RSD in $\gamma_t$: effect in DESC-CCL | 0.076 | 0 | 0.076 | 0 | 13 |
| DESC-CCL $C_{gs}$ in Limber vs CoCoA | 0.708 | 0.058 | 0.580 | 0.096 | 0 |
| DESC-CCL no RSD vs CoCoA | 6.982 | 0.058 | 0.269 | 7.207 | 13 |
| DESC-CCL RSD also in $\gamma_t$ vs CoCoA | 0.258 | 0.058 | 0.253 | 0.096 | 13 |
| DESC-CCL flat-sky FFTLog, bin centers vs CoCoA | 22.602 | 19.417 | 58.273 | 13.693 | 13 |
| DESC-CCL full-sky, bin centers vs CoCoA | 22.771 | 28.309 | 84.600 | 13.843 | 13 |
| DESC-CCL the CCL-benchmark modeling vs CoCoA | 55.938 | 13.128 | 41.945 | 64.066 | 0 |

| DESC-CCL run | pair | type | call | error |
|---|---|---|---|---|
| fid_ccl_numerics | l7-s2 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| fid_flat | l7-s2 | NG | correlation | Error CCL_ERROR_SPLINE: ccl_correlation.c: ccl_tracer_corr_fftlog(): failed to create spline |
| fid_norsd | l7-s2 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| fid_points | l7-s2 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| fid_ref | l7-s2 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| fid_rsd_gs | l7-s2 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| fid_separable | l7-s2 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| ns_hi_ref | l7-s2 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| ns_lo_ref | l7-s2 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| omm_hi_ref | l7-s2 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |
| omm_lo_ref | l7-s2 | NG | correlation | Error CCL_ERROR_MEMORY: ccl_correlation.c: ccl_tracer_corr_legendre(): ran out of memory |

