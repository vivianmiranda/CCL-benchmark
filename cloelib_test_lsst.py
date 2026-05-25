"""cloelib analog of the LSST Y1 CCL 3x2pt benchmark.

What's timed: only the cloelib pipeline (tracer construction, get_Cl,
get_xi, flattening). CAMBBackground + CAMBNonLinearPerturbations
construction is done OUTSIDE the timed region.

JAX is pinned to CPU (no GPU usage) via JAX_PLATFORMS=cpu set before
import.

═════════════════════════════════════════════════════════════════════════
COMPAT NOTES — what's missing relative to the LSST CCL script
═════════════════════════════════════════════════════════════════════════

1. IA model: TATT → NLA.
   cloelib's ShearTracer only supports the NLA model via (AIA, CIA, EtaIA).
   There is no `translate_IA_norm`, no `PTIntrinsicAlignmentTracer`, and
   no pathway for A2/eta2 (tidal-torquing) or A1δ (density-tidal) terms.
   We map A1→AIA, eta1→EtaIA; A2, eta2 are dropped.

2. Galaxy clustering w(θ): non-Limber FKEM → Limber.
   cloelib's AngularTwoPoint.get_Cl is Limber-only. The CCL LSST script
   used FKEM with l_limber=100, fkem_Nchi=500 for the lens autos.

3. ξ projection: flat-sky FFTLog → curved-sky Wigner-d.
   The LSST CCL script uses `method='FFTLog'` (flat-sky, point evaluation
   at scalar θ — no `theta_max`, no bin averaging). cloelib uses
   `AngularCorrelationFunctionWigner` (curved-sky Wigner small-d, also
   point evaluation). The geometries differ by O((θ/rad)²) at large θ —
   relevant for θ up to 900' (=15°) at the wide end.

4. Cosmology parameterisation: σ8 → As.
   CAMBBackground takes As, not σ8. Fiducial As=2.1e-9 anchors σ8≈0.84.
   We perturb As by 10% rather than σ8 by 5% (σ8 ∝ √As).

5. No PT pipeline.
   Linear bias only, through cloelib's `galaxy_bias_model='per_bin'`.

═════════════════════════════════════════════════════════════════════════
"""

import os
import time
import numpy as np

# Pin JAX to CPU + float64 BEFORE importing jax.
os.environ["JAX_PLATFORMS"] = "cpu"
os.environ["JAX_ENABLE_X64"] = "True"

import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

assert jax.default_backend() == "cpu", (
  f"JAX backend is {jax.default_backend()}, not CPU — "
  "check JAX_PLATFORMS env var or jaxlib install"
)

from cloelib.cosmology.camb_cosmology import (
  CAMBBackground,
  CAMBNonLinearPerturbations,
)
from cloelib.observables.photo import ShearTracer, PositionsTracer
from cloelib.summary_statistics.angular_two_point import AngularTwoPoint
from cloelib.summary_statistics.angular_correlation_function_wigner import (
  AngularCorrelationFunctionWigner,
)

# ── Configuration ──────────────────────────────────────────────
nz_lens_file = "lsst_y1_lens.nz"
nz_source_file = "lsst_y1_source.nz"
lens_ntomo = 5
source_ntomo = 5
n_theta = 26
theta_min_arcmin = 2.5
theta_max_arcmin = 900.0

# No GGL exclusions in the LSST setup
ggl_excl_set = set()

n_warmup = 1
n_bench = 10
rng_seed = 42

# Fiducial cosmology (CAMBBackground takes As, not sigma8)
fid_params = dict(
  H0=67.0, Omega_b0=0.045, Omega_cdm0=0.27, Omega_k0=0.0,
  As=2.1e-9, ns=0.96, mnu=0.0, N_mnu=0,
  w0=-1.0, wa=0.0, gamma_MG=0.545,
)

# Fiducial NLA (mapped from TATT: A1→AIA, eta1→EtaIA; A2/eta2 dropped)
fid_nla = dict(AIA=0.7, EtaIA=-1.7)
CIA = 0.0134

fid_bias = {f'b{i}': 1.2 + 0.1 * i for i in range(lens_ntomo)}
dz_sigma = 0.01

cosmo_perturb_scales = dict(
  H0=0.05, Omega_b0=0.05, Omega_cdm0=0.05, As=0.10, ns=0.05,
)

# ── Read n(z) once ────────────────────────────────────────────
def load_nz(filename, ntomo):
  data = np.loadtxt(filename)
  z = data[:, 0]
  nz = np.stack([data[:, i + 1] for i in range(ntomo)])
  return z, nz

z_lens_raw, nz_lens_raw = load_nz(nz_lens_file, lens_ntomo)
z_source_raw, nz_source_raw = load_nz(nz_source_file, source_ntomo)

# ── Common equispaced z-grid ──────────────────────────────────
# AngularTwoPoint.get_Cl does dz = z[1] - z[0] → must be equispaced and
# must not contain z=0. CAMB caps at 256 redshifts internally (z=0 is
# auto-prepended) so keep N_z ≤ 255.
N_z = 250
z_grid_np = np.linspace(0.005, 3.0, N_z)
z_grid_j = jnp.asarray(z_grid_np)

def interp_normalise(z_src, nz_src, z_target):
  ntomo = nz_src.shape[0]
  out = np.zeros((ntomo, len(z_target)))
  for i in range(ntomo):
    out[i] = np.interp(z_target, z_src, nz_src[i], left=0.0, right=0.0)
    norm = np.trapezoid(out[i], z_target)
    if norm > 0:
      out[i] = out[i] / norm
  return out

nz_lens = interp_normalise(z_lens_raw, nz_lens_raw, z_grid_np)
nz_source = interp_normalise(z_source_raw, nz_source_raw, z_grid_np)
dndz_l = jnp.asarray(nz_lens)
dndz_s = jnp.asarray(nz_source)

# ── ell, theta, k grids ───────────────────────────────────────
# Match the LSST CCL script's 3-segment ell construction
ell_low = np.arange(2, 50)
ell_mid = np.unique(np.geomspace(50, 3000, 150).astype(int))
ell_high = np.unique(np.geomspace(3000, 30000, 150).astype(int))
ell = np.unique(np.concatenate([ell_low, ell_mid, ell_high]))
ell_j = jnp.asarray(ell)

# theta as direct points (not edges) — matches the CCL FFTLog call
theta_arcmin = np.logspace(
  np.log10(theta_min_arcmin), np.log10(theta_max_arcmin), n_theta,
)
theta_rad = jnp.deg2rad(jnp.asarray(theta_arcmin / 60.0))

ks_j = jnp.logspace(-4, 2, 256)

# ── CAMB setup (NOT timed) ────────────────────────────────────
def make_perturbations(params):
  bg = CAMBBackground(**params)
  pert = CAMBNonLinearPerturbations(
    background=bg, linearperturbations=None, redshifts=z_grid_np,
  )
  return pert

# ── cloelib 3x2pt evaluation (TIMED) ──────────────────────────
def compute_3x2pt_cloelib(pert, nla, bias, dz_lens, dz_source):
  shear_nuisance = dict(AIA=nla['AIA'], CIA=CIA, EtaIA=nla['EtaIA'])
  for i in range(source_ntomo):
    shear_nuisance[f'multiplicative_bias_{i+1}'] = 0.0
    shear_nuisance[f'dz_shear_{i+1}'] = float(dz_source[i])
    shear_nuisance[f'width_shear_{i+1}'] = 1.0

  pos_nuisance = {}
  for i in range(lens_ntomo):
    pos_nuisance[f'magnification_bias_{i+1}'] = 0.0
    pos_nuisance[f'dz_pos_{i+1}'] = float(dz_lens[i])
    pos_nuisance[f'width_pos_{i+1}'] = 1.0
    pos_nuisance[f'b1_photo_bin{i}'] = float(bias[f'b{i}'])

  shear_tracer = ShearTracer(
    perturbations=pert, dndz=dndz_s, z=z_grid_j, nuisance_params=shear_nuisance,
  )
  pos_tracer = PositionsTracer(
    perturbations=pert, dndz=dndz_l, z=z_grid_j,
    galaxy_bias_model='per_bin', nuisance_params=pos_nuisance,
  )

  ap_she_she = AngularTwoPoint(shear_tracer, shear_tracer)
  ap_pos_she = AngularTwoPoint(pos_tracer, shear_tracer)
  ap_pos_pos = AngularTwoPoint(pos_tracer, pos_tracer)

  acf_shear = AngularCorrelationFunctionWigner(ap_she_she, ell_j, ks_j)
  acf_ggl   = AngularCorrelationFunctionWigner(ap_pos_she, ell_j, ks_j)
  acf_wt    = AngularCorrelationFunctionWigner(ap_pos_pos, ell_j, ks_j)

  xi_shear = acf_shear.get_xi(theta_rad)
  xi_ggl   = acf_ggl.get_xi(theta_rad)
  xi_wt    = acf_wt.get_xi(theta_rad)

  # Block-concretise (force JAX dispatch to complete before timer stops)
  blocks = []
  for i in range(1, source_ntomo + 1):
    for j in range(i, source_ntomo + 1):
      blocks.append(np.asarray(xi_shear[('SHE', 'SHE', i, j)].array[0, 0, :]))
  for i in range(1, source_ntomo + 1):
    for j in range(i, source_ntomo + 1):
      blocks.append(np.asarray(xi_shear[('SHE', 'SHE', i, j)].array[1, 1, :]))
  for i in range(1, lens_ntomo + 1):
    for j in range(1, source_ntomo + 1):
      if (i, j) in ggl_excl_set:
        continue
      blocks.append(np.asarray(xi_ggl[('POS', 'SHE', i, j)].array[0, :]))
  for i in range(1, lens_ntomo + 1):
    blocks.append(np.asarray(xi_wt[('POS', 'POS', i, i)].array))
  return np.concatenate(blocks)


# ── Generate perturbed cosmologies ────────────────────────────
rng = np.random.default_rng(rng_seed)

def perturb_cosmo(params):
  out = dict(params)
  for k, scale in cosmo_perturb_scales.items():
    out[k] = params[k] * (1.0 + scale * rng.standard_normal())
  return out

def perturb_dict(d, scale=0.05):
  return {k: v * (1.0 + scale * rng.standard_normal()) for k, v in d.items()}

cosmo_list  = [perturb_cosmo(fid_params)  for _ in range(n_warmup + n_bench)]
nla_list    = [perturb_dict(fid_nla)      for _ in range(n_warmup + n_bench)]
bias_list   = [perturb_dict(fid_bias)     for _ in range(n_warmup + n_bench)]
dz_lens_l   = [rng.normal(0.0, dz_sigma, lens_ntomo)
               for _ in range(n_warmup + n_bench)]
dz_source_l = [rng.normal(0.0, dz_sigma, source_ntomo)
               for _ in range(n_warmup + n_bench)]

# ── Warmup ─────────────────────────────────────────────────────
print(f"JAX backend: {jax.default_backend()} (CPU pinned)")
print(f"Setup: {lens_ntomo} lens × {source_ntomo} source bins, "
      f"{n_theta} θ in [{theta_min_arcmin}, {theta_max_arcmin}] arcmin, "
      f"ell ∈ [{ell.min()}, {ell.max()}]")
print(f"Warmup ({n_warmup} run(s), not timed)...")
for k in range(n_warmup):
  pert = make_perturbations(cosmo_list[k])
  _ = compute_3x2pt_cloelib(pert, nla_list[k], bias_list[k],
                            dz_lens_l[k], dz_source_l[k])
  print(f"  warmup {k+1}/{n_warmup} done")

# ── Benchmark — CAMB time tracked separately ──────────────────
print(f"\nBenchmark ({n_bench} cosmologies, cloelib-only timing)...")
camb_times = []
cloelib_times = []
for k in range(n_bench):
  params = cosmo_list[n_warmup + k]
  nla    = nla_list[n_warmup + k]
  bias   = bias_list[n_warmup + k]
  dz_l   = dz_lens_l[n_warmup + k]
  dz_s   = dz_source_l[n_warmup + k]

  t0 = time.perf_counter()
  pert = make_perturbations(params)
  t_camb = time.perf_counter() - t0

  t0 = time.perf_counter()
  dv = compute_3x2pt_cloelib(pert, nla, bias, dz_l, dz_s)
  t_cl = time.perf_counter() - t0

  camb_times.append(t_camb)
  cloelib_times.append(t_cl)
  print(f"  [{k+1:2d}/{n_bench}]  cloelib={t_cl:6.2f}s  (CAMB={t_camb:5.2f}s)  "
        f"Om_c={params['Omega_cdm0']:.4f}, As={params['As']:.3e}, "
        f"AIA={nla['AIA']:.3f}, b0={bias['b0']:.3f}")

camb_times = np.array(camb_times)
cloelib_times = np.array(cloelib_times)

print(f"\n  cloelib pipeline (TIMED):")
print(f"    mean:   {cloelib_times.mean():.2f}s")
print(f"    std:    {cloelib_times.std():.2f}s")
print(f"    min:    {cloelib_times.min():.2f}s")
print(f"    max:    {cloelib_times.max():.2f}s")
print(f"    total:  {cloelib_times.sum():.2f}s for {n_bench} evals")
print(f"\n  CAMB setup (NOT in cloelib total, for reference):")
print(f"    mean:   {camb_times.mean():.2f}s")
print(f"    total:  {camb_times.sum():.2f}s")
print(f"\n  datavec length: {len(dv)}")
