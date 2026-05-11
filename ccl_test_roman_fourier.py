import time
import numpy as np
import pyccl as ccl
from pyccl import nl_pt as pt

# ── Configuration ──────────────────────────────────────────────
nz_lens_file = "roman_example1.nz"
nz_source_file = "roman_example1.nz"
lens_ntomo = 8
source_ntomo = 8

# Fourier-space binning
n_cl = 15
l_min = 30
l_max = 4000          # for clustering & GGL
l_max_shear = 4000    # for cosmic shear

# GGL pairs to exclude (format: [lens_bin, source_bin])
ggl_exclude = [[6, 0], [7, 0], [7, 1]]

# Benchmark settings
n_warmup = 1
n_bench = 10
rng_seed = 42

# Fiducial cosmology
fid_params = dict(
  Omega_c=0.27, 
  Omega_b=0.045, 
  h=0.67,
  sigma8=0.83, 
  n_s=0.96,
)

# Fiducial TATT parameters (perturbed in benchmark loop)
fid_tatt = dict(A1=0.7, A2=-1.36, eta1=-1.7, eta2=-2.5)
z_pivot = 0.62

# Galaxy bias (fiducial — perturbed in benchmark loop)
fid_bias = {f'b{i}': 1.2 + 0.1 * i for i in range(lens_ntomo)}

# Photo-z shifts: fiducial = 0, perturbed additively (multiplicative
# perturbation of 0 doesn't do anything). sigma ~ 0.01 is realistic
# for stage-IV surveys.
dz_sigma = 0.01

# z grid for IA
z_ia = np.linspace(0.0, 3.0, 200)

# ── Read n(z) once (doesn't depend on cosmology) ──────────────
def load_nz(filename, ntomo):
  data = np.loadtxt(filename)
  z = data[:, 0]
  nz = [data[:, i + 1] for i in range(ntomo)]
  return z, nz

z_lens, nz_lens = load_nz(nz_lens_file, lens_ntomo)
z_source, nz_source = load_nz(nz_source_file, source_ntomo)

# ── Ell bins (log-spaced) ─────────────────────────────────────
ell_shear = np.geomspace(l_min, l_max_shear, n_cl)
ell_other = np.geomspace(l_min, l_max, n_cl)
print(f"ell_shear: {ell_shear[0]:.1f} to {ell_shear[-1]:.1f}, {n_cl} pts")
print(f"ell_other: {ell_other[0]:.1f} to {ell_other[-1]:.1f}, {n_cl} pts")

ggl_set = {tuple(p) for p in ggl_exclude}

# ── Single 3x2pt evaluation (Fourier space) ───────────────────
def compute_3x2pt(params, tatt, bias, dz_lens, dz_source):
  """Compute Fourier-space 3x2pt data vector."""

  cosmo = ccl.Cosmology(
    **params,
    transfer_function='eisenstein_hu',
  )

  galaxy_bias = [bias[f'b{i}'] for i in range(lens_ntomo)]

  # TATT: A → c
  A1_z = tatt['A1'] * ((1.0 + z_ia) / (1.0 + z_pivot))**tatt['eta1']
  A2_z = tatt['A2'] * ((1.0 + z_ia) / (1.0 + z_pivot))**tatt['eta2']
  c1, c2, cdelta = pt.translate_IA_norm(
    cosmo, z=z_ia, a1=A1_z, a1delta=A1_z, a2=A2_z
  )

  # PT tracers + calculator
  ptt_m = pt.PTMatterTracer()
  ptt_ia = pt.PTIntrinsicAlignmentTracer(
    c1=(z_ia, c1), c2=(z_ia, c2), cdelta=(z_ia, cdelta),
  )
  a_arr = 1.0 / (1.0 + np.linspace(0, 3.0, 50)[::-1])
  ptc = pt.EulerianPTCalculator(
    with_NC=True, with_IA=True,
    log10k_min=-4, log10k_max=2, nk_per_decade=20,
    cosmo=cosmo, a_arr=a_arr,
  )
  ptc.update_ingredients(cosmo)

  pk_mm = ptc.get_biased_pk2d(ptt_m, tracer2=ptt_m)
  pk_mi = ptc.get_biased_pk2d(ptt_m, tracer2=ptt_ia)
  pk_ii = ptc.get_biased_pk2d(ptt_ia, tracer2=ptt_ia)

  # Angular tracers (with photo-z shifts: z -> z - dz, valid z >= 0)
  source_L, source_IA = [], []
  for i in range(source_ntomo):
    z_shifted = z_source - dz_source[i]
    mask = z_shifted >= 0.0
    dndz_src = (z_shifted[mask], nz_source[i][mask])
    source_L.append(ccl.WeakLensingTracer(
      cosmo, dndz=dndz_src,
      has_shear=True, ia_bias=None,
    ))
    source_IA.append(ccl.WeakLensingTracer(
      cosmo, dndz=dndz_src,
      has_shear=False,
      ia_bias=(z_shifted[mask], np.ones(mask.sum())),
      use_A_ia=True,
    ))
  lens_tracers = []
  for i in range(lens_ntomo):
    z_shifted = z_lens - dz_lens[i]
    mask = z_shifted >= 0.0
    lens_tracers.append(ccl.NumberCountsTracer(
      cosmo, has_rsd=False,
      dndz=(z_shifted[mask], nz_lens[i][mask]),
      bias=(z_shifted[mask], galaxy_bias[i] * np.ones(mask.sum())),
    ))

  # Shear C_ell with TATT (l_min >= 30, Limber is accurate)
  cls_shear = {}
  for i in range(source_ntomo):
    for j in range(i, source_ntomo):
      cl_gg = ccl.angular_cl(cosmo, source_L[i], source_L[j],
                             ell_shear, p_of_k_a=pk_mm)
      cl_gi = ccl.angular_cl(cosmo, source_L[i], source_IA[j],
                             ell_shear, p_of_k_a=pk_mi)
      cl_ig = ccl.angular_cl(cosmo, source_IA[i], source_L[j],
                             ell_shear, p_of_k_a=pk_mi)
      cl_ii = ccl.angular_cl(cosmo, source_IA[i], source_IA[j],
                             ell_shear, p_of_k_a=pk_ii)
      cls_shear[(i, j)] = cl_gg + cl_gi + cl_ig + cl_ii

  # GGL C_ell with TATT
  cls_ggl = {}
  for i in range(lens_ntomo):
    for j in range(source_ntomo):
      cl_gG = ccl.angular_cl(cosmo, lens_tracers[i], source_L[j],
                             ell_other, p_of_k_a=pk_mm)
      cl_gI = ccl.angular_cl(cosmo, lens_tracers[i], source_IA[j],
                             ell_other, p_of_k_a=pk_mi)
      cls_ggl[(i, j)] = cl_gG + cl_gI

  # Clustering C_ell, autos only (Roman convention)
  cls_clustering = {}
  for i in range(lens_ntomo):
    cls_clustering[(i, i)] = ccl.angular_cl(
      cosmo, lens_tracers[i], lens_tracers[i], ell_other,
      p_of_k_a=pk_mm,
    )

  # Flatten data vector: C_ell directly (no real-space transform)
  datavec = np.concatenate([
    *[cls_shear[(i, j)] for i in range(source_ntomo)
      for j in range(i, source_ntomo)],
    *[cls_ggl[(i, j)] for i in range(lens_ntomo)
      for j in range(source_ntomo) if (i, j) not in ggl_set],
    *[cls_clustering[(i, i)] for i in range(lens_ntomo)],
  ])
  return datavec

# ── Generate perturbed cosmologies + nuisance ────────────────
rng = np.random.default_rng(rng_seed)

def perturb(params, scale=0.05):
  """Multiply each param by (1 + N(0, scale))."""
  return {k: v * (1.0 + scale * rng.standard_normal())
          for k, v in params.items()}

cosmo_list = [perturb(fid_params) for _ in range(n_warmup + n_bench)]
tatt_list = [perturb(fid_tatt) for _ in range(n_warmup + n_bench)]
bias_list = [perturb(fid_bias) for _ in range(n_warmup + n_bench)]
dz_lens_list = [rng.normal(0.0, dz_sigma, lens_ntomo)
                for _ in range(n_warmup + n_bench)]
dz_source_list = [rng.normal(0.0, dz_sigma, source_ntomo)
                  for _ in range(n_warmup + n_bench)]

# ── Warmup ─────────────────────────────────────────────────────
print(f"\nWarmup ({n_warmup} run(s), not timed)...")
for k in range(n_warmup):
  _ = compute_3x2pt(cosmo_list[k], tatt_list[k], bias_list[k],
                    dz_lens_list[k], dz_source_list[k])
  print(f"  warmup {k+1}/{n_warmup} done")

# ── Benchmark ──────────────────────────────────────────────────
print(f"\nBenchmark ({n_bench} cosmologies)...")
times = []
for k in range(n_bench):
  params = cosmo_list[n_warmup + k]
  tatt = tatt_list[n_warmup + k]
  bias = bias_list[n_warmup + k]
  dz_l = dz_lens_list[n_warmup + k]
  dz_s = dz_source_list[n_warmup + k]
  t0 = time.perf_counter()
  dv = compute_3x2pt(params, tatt, bias, dz_l, dz_s)
  dt = time.perf_counter() - t0
  times.append(dt)
  print(f"  [{k+1:2d}/{n_bench}]  {dt:6.2f}s  "
        f"(Om_c={params['Omega_c']:.4f}, "
        f"sig8={params['sigma8']:.4f}, "
        f"A1={tatt['A1']:.3f}, "
        f"b0={bias['b0']:.3f}, "
        f"dz_s0={dz_s[0]:+.4f})")

times = np.array(times)
print(f"\n  mean:   {times.mean():.2f}s")
print(f"  std:    {times.std():.2f}s")
print(f"  min:    {times.min():.2f}s")
print(f"  max:    {times.max():.2f}s")
print(f"  total:  {times.sum():.2f}s for {n_bench} evals")

# Data vector summary
n_shear = source_ntomo * (source_ntomo + 1) // 2
n_ggl = lens_ntomo * source_ntomo - len(ggl_exclude)
n_wth = lens_ntomo

print(f"\n  shear C_ell:      {n_shear} pairs x {n_cl} = {n_shear * n_cl}")
print(f"  GGL C_ell:        {n_ggl} pairs x {n_cl} = {n_ggl * n_cl}"
      f"  ({len(ggl_exclude)} excluded)")
print(f"  clustering C_ell: {n_wth} pairs x {n_cl} = {n_wth * n_cl}  (auto only)")
print(f"  datavec length:   {len(dv)}")