"""CoCoA side of the CoCoA vs DESC-CCL comparison (cocoa env).

python cocoa_export.py <project> <model> <out.npz>
project: lsst_y1 | roman_real (EXAMPLE_EVALUATE2, 3x2pt, NLA)
model:   fid | omm_lo | omm_hi | ns_lo | ns_hi

Photo-z shifts and shear calibration are set to zero (magnification,
b2 and point masses already are), so both codes see the same n(z) and no
calibration. Saves the full (unmasked) theory vector, the masked inverse
covariance and mask, the block sizes, the theta binning, the CAMB linear
and nonlinear P(k) (Mpc^3, k in 1/Mpc) cosmolike was given, and every
input CCL needs."""
import os, sys, json, re
import numpy as np
from cobaya.model import get_model
from cobaya.yaml import yaml_load_file

project, mname, out = sys.argv[1:4]
MODELS = {"fid": {}, "omm_lo": {"omegam": 0.25}, "omm_hi": {"omegam": 0.35},
          "ns_lo": {"ns": 0.92}, "ns_hi": {"ns": 1.01}}
R = os.environ["ROOTDIR"]; os.chdir(R)
info = yaml_load_file("./projects/%s/EXAMPLE_EVALUATE2.yaml" % project)
ov = dict(info["sampler"]["evaluate"]["override"]); info.pop("sampler"); info.pop("output", None)
info["timing"] = False
name = list(info["likelihood"])[0]
lb = info["likelihood"][name]
lb["IA_model"] = 0
# optional likelihood overrides (JSON), e.g. a high-accuracy CoCoA run
lb.update(json.loads(os.environ.get("COCOA_OVR", "{}")))
# mode "dv": a scratch copy of the data folder whose dataset keeps every
# point (mask_file = the project's ones.mask), so the theory vector comes
# out full length; mode "cov": the shipped dataset, for the masked inverse
# covariance (re-initializing the data in one process is not supported)
mode = os.environ.get("EXPORT_MODE", "dv")
if mode == "dv":
  src = os.path.join(R, lb["path"]); dst = os.path.join(os.path.dirname(os.path.abspath(out)), "data_" + project)
  os.makedirs(dst, exist_ok=True)
  for fn in os.listdir(src):
    t = os.path.join(dst, fn)
    if not os.path.lexists(t): os.symlink(os.path.join(src, fn), t)
  ds_name = lb["data_file"]; ds_path = os.path.join(dst, "ones_" + ds_name)
  with open(os.path.join(src, ds_name)) as fin, open(ds_path + ".tmp", "w") as fo:
    for line in fin:
      fo.write("mask_file = ones.mask\n" if line.strip().startswith("mask_file") else line)
  os.replace(ds_path + ".tmp", ds_path)
  lb["path"] = dst; lb["data_file"] = "ones_" + ds_name
model = get_model(info)
lik = model.likelihood[name]
import importlib
ci = importlib.import_module("cosmolike_%s_interface" % project)
point = {p: float(ov[p]) for p in model.parameterization.sampled_params() if p in ov}
for p in model.parameterization.sampled_params():
  if p not in point:
    ref = model.info()["params"][p].get("ref")
    point[p] = float(ref["loc"] if isinstance(ref, dict) else ref)
for p in point:
  if re.search(r"_DZ_|_M[0-9]+$|_PM[0-9]+$|_BMAG_", p): point[p] = 0.0
point.update(MODELS[mname])
post = model.logposterior(point, cached=False)
inputs = model.parameterization.to_input(point)

icov = np.array(ci.get_inv_cov_masked()); mask = np.array(ci.get_mask())
sizes = list(ci.compute_data_vector_3x2pt_real_sizes())
dv = np.array(ci.compute_data_vector_masked())

# harmonic-space spectra at a few multipoles (Limber, as cosmolike's
# *_tomo_limber bindings compute them), for the C_l-level comparison
ell_h = np.geomspace(20.0, 5000.0, 24)
cl_ss = np.array(ci.C_ss_tomo_limber(l=ell_h)); cl_gs = np.array(ci.C_gs_tomo_limber(l=ell_h))
cl_gg = np.array(ci.C_gg_tomo_limber(l=ell_h))
# the P(k) cosmolike was given (CAMB through cobaya, same extrapolation)
prov = model.provider
h = prov.get_param("H0")/100.0
kw = dict(extrap_kmin=1e-6, extrap_kmax=2.5e2*lik.accuracyboost)
PKL = prov.get_Pk_interpolator(("delta_tot", "delta_tot"), nonlinear=False, **kw)
PKN = prov.get_Pk_interpolator(("delta_tot", "delta_tot"), nonlinear=True, **kw)
zg = np.linspace(0.0, 6.0, 301); kg = np.logspace(-5, np.log10(200.0), 700)
lnPL = PKL.logP(zg, kg); lnPN = PKN.logP(zg, kg)
z1 = np.array(lik.z_interp_1D); z1 = z1[z1 <= 50.0]
chi1 = np.array(prov.get_comoving_radial_distance(z1))
chi = np.interp(zg, z1, chi1)
D0 = np.sqrt(PKL.P(zg, 5e-4)/PKL.P(0, 5e-4))
# the growth table over CAMB's full z range: CCL's RSD kernel asks for
# D(a) beyond the P(k) table (z = 6) for the high-z lens bins
zD = np.concatenate((np.linspace(0.0, 6.0, 301), np.linspace(6.05, 49.0, 400)))
DD = np.sqrt(PKL.P(zD, 5e-4)/PKL.P(0, 5e-4))

ds = {}
for line in open(os.path.join(lik.path, lik.data_file)):
  if "=" in line and not line.strip().startswith("#"):
    k, v = line.split("=", 1); ds[k.strip()] = v.strip()
pars = {p: float(inputs[p]) for p in lik.input_params if p in inputs}
for p in ("H0", "omegam", "omegab", "ns", "As", "mnu", "w", "wa", "omnuh2"):
  try: pars[p] = float(prov.get_param(p))
  except Exception: pass
np.savez(out, dv=dv, icov=icov, mask=mask, sizes=sizes, chi2=-2.0*post.loglikes[0],
         ell_h=ell_h, cl_ss=cl_ss, cl_gs=cl_gs, cl_gg=cl_gg, zg=zg, kg=kg, lnPL=lnPL, lnPN=lnPN, chi=chi, D0=D0, zD=zD, DD=DD, z1=z1, chi1=chi1,
         meta=json.dumps(dict(project=project, model=mname, point=point, pars=pars,
                              dataset=ds, path=lik.path, rootdir=R, ntheta=int(lik.ntheta),
                              theta_min=float(lik.theta_min_arcmin), theta_max=float(lik.theta_max_arcmin),
                              lens_ntomo=int(lik.lens_ntomo), source_ntomo=int(lik.source_ntomo),
                              accuracyboost=float(lik.accuracyboost),
                              adopt_limber_gs=getattr(lik, "adopt_limber_gs", None),
                              adopt_limber_gg=getattr(lik, "adopt_limber_gg", None))))
print("EXPORT %s %s mode %s chi2 %.4f dv %d kept %d sizes %s" % (project, mname, mode, -2*post.loglikes[0], len(dv), int(mask.sum()), sizes), flush=True)
