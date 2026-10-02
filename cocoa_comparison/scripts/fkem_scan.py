"""CCL FKEM robustness scan on one lens-source pair: the non-Limber/Limber
ratio of C_gs for several fkem_chi_min and fkem_Nchi (CCL's documented
knobs), on CoCoA's inputs (same CosmologyCalculator as ccl_compute.py).
python fkem_scan.py <cocoa export .npz> <lens> <source>"""
import json, os, sys
import numpy as np
if not hasattr(np, "trapz"): np.trapz = np.trapezoid
sys.argv = [sys.argv[0], sys.argv[1], os.environ.get("VARIANT", "matched"), "/dev/null"] + sys.argv[2:]
lens, srcb = int(sys.argv[4]), int(sys.argv[5])
H = os.path.dirname(os.path.abspath(__file__))
code = open(os.path.join(H, "ccl_compute.py")).read()
code = code[:code.index("# ---- harmonic mode")]          # the cosmology and tracer setup only
exec(code)
ells = np.array([10., 20., 35., 50., 75., 100., 149.])
lim = ccl.angular_cl(cosmo, lens_t[lens], src_t[srcb], ells)
print("pair (%d,%d): C_gs non-Limber/Limber - 1 at l =" % (lens, srcb), ells.astype(int))
for chimin in (None, 1.0, 10.0, 50.0, 200.0):
  for nchi in (None, 2000, 8000):
    kw = dict(l_limber=1000, non_limber_integration_method="FKEM")
    if chimin is not None: kw["fkem_chi_min"] = chimin
    if nchi is not None: kw["fkem_Nchi"] = nchi
    try:
      nl = ccl.angular_cl(cosmo, lens_t[lens], src_t[srcb], ells, **kw)
      print("  chi_min %-6s Nchi %-5s" % (chimin, nchi), np.round(nl/lim - 1, 4), flush=True)
    except Exception as e:
      print("  chi_min %-6s Nchi %-5s FAIL %s" % (chimin, nchi, str(e)[:60]), flush=True)
