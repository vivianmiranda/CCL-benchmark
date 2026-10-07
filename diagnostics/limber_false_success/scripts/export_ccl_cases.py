"""Export shared CCL inputs for a CoCoA quadrature-only diagnostic.

This is not a replacement for either code's normal cosmology pipeline.
Each code will be tested against its own finer quadrature, with fixed inputs.
"""
import json
from pathlib import Path
import sys
import numpy as np
import pyccl as ccl

from paths import output_dir, INPUTS, VENDOR
root = output_dir()
sys.path.insert(0, str(VENDOR))
from validate_limber_bad_point import make_cosmology, make_tracers

frozen = np.load(INPUTS / 'limber_outlier_tomography.npz')
z = frozen['z']
nz = [frozen[f'nz_{i}'] for i in range(5)]
np.savetxt(root / 'frozen_nz.txt', np.column_stack([z, *nz]))
w_root = json.loads((root / 'root_spike.json').read_text())['w']
cases = []
h = .6727
k = np.geomspace(1.e-5, 1.e3, 1500)
z_power = np.linspace(0., 4.4, 221)
z_distance = np.linspace(0., 50.1, 10021)
z_growth = np.linspace(0., 4.4, 881)
ell = np.unique(np.r_[np.geomspace(20., 2000., 25), 355.655882])
for index, w in enumerate([w_root-5.e-5, w_root, w_root+5.e-5]):
    cosmo = make_cosmology(w)
    tracers = make_tracers(cosmo, z, nz, frozen['bin_centers'])
    bias = 1.05/ccl.growth_factor(cosmo, 1./(1.+frozen['bin_centers']))
    linear = np.array([ccl.linear_matter_power(cosmo, k, 1./(1.+zz))
                       for zz in z_power])*h**3
    nonlinear = np.array([ccl.nonlin_matter_power(cosmo, k, 1./(1.+zz))
                          for zz in z_power])*h**3
    growth = ccl.growth_factor(cosmo, 1./(1.+z_growth))*(1.+z_growth)
    growth /= growth[-1]
    tables = dict(log10k_2D=np.log10(k/h), z_2D=z_power,
                  lnP_linear=np.log(linear).ravel(order='F'),
                  lnP_nonlinear=np.log(nonlinear).ravel(order='F'),
                  lnP_linear_cb=np.log(linear).ravel(order='F'),
                  G=growth, z_G=z_growth, z_1D=z_distance,
                  chi=ccl.comoving_radial_distance(cosmo, 1./(1.+z_distance))*h,
                  omegan2=np.array(0.))
    np.savez(root/f'shared_input_{index}.npz', **tables, bias=bias, ell=ell,
             w=np.array(w))
    if index == 1:
        np.savez_compressed(root/'plot_geometry.npz',
                            z_1D=tables['z_1D'], chi=tables['chi'])
    values = {method:float(ccl.angular_cl(cosmo, tracers[3], tracers[4],
              355.655882, l_limber=-1, limber_integration_method=method))
              for method in ('qag_quad', 'spline')}
    cases.append(dict(index=index, w=w, values=values))
    print(cases[-1], flush=True)
(root/'shared_cases.json').write_text(json.dumps(cases, indent=2)+'\n')
