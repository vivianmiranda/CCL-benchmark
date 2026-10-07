"""Scan the pinned PR tomography with the unmodified installed CCL API."""
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import numpy as np
import pyccl as ccl

from paths import output_dir, INPUTS, VENDOR
root = output_dir()
signal.alarm(600)
frozen = INPUTS / 'limber_outlier_tomography.npz'
data = np.load(frozen)
z = data['z']
centers = data['bin_centers']
pairs = np.array([(i, j) for i in range(5) for j in range(i, 5)])
ell = np.unique(np.concatenate((np.geomspace(20., 2000., 25), [355.655882])))
w_values = np.unique(np.concatenate((np.linspace(-1.08, -.94, 29), [-1.03895, -1.0389, -1.03885])))
values = np.empty((len(w_values), 2, len(pairs), len(ell)))
for iw, w in enumerate(w_values):
    cosmo = ccl.Cosmology(Omega_c=.2664, Omega_b=.0492, h=.6727,
                         n_s=.9645, sigma8=.831, w0=float(w), wa=0.,
                         m_nu=0., mass_split='equal')
    bias = 1.05 / ccl.growth_factor(cosmo, 1. / (1. + centers))
    tracers = []
    for i in range(5):
        tracer = ccl.NumberCountsTracer(cosmo, has_rsd=False,
            dndz=(z, data[f'nz_{i}']), bias=(z, np.full_like(z, bias[i])))
        tracers.append(tracer)
    for ip, (i, j) in enumerate(pairs):
        for im, method in enumerate(('qag_quad', 'spline')):
            values[iw, im, ip] = ccl.angular_cl(cosmo, tracers[i], tracers[j],
                ell, l_limber=-1, limber_integration_method=method)
    residual = values[iw, 0] / values[iw, 1] - 1
    print(f'w={w:.8f}, max QAG/spline-1 = {np.max(np.abs(residual)):.8g}', flush=True)
np.savez(root / 'ccl_scan.npz', w=w_values, ell=ell, pairs=pairs, values=values)
residual = values[:, 0] / values[:, 1] - 1
where = np.unravel_index(np.argmax(np.abs(residual)), residual.shape)
iw, ip, il = where
report = {
    'ccl_version': ccl.__version__, 'ccl_file': ccl.__file__,
    'tomography_sha256': hashlib.sha256(frozen.read_bytes()).hexdigest(),
    'pr_sha': '8978d57e02e07e6557428ff50e134f3f21d70ab9',
    'counts': {'w': len(w_values), 'pairs': len(pairs), 'ell': len(ell)},
    'gsl': {name:getattr(ccl.gsl_params,name) for name in ('INTEGRATION_LIMBER_EPSREL','INTEGRATION_LIMBER_GAUSS_KRONROD_POINTS','N_ITERATION')},
    'maximum_relative_method_difference': float(abs(residual[where])),
    'worst_w': float(w_values[iw]), 'worst_pair': pairs[ip].tolist(),
    'worst_ell': float(ell[il]), 'worst_qag': float(values[iw, 0, ip, il]),
    'worst_spline': float(values[iw, 1, ip, il]),
}
(root / 'ccl_scan.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2), flush=True)
