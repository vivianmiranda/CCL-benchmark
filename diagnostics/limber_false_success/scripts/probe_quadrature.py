"""Apply native GSL rules to CCL's actual kernels and power interpolator.

This diagnostic uses the tracer's native chi bounds, not the first positive
sample of n(z), because those intervals differ in the frozen PR example.
No CCL function is replaced and no production setting is changed.
"""
import ctypes as ct
import ctypes.util
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
sys.path.insert(0, str(VENDOR))
import validate_limber_bad_point as diagnostic
signal.alarm(600)
data = np.load(INPUTS / 'limber_outlier_tomography.npz')
z = data['z']
bins = [data[f'nz_{i}'] for i in range(5)]
gsl_path = os.environ.get('GSL_LIBRARY')
if not gsl_path:
    candidates = [Path(sys.base_prefix) / 'lib' / name
                  for name in ('libgsl.dylib', 'libgsl.so')]
    gsl_path = next((str(path) for path in candidates if path.exists()),
                    ctypes.util.find_library('gsl'))
if not gsl_path:
    raise RuntimeError('Set GSL_LIBRARY to the GSL library used by your CCL build')
gsl = ct.CDLL(gsl_path)
callback_type = ct.CFUNCTYPE(ct.c_double, ct.c_double, ct.c_void_p)
class Function(ct.Structure):
    """GSL's scalar callback and caller-owned context pointer."""
    _fields_ = [('function', callback_type), ('params', ct.c_void_p)]

void = ct.c_void_p
size = ct.c_size_t
real = ct.c_double
real_pointer = ct.POINTER(real)
function_pointer = ct.POINTER(Function)
gsl.gsl_set_error_handler_off()
for prefix in ('gsl_integration_workspace', 'gsl_integration_cquad_workspace',
               'gsl_integration_glfixed_table'):
    allocate = getattr(gsl, prefix + '_alloc')
    allocate.argtypes = [size]
    allocate.restype = void
    release = getattr(gsl, prefix + '_free')
    release.argtypes = [void]
gsl.gsl_integration_qag.argtypes = [function_pointer, real, real, real, real,
                                    size, ct.c_int, void, real_pointer, real_pointer]
gsl.gsl_integration_cquad.argtypes = [function_pointer, real, real, real, real,
                                      void, real_pointer, real_pointer, ct.POINTER(size)]
gsl.gsl_integration_glfixed_point.argtypes = [real, real, size, real_pointer,
                                             real_pointer, void]


def probe(w, pair, ell):
    """Check one native CCL integral with independent GSL subdivisions."""
    cosmo = diagnostic.make_cosmology(w)
    tracers = diagnostic.make_tracers(cosmo, z, bins, data['bin_centers'])
    first, second = tracers[pair[0]], tracers[pair[1]]
    power = cosmo.get_nonlin_power()
    near = max(first._trc[0].chi_min, second._trc[0].chi_min)
    far = min(first._trc[0].chi_max, second._trc[0].chi_max)
    lo = float(np.log(max(ccl.spline_params.K_MIN, (ell+.5)/far)))
    hi = float(np.log(min(ccl.spline_params.K_MAX, 2*(ell+.5)/near)))
    errors = []
    calls = []
    def integrand(logk, unused):
        """Evaluate exactly the CCL density-kernel product, including P(k,a)."""
        calls.append(logk)
        try:
            return diagnostic.limber_integrand(logk, ell, cosmo, first, second,
                                                power) / (ell+.5)
        except Exception as error:
            errors.append(str(error))
            return float('nan')
    callback = callback_type(integrand)
    function = Function(callback, None)
    qag = []
    for parts in (1, 2, 4, 8):
        edges = np.linspace(lo, hi, parts+1)
        total = 0.; error_total = 0.; status = []
        calls.clear()
        workspace = gsl.gsl_integration_workspace_alloc(1000)
        for left, right in zip(edges[:-1], edges[1:]):
            result, error = real(), real()
            code = gsl.gsl_integration_qag(ct.byref(function), left, right,
                0., 1.e-4, 1000, 4, workspace, ct.byref(result), ct.byref(error))
            total += result.value; error_total += error.value; status.append(code)
        gsl.gsl_integration_workspace_free(workspace)
        qag.append(dict(parts=parts, value=total, estimated_error=error_total,
                        status=status, function_calls=len(calls)))
    workspace = gsl.gsl_integration_cquad_workspace_alloc(1000)
    result, error, count = real(), real(), size()
    code = gsl.gsl_integration_cquad(ct.byref(function), lo, hi, 0., 1.e-8,
        workspace, ct.byref(result), ct.byref(error), ct.byref(count))
    gsl.gsl_integration_cquad_workspace_free(workspace)
    cquad = dict(value=result.value, estimated_error=error.value,
                 status=code, function_calls=count.value)
    fixed = []
    a_min = float(ccl.scale_factor_of_chi(cosmo, (ell+.5)/np.exp(lo)))
    a_max = float(ccl.scale_factor_of_chi(cosmo, (ell+.5)/np.exp(hi)))
    for nodes in (64, 96, 128, 256, 512, 1024):
        rule = gsl.gsl_integration_glfixed_table_alloc(nodes)
        x, weight = real(), real()
        totals = []
        for coordinate in ('logk', 'a'):
            total = 0.
            for node in range(nodes):
                left, right = (lo, hi) if coordinate == 'logk' else (a_min, a_max)
                gsl.gsl_integration_glfixed_point(left, right, node,
                    ct.byref(x), ct.byref(weight), rule)
                if coordinate == 'logk':
                    value = integrand(x.value, None)
                else:
                    a = x.value
                    chi = float(ccl.comoving_radial_distance(cosmo, a))
                    logk = float(np.log((ell+.5)/chi))
                    jacobian = 2997.92458 / (.6727*a*a*ccl.h_over_h0(cosmo,a)*chi)
                    value = integrand(logk, None)*jacobian
                total += weight.value*value
            totals.append(float(total))
        fixed.append(dict(nodes=nodes, logk=totals[0], a=totals[1]))
        gsl.gsl_integration_glfixed_table_free(rule)
    if errors:
        raise RuntimeError(errors)
    native = {}
    for method in ('qag_quad','spline'):
        native[method] = float(ccl.angular_cl(cosmo, first, second, ell,
            l_limber=-1, limber_integration_method=method))
    plot_logk = np.linspace(lo, hi, 2001)
    plot_y = np.array([integrand(float(x), None) for x in plot_logk])
    record = dict(w=w, pair=pair, ell=ell, bounds=[lo,hi], native=native,
                  qag=qag, cquad=cquad, fixed=fixed)
    return record, plot_logk, plot_y

if __name__ == '__main__':
    cases = [(-1.0389, [3,4],355.655882), (-.94,[1,3],632.4555320336759)]
    records = []
    for number, (w,pair,ell) in enumerate(cases):
        record, logk, integrand = probe(w,pair,ell)
        records.append(record)
        np.savez(root/f'ccl_integrand_{number}.npz',logk=logk,integrand=integrand)
        (root/'gsl_probe.json').write_text(json.dumps(records,indent=2)+'\n')
        print(json.dumps(record,indent=2),flush=True)
