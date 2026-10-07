"""Re-evaluate the saved failing CCL integral without modifying CCL or GSL.

The last digits of the accidental crossing depend on the installed build.
Use --w to try a root from find_first_panel.py on a different installation.
"""
import argparse
import ctypes as ct
import hashlib
import json
import platform
import sys

import numpy as np
import pyccl as ccl

from probe_quadrature import probe, root, gsl, gsl_path, INPUTS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--w', type=float, default=-1.0389272675840322)
    args = parser.parse_args()
    output = root / 'root_spike.json'
    if output.exists():
        raise FileExistsError(f'{output} exists; choose a fresh QAG_DIAGNOSTIC_OUTPUT')
    record, logk, integrand = probe(args.w, [3, 4], 355.655882)
    output.write_text(json.dumps(record, indent=2) + '\n')
    np.savez(root/'root_spike_integrand.npz', logk=logk, integrand=integrand)
    provenance = {
        'ccl_version': ccl.__version__, 'ccl_file': ccl.__file__,
        'gsl_version': ct.c_char_p.in_dll(gsl, 'gsl_version').value.decode(),
        'gsl_library': gsl_path, 'python': sys.version, 'platform': platform.platform(),
        'pr_head': '8978d57e02e07e6557428ff50e134f3f21d70ab9',
        'tomography_sha256': hashlib.sha256(
            (INPUTS/'limber_outlier_tomography.npz').read_bytes()).hexdigest(),
        'settings': {name: getattr(ccl.gsl_params, name) for name in (
            'INTEGRATION_LIMBER_EPSREL', 'INTEGRATION_LIMBER_GAUSS_KRONROD_POINTS',
            'N_ITERATION')},
    }
    (root/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    error = record['native']['qag_quad'] / record['cquad']['value'] - 1
    print(json.dumps({'w': args.w, 'native_QAG_over_CQUAD_minus_one': error,
                      'QAG_status': record['qag'][0]['status'],
                      'QAG_calls': record['qag'][0]['function_calls'],
                      'output': str(output)}, indent=2))


if __name__ == '__main__':
    main()
