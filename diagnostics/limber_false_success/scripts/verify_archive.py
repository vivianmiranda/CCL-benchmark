"""Audit saved results and README numbers; no CCL/CoCoA calculation is run."""
import hashlib
import json
import numpy as np

from collect_cocoa import collect
from paths import PACKAGE, RESULTS, INPUTS


def main():
    manifest = json.loads((PACKAGE/'archive_manifest.json').read_text())
    for name, digest in manifest['files_sha256'].items():
        assert hashlib.sha256((PACKAGE/name).read_bytes()).hexdigest() == digest, name
    summary = json.loads((RESULTS/'study_summary.json').read_text())
    assert hashlib.sha256((INPUTS/'limber_outlier_tomography.npz').read_bytes()).hexdigest() == summary['provenance']['tomography_sha256']
    record = json.loads((RESULTS/'root_spike.json').read_text())
    assert record['qag'][0]['status'] == [0]
    assert record['qag'][0]['function_calls'] == 41
    np.testing.assert_allclose(record['qag'][0]['value'], record['native']['qag_quad'], rtol=1e-14)
    error = record['native']['qag_quad']/record['cquad']['value']-1
    np.testing.assert_allclose(error, summary['reproduced_case']['native_qag_fractional_error_vs_cquad'], rtol=1e-9)
    measured = collect(RESULTS)
    saved = json.loads((RESULTS/'cocoa_convergence.json').read_text())
    for actual, expected in zip(measured, saved):
        assert actual['case'] == expected['case']
        for a, e in zip(actual['levels'], expected['levels']):
            for key in a:
                np.testing.assert_allclose(a[key], e[key], rtol=1e-9, atol=3e-14,
                                           err_msg=f'case {actual["case"]}: {key}')
    max_cov = max(r['levels'][0]['covariance_gg_max_abs_fractional_change_vs_level4'] for r in measured)
    max_dv = max(r['levels'][0]['data_vector_autos_max_abs_fractional_change_vs_level4'] for r in measured)
    max_refine = max(r['levels'][3]['covariance_gg_max_abs_fractional_change_vs_level4'] for r in measured)
    for value, key in [(max_cov, 'covariance_default_max_abs_fractional_change'),
                       (max_dv, 'data_vector_autos_default_128_vs_1024_max_abs_fractional_change'),
                       (max_refine, 'covariance_512_vs_1024_max_abs_fractional_change')]:
        np.testing.assert_allclose(value, summary['cocoa'][key], rtol=1e-9, atol=3e-14)
    nodes = np.load(RESULTS/'cocoa_case1_level0_nodes.npz')
    original = np.load(RESULTS/'cocoa_case1_level0.npz')
    for key in original.files:
        np.testing.assert_array_equal(nodes[key], original[key])
    canonical, _ = np.polynomial.legendre.leggauss(96)
    a = nodes['geometry'][0]
    for left, right in zip(nodes['a_edges'][:-1], nodes['a_edges'][1:]):
        selected = a[(a > left) & (a < right)]
        assert len(selected) == 96
        np.testing.assert_allclose(np.sort(selected), left+(canonical+1)*(right-left)/2,
                                   rtol=0, atol=3e-16)
    print(json.dumps({'verified_archive_files': len(manifest['files_sha256']),
        'native_QAG_error_percent': 100*error,
        'CoCoA_covariance_spectra_default_max_change_percent': 100*max_cov,
        'CoCoA_data_vector_autos_default_max_change_percent': 100*max_dv,
        'CoCoA_covariance_512_vs_1024_max_fractional_change': max_refine,
        'actual_CoCoA_geometry_nodes_verified': len(a)}, indent=2))


if __name__ == '__main__':
    main()
