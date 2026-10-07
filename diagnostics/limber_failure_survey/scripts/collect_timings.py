"""Verify and archive the native CCL QAG/CQUAD timing measurements."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from common import digest, write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--destination',type=Path,required=True)
    a=p.parse_args()
    if a.destination.exists():raise FileExistsError(a.destination)
    folders=['baseline_1thread_verified','paired_1thread_verified',
             'baseline_6threads_verified','paired_6threads_verified',
             'baseline_narrow_verified','paired_narrow_verified']
    files=[]; rows=[]; records={}
    for folder in folders:
        root=a.source/folder
        record=json.loads((root/'timing.json').read_text());records[folder]=record
        assert record['status']=='completed'
        files.append(root/'timing.json')
        for case in record['cases']:
            assert case['verified_omp_threads_for_timing']==record['threads']
            path=root/f'{case["case"]["id"]}.npz'
            assert digest(path)==case['arrays_sha256']
            files.append(path)
            if record['baseline_only']:continue
            for name,w in case['workloads'].items():
                assert w['baseline_qag_bitwise_equal']
                for method,samples in w['samples'].items():
                    assert len(samples)==9
                    for s in samples:
                        assert not s['counters']['workspace_allocations']
                        assert not s['counters']['unexpected_callbacks']
                        if method=='cquad_diagnostic':
                            assert not s['counters']['cquad_failures'] and not s['warnings']
                rows.append({'case':case['case'],'threads':record['threads'],'workload':name,
                    'QAG_ms':1000*w['summary']['qag_quad']['public_seconds_per_workload']['mean'],
                    'QAG_scatter_ms':1000*w['summary']['qag_quad']['public_seconds_per_workload']['sample_scatter'],
                    'CQUAD_ms':1000*w['summary']['cquad_diagnostic']['public_seconds_per_workload']['mean'],
                    'CQUAD_scatter_ms':1000*w['summary']['cquad_diagnostic']['public_seconds_per_workload']['sample_scatter'],
                    'CQUAD_over_QAG':w['cquad_over_qag_ratio']['public_seconds_per_workload'],
                    'error_percent':{k:100*v for k,v in w['max_fractional_error_vs_strict_reference'].items()},
                    'reference_scope':w['reference_scope'], 'unmodified_control_ratio':w['unmodified_qag_control']['instrumented_over_unmodified_mean_ratio']})
    for ident in ['u000__pr_original','known__pr_original']:
        one=np.load(a.source/'paired_1thread_verified'/f'{ident}.npz')
        six=np.load(a.source/'paired_6threads_verified'/f'{ident}.npz')
        for key in six.files:assert np.array_equal(one[key],six[key]),(ident,key)
    a.destination.mkdir(parents=True)
    for src in files:
        dst=a.destination/src.relative_to(a.source);dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src,dst);assert digest(src)==digest(dst)
    write_json(a.destination/'summary.json',{'rows':rows,'one_vs_six_thread_spectra_bitwise_equal':True,
        'scope':'Warmed native CCL angular_cl; identical requested relative tolerance 1e-4. Setup/CAMB excluded. Isolated CQUAD substitution with reused per-thread workspace; no Python integrand callbacks.',
        'machine':'Apple M2 Pro, macOS 13.7.5. Numerical/compiler jobs sequential; ordinary desktop background remained active.',
        'limits':'Nine interleaved batches per method. Not full data-vector, covariance or MCMC timings. Cost at the same requested tolerance, not matched achieved accuracy when QAG fails. Additional grid reference errors use an internal strict-CQUAD control.',
        'files':{str(f.relative_to(a.destination)):digest(f) for f in a.destination.rglob('*') if f.is_file()},
        'discarded_harness_attempt':'Original point timing aborted because a thread-local workspace was allocated during a measured batch; preserved locally. Accepted runs set and verify the native OpenMP thread count after setup.'})
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
