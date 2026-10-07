"""Native-C QAG/CQUAD diagnostic, with a separate no-preload baseline.

This does not add CQUAD to the installed CCL API. A process-local macOS
shim changes only the GSL call during an explicitly armed angular_cl call.
The integrand remains the installed CCL C callback. No Python callbacks
occur inside either quadrature. Run only on a quiet machine.

First run --baseline-only without DYLD_INSERT_LIBRARIES. Then run this
script in a new process with --shim and --baseline pointing to that folder.
Both stages use fresh --output folders. Parent controls compilation/runs.
"""
import argparse
import ctypes as ct
import json
import os
from pathlib import Path
import platform
import statistics
import sys
import time
import traceback
import warnings

COUNTERS = ('qag_calls', 'cquad_calls', 'unexpected_callbacks',
            'workspace_allocations', 'kernel_nanoseconds', 'qag_failures',
            'cquad_failures', 'fallback_cquad_calls')


class Shim:
    def __init__(self, path):
        self.path = path.resolve()
        inserted = [Path(x).resolve() for x in os.environ.get('DYLD_INSERT_LIBRARIES', '').split(':') if x]
        if self.path not in inserted:
            raise RuntimeError('The shim must be preloaded before Python starts')
        self.lib = ct.CDLL(str(self.path))
        self.lib.limber_timing_ready.restype = ct.c_int
        self.lib.limber_timing_begin.argtypes = [ct.c_int, ct.c_double]
        self.lib.limber_timing_begin.restype = ct.c_int
        self.lib.limber_timing_end.restype = None
        self.lib.limber_timing_counter.argtypes = [ct.c_int]
        self.lib.limber_timing_counter.restype = ct.c_uint64
        self.lib.limber_timing_status_count.argtypes = [ct.c_int, ct.c_int]
        self.lib.limber_timing_status_count.restype = ct.c_uint64
        if not self.lib.limber_timing_ready():
            raise RuntimeError('Cannot resolve original native GSL functions safely')

    def begin(self, mode, epsrel=0.):
        if self.lib.limber_timing_begin(mode, epsrel):
            raise RuntimeError('Cannot arm the isolated native timing shim')

    def end(self):
        self.lib.limber_timing_end()

    def counters(self):
        result = {key: int(self.lib.limber_timing_counter(i)) for i, key in enumerate(COUNTERS)}
        result['status_counts'] = {
            name: {str(status): int(self.lib.limber_timing_status_count(method, status))
                   for status in range(64) if self.lib.limber_timing_status_count(method, status)}
            for method, name in enumerate(('qag', 'cquad'))}
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--case', action='append', required=True,
                        help='Repeat for ordinary, known failing, or synthetic-input cases')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--shim', type=Path)
    parser.add_argument('--repetitions', type=int, default=9)
    parser.add_argument('--point-batch', type=int, default=256)
    parser.add_argument('--grid-batch', type=int, default=1)
    parser.add_argument('--workload', choices=['both', 'point', 'all_pairs_26ell'], default='both')
    parser.add_argument('--threads', type=int, choices=range(1, 9), default=1)
    args = parser.parse_args()
    if args.repetitions < 3 or min(args.point_batch, args.grid_batch) < 1:
        parser.error('Use at least three repetitions and positive batch sizes')
    if args.baseline_only:
        if args.shim or args.baseline or os.environ.get('DYLD_INSERT_LIBRARIES'):
            parser.error('Baseline must run without a shim or DYLD_INSERT_LIBRARIES')
    elif args.shim is None or args.baseline is None:
        parser.error('Instrumented run requires both --shim and --baseline')
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    os.environ['OMP_NUM_THREADS'] = str(args.threads)
    os.environ['OMP_DYNAMIC'] = 'FALSE'
    for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[key] = '1'
    import numpy as np
    import pyccl as ccl
    from common import digest, read_manifest, selected_case, write_json
    from ccl_helpers import make_cosmology, make_tracers, warning_records

    study = args.study.resolve()
    manifest, config = read_manifest(study)
    shim = None if args.baseline_only else Shim(args.shim)
    source = Path(__file__).resolve()
    binary = Path(ccl.__file__).parent / '_ccllib.so'
    omp = ct.CDLL(str(binary))
    omp.omp_set_num_threads.argtypes = [ct.c_int]
    omp.omp_get_max_threads.restype = ct.c_int
    omp.omp_set_num_threads(args.threads)
    if float(ccl.gsl_params.INTEGRATION_LIMBER_EPSREL) != 1e-4:
        raise RuntimeError('This timing protocol expects the installed 1e-4 Limber default')
    metadata = {
        'schema': 'native-ccl-integrator-timing-v1',
        'status': 'started', 'baseline_only': args.baseline_only,
        'ccl_version': ccl.__version__, 'ccl_binary_sha256': digest(binary),
        'python': sys.version, 'platform': platform.platform(),
        'script_sha256': digest(source),
        'shim_source_sha256': digest(source.with_name('native_timing_shim.c')),
        'shim_binary_sha256': digest(args.shim) if args.shim else None,
        'config_sha256': manifest['config_sha256'], 'threads': args.threads,
        'blas_threads': 1, 'repetitions': args.repetitions,
        'native_default_epsrel': 1e-4, 'strict_reference_epsrel': 1e-8,
        'scope': 'Warmed native CCL kernels; no Python integrand callbacks. Setup/CAMB excluded from both timing scopes.',
        'public_call_scope': 'Python angular_cl and CCL tracer collections, per-call QAG workspaces, and native projection. In CQUAD mode the outer QAG workspace is unused but still allocated.',
        'kernel_scope': 'Sum of native GSL integration-call durations, excluding workspace allocation. Includes native CQUAD fallbacks in QAG mode. With multiple threads this sum is not wall time.',
        'cquad_scope': 'Isolated native CQUAD diagnostic using cached per-thread workspace, not a released CCL3.3.3 API option or a benchmark of PR1313 split-QAG mitigation.',
        'instrumentation_control': 'Separate-process warmed public QAG repetitions without preloading. Instrumented/unmodified timing ratio includes process-to-process variation; no correction is subtracted from the paired QAG/CQUAD ratio.',
        'selected_workload': args.workload, 'cases': []}
    journal = args.output / 'timing.json'
    write_json(journal, metadata)
    try:
        for ident in args.case:
            case = selected_case(manifest, ident)
            input_path = study / 'inputs' / f'{case["family"]}.npz'
            if digest(input_path) != manifest['families'][case['family']]['npz_sha256']:
                raise ValueError('Input fingerprint mismatch')
            data = np.load(input_path)
            start = time.perf_counter()
            cosmo = make_cosmology(case)
            cosmo.compute_distances()
            cosmo.compute_nonlin_power()
            tracers, bias = make_tracers(cosmo, data['z'], data['nz'], data['centers'])
            setup = time.perf_counter() - start
            # CAMB may change the shared OpenMP runtime during setup. Restore
            # and verify the requested setting before any timed CCL call.
            setup_threads = omp.omp_get_max_threads()
            omp.omp_set_num_threads(args.threads)
            if omp.omp_get_max_threads() != args.threads:
                raise RuntimeError('Could not establish the requested native OpenMP thread count')
            point_pair = tuple(config['primary_pair_zero_based'])
            pairs = [(i, j) for i in range(5) for j in range(i, 5)]
            grid = np.geomspace(30., 3000., 26)
            workloads = {
                'point': ([(point_pair[0], point_pair[1])], np.array([config['ell']]), args.point_batch),
                'all_pairs_26ell': (pairs, grid, args.grid_batch)}
            if args.workload != 'both':
                workloads = {args.workload: workloads[args.workload]}
            case_record = {'case': case, 'input_sha256': digest(input_path),
                           'setup_seconds_not_benchmarked': setup, 'bias': bias.tolist(),
                           'omp_threads_after_setup_before_restore': setup_threads,
                           'verified_omp_threads_for_timing': omp.omp_get_max_threads(), 'workloads': {}}
            metadata['cases'].append(case_record)
            archive = {}
            base_data = None
            if not args.baseline_only:
                baseline_meta = json.loads((args.baseline / 'timing.json').read_text())
                if not baseline_meta['baseline_only'] or baseline_meta['status'] != 'completed':
                    raise ValueError('Need a completed no-preload baseline')
                for key in ('config_sha256', 'ccl_binary_sha256', 'threads', 'repetitions', 'selected_workload'):
                    if baseline_meta[key] != metadata[key]:
                        raise ValueError(f'Baseline mismatch: {key}')
                base_case = next(x for x in baseline_meta['cases'] if x['case']['id'] == ident)
                if base_case['case'] != case or base_case['input_sha256'] != digest(input_path):
                    raise ValueError('Baseline case mismatch')
                base_data = np.load(args.baseline / f'{ident}.npz')
                case_record['baseline_json_sha256'] = digest(args.baseline / 'timing.json')
                case_record['baseline_arrays_sha256'] = digest(args.baseline / f'{ident}.npz')

            for name, (pair_list, ell, batch) in workloads.items():
                def calculation():
                    return np.array([ccl.angular_cl(cosmo, tracers[i], tracers[j], ell,
                                    l_limber=-1, limber_integration_method='qag_quad')
                                    for i, j in pair_list])

                def evaluate(mode, count=1, override=0.):
                    problem = None
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter('always')
                        if shim:
                            shim.begin(mode, override)
                        try:
                            started = time.perf_counter_ns()
                            for _ in range(count):
                                result = calculation()
                            wall = time.perf_counter_ns() - started
                        except Exception as error:
                            problem = error
                            wall = time.perf_counter_ns() - started
                        finally:
                            if shim:
                                shim.end()
                    stats = shim.counters() if shim else {}
                    evidence = {'case': ident, 'workload': name, 'mode': mode,
                        'epsrel_override': override, 'batch': count, 'counters': stats,
                        'warnings': warning_records(caught)}
                    if problem:
                        metadata['last_failed_evaluation'] = dict(evidence, exception=repr(problem))
                        write_json(journal, metadata)
                        raise problem
                    if shim and (stats['unexpected_callbacks'] or not stats['qag_calls'] + stats['cquad_calls']):
                        metadata['last_failed_evaluation'] = evidence
                        raise RuntimeError(f'Shim did not isolate the native CCL callback: {stats}')
                    if shim:
                        expected_calls = len(pair_list)*len(ell)*count
                        actual_calls = stats['qag_calls'] if mode == 1 else stats['cquad_calls']
                        if actual_calls != expected_calls:
                            raise RuntimeError(f'Expected {expected_calls} native calls, got {actual_calls}')
                        if mode == 2 and (stats['cquad_failures'] or caught):
                            metadata['last_failed_evaluation'] = evidence
                            raise RuntimeError('CQUAD returned a failure or warning')
                    if not np.all(np.isfinite(result)) or np.any(result <= 0):
                        raise ValueError('Non-finite or non-positive density spectrum')
                    return result, {'public_seconds_per_workload': wall / 1e9 / count,
                        'kernel_seconds_per_workload': stats.get('kernel_nanoseconds', 0) / 1e9 / count,
                        'batch': count, 'counters': stats, 'warnings': warning_records(caught)}

                if args.baseline_only:
                    value, details = evaluate(0)
                    samples = []
                    for repetition in range(args.repetitions):
                        repeated_value, sample = evaluate(0, batch)
                        if not np.array_equal(repeated_value, value):
                            raise ValueError('Unmodified QAG changed between repetitions')
                        sample['repetition'] = repetition
                        samples.append(sample)
                    public_times = [s['public_seconds_per_workload'] for s in samples]
                    archive[name] = value
                    archive[name + '_ell'] = ell
                    archive[name + '_pairs'] = np.array(pair_list)
                    case_record['workloads'][name] = {
                        'baseline_values_only': False,
                        'scope': 'Warmed unmodified public angular_cl, no preloaded shim',
                        'warmup': details, 'samples': samples,
                        'summary': {'public_seconds_per_workload': {
                            'mean': statistics.mean(public_times),
                            'sample_scatter': statistics.stdev(public_times),
                            'median': statistics.median(public_times)}}}
                    write_json(journal, metadata)
                    np.savez_compressed(args.output/f'{ident}.npz', **archive)
                    continue
                for key, value in ((name + '_ell', ell), (name + '_pairs', np.array(pair_list))):
                    if not np.array_equal(base_data[key], value):
                        raise ValueError('Baseline estimator coordinates differ')
                baseline_samples = base_case['workloads'][name]['samples']
                if len(baseline_samples) != args.repetitions or any(s['batch'] != batch for s in baseline_samples):
                    raise ValueError('No-preload QAG repetitions/batch differ from the instrumented protocol')
                # Mode1 first establishes the native callback address guard.
                qag, warm_qag = evaluate(1)
                if not np.array_equal(qag, base_data[name]):
                    raise ValueError('Instrumented QAG is not bitwise equal to unmodified CCL')
                cquad, warm_cquad = evaluate(2)
                reference, reference_check = evaluate(2, override=1e-8)
                case_record['workloads'][name] = {'warm_qag': warm_qag,
                    'warm_cquad': warm_cquad, 'strict_reference': reference_check}
                if reference_check['counters']['cquad_failures'] or reference_check['warnings']:
                    raise ValueError('Strict native CQUAD reference did not return cleanly')
                accepted_control = None
                if name == 'point':
                    summary_path = study/'summary.json'
                    saved = json.loads(summary_path.read_text())
                    accepted = next(r for r in saved['rows'] if r['case']['id'] == ident)
                    if not accepted['ccl_reference_valid']:
                        raise ValueError('The independent accuracy reference is not validated')
                    difference = abs(float(reference[0,0])/accepted['common_CCL_CQUAD_reference_value']-1)
                    if difference > config['reference_consistency_fractional_tolerance']:
                        raise ValueError(f'Native timing reference differs from accuracy archive: {difference}')
                    accepted_control = {'summary_sha256': digest(summary_path),
                        'reference_value': accepted['common_CCL_CQUAD_reference_value'],
                        'fractional_difference': difference, 'basis': accepted['ccl_reference_basis']}
                expected = {'qag_quad': qag, 'cquad_diagnostic': cquad}
                records = {key: [] for key in expected}
                for repetition in range(args.repetitions):
                    order = list(expected) if repetition % 2 == 0 else list(reversed(expected))
                    for method in order:
                        value, sample = evaluate(1 if method == 'qag_quad' else 2, batch)
                        if not np.array_equal(value, expected[method]):
                            raise ValueError(f'{method} result changed between repetitions')
                        if sample['counters']['workspace_allocations']:
                            metadata['last_failed_evaluation'] = dict(case=ident, workload=name, method=method, sample=sample)
                            raise RuntimeError('CQUAD workspace allocation occurred inside a measured repetition')
                        sample['repetition'] = repetition
                        sample['position_in_interleaving'] = order.index(method)
                        records[method].append(sample)
                report = {'shape': list(reference.shape), 'multipoles': ell.tolist(),
                    'pairs': pair_list, 'baseline_qag_bitwise_equal': True,
                    'warm_qag': warm_qag, 'warm_cquad': warm_cquad,
                    'strict_reference': reference_check, 'samples': records,
                    'accepted_accuracy_reference_check': accepted_control,
                    'reference_scope': ('Verified against the accepted accuracy-study reference' if name == 'point'
                        else 'Internal strict-CQUAD comparison at additional multipoles; not independently certified convergence'),
                    'max_fractional_error_vs_strict_reference': {
                        method: float(np.max(np.abs(value / reference - 1)))
                        for method, value in expected.items()}, 'summary': {}}
                for method, samples in records.items():
                    report['summary'][method] = {}
                    for metric in ('public_seconds_per_workload', 'kernel_seconds_per_workload'):
                        numbers = [s[metric] for s in samples]
                        report['summary'][method][metric] = {
                            'mean': statistics.mean(numbers),
                            'sample_scatter': statistics.stdev(numbers),
                            'median': statistics.median(numbers)}
                report['cquad_over_qag_ratio'] = {
                    metric: report['summary']['cquad_diagnostic'][metric]['mean'] /
                            report['summary']['qag_quad'][metric]['mean']
                    for metric in ('public_seconds_per_workload', 'kernel_seconds_per_workload')}
                unmodified = base_case['workloads'][name]['summary']['public_seconds_per_workload']
                report['unmodified_qag_control'] = {
                    'public_seconds_per_workload': unmodified,
                    'repetitions': len(base_case['workloads'][name]['samples']),
                    'batch': base_case['workloads'][name]['samples'][0]['batch'],
                    'instrumented_over_unmodified_mean_ratio':
                        report['summary']['qag_quad']['public_seconds_per_workload']['mean'] / unmodified['mean'],
                    'interpretation': 'Instrumentation-overhead control from a separate fresh process; includes process-to-process timing variation. Primary method comparison remains the interleaved instrumented QAG/CQUAD ratio.'}
                report['interpretation'] = 'Ratio is cost at the same requested tolerance; when QAG is inaccurate it is not a speed comparison at matched achieved accuracy.'
                case_record['workloads'][name] = report
                archive[name + '_qag'] = qag
                archive[name + '_cquad'] = cquad
                archive[name + '_reference'] = reference
                archive[name + '_ell'] = ell
                archive[name + '_pairs'] = np.array(pair_list)
                np.savez_compressed(args.output/f'{ident}.npz', **archive)
                write_json(journal, metadata)
            arrays = args.output / f'{ident}.npz'
            np.savez_compressed(arrays, **archive)
            case_record['arrays_sha256'] = digest(arrays)
            write_json(journal, metadata)
        metadata['status'] = 'completed'
    except Exception as error:
        metadata.update(status='failed', exception=repr(error), traceback=traceback.format_exc())
        raise
    finally:
        if shim:
            shim.end()
        write_json(journal, metadata)


if __name__ == '__main__':
    main()
