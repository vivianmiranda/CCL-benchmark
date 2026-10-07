"""Check narrow-overlap references with physical-support CQUAD subdivisions.

Every subinterval uses native GSL CQUAD and the same native CCL integrand.
The original whole-domain limits remain; no tail is discarded. This is a
reference diagnostic, never a production-runtime benchmark.
"""
import argparse
import json
import math
from pathlib import Path
import time
import traceback
import warnings

import numpy as np
import pyccl as ccl
from common import digest, read_manifest, selected_case, write_json, source_hashes
from ccl_helpers import GSLRules, make_cosmology, make_tracers, warning_records


def partition(context, cosmo, ell, mean, sigma, sigmas):
    redshifts = np.array([mean+s*sigma for s in sigmas])
    if np.any(redshifts <= 0):
        raise ValueError('The requested support partitions reach nonpositive redshift')
    chi = ccl.comoving_radial_distance(cosmo, 1/(1+redshifts))
    logk = np.log((ell+.5)/chi)
    internal = logk[(logk > context['lo']) & (logk < context['hi'])]
    edges = np.unique(np.r_[context['lo'], internal, context['hi']])
    if len(edges) != len(sigmas)+2:
        raise ValueError('A narrow-population support boundary lies outside the native integration interval')
    return edges, {'redshifts': redshifts.tolist(), 'logk': logk.tolist(), 'sigma_offsets': list(sigmas)}


def integrate_support(rules, context, edges):
    controls = []
    for left, right in zip(edges[:-1], edges[1:]):
        interval = dict(context, lo=float(left), hi=float(right))
        record = rules.integrate(interval, 'cquad', 1e-8)
        record['logk_bounds'] = [float(left), float(right)]
        controls.append(record)
    successful = all(r['value'] is not None and np.isfinite(r['value']) and
                     not any(r['status']) and not r['callback_errors'] for r in controls)
    value = math.fsum(r['value'] for r in controls) if all(r['value'] is not None for r in controls) else None
    return {'method': 'native GSL CQUAD on physically partitioned native CCL integrand',
            'epsrel_per_interval': 1e-8, 'native_domain_retained': True,
            'edges_logk': edges.tolist(), 'intervals': controls, 'value': value,
            'successful': bool(successful and value is not None and value > 0),
            'calls': sum(r['calls'] for r in controls)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--case', required=True)
    args = parser.parse_args(); output = args.output.resolve()
    manifest, config = read_manifest(output); case = selected_case(manifest, args.case)
    definition = manifest['families'][case['family']]['definition']
    if definition['kind'] != 'positive_mixture':
        raise ValueError('Physical narrow-population controls require an explicit positive-mixture input')
    original_path = output/'ccl'/f'{case["id"]}.json'
    original = json.loads(original_path.read_text())
    if original['status'] != 'completed' or original['case'] != case:
        raise ValueError('The native CCL case must finish first')
    destination = output/'support_references'/f'{case["id"]}.json'
    if destination.exists():
        raise FileExistsError(destination)
    supplemental_path = output/'reference_supplements'/f'{case["id"]}.json'
    supplement = json.loads(supplemental_path.read_text()) if supplemental_path.exists() else None
    if supplement and supplement['original_result_sha256'] != digest(original_path):
        raise ValueError('Existing supplement no longer matches its original record')
    start = time.perf_counter()
    record = {'case': case, 'status': 'started', 'original_result_sha256': digest(original_path),
              'config_sha256': manifest['config_sha256'], 'script_sha256': source_hashes(),
              'input_definition': definition, 'references': [], 'elapsed_is_not_benchmark': True}
    write_json(destination, record)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            inp = output/'inputs'/f'{case["family"]}.npz'
            fine_path = output/'inputs'/f'{case["family"]}_refined.npz'
            family = manifest['families'][case['family']]
            if digest(inp) != family['npz_sha256'] or digest(fine_path) != family['refined_npz_sha256']:
                raise ValueError('Saved input grid changed')
            data, fine = np.load(inp), np.load(fine_path)
            cosmo = make_cosmology(case)
            tracers, _ = make_tracers(cosmo, data['z'], data['nz'], data['centers'])
            fine_tracers, _ = make_tracers(cosmo, fine['z'], fine['nz'], fine['centers'])
            rules = GSLRules(); mean, sigma = definition['shared_mean_z'], definition['shared_sigma_z']
            tol = config['reference_consistency_fractional_tolerance']
            for old in original['references']:
                pair = old['pair']; context = rules.context(cosmo, tracers, pair, config['ell'])
                edges, physical = partition(context, cosmo, config['ell'], mean, sigma, [-6, 0, 6])
                dense_edges, dense_physical = partition(context, cosmo, config['ell'], mean, sigma, [-6, -3, -1, 0, 1, 3, 6])
                whole = integrate_support(rules, context, edges)
                refined = integrate_support(rules, context, dense_edges)
                fc = rules.context(cosmo, fine_tracers, pair, config['ell'])
                fe, _ = partition(fc, cosmo, config['ell'], mean, sigma, [-6, -3, -1, 0, 1, 3, 6])
                fine_control = integrate_support(rules, fc, fe)
                accepted = next((r for r in (supplement or {}).get('references', []) if r['pair'] == pair), old)
                original_value = accepted.get('reference_value', accepted.get('cquad', {}).get('value'))
                value = refined['value']
                changes = {}
                for key, item in [('coarse_support', whole), ('fine_input', fine_control)]:
                    changes[key] = abs(item['value']/value-1) if item['successful'] and refined['successful'] else None
                changes['original_reference'] = abs(original_value/value-1) if original_value is not None and value is not None and value > 0 else None
                valid = bool(whole['successful'] and refined['successful'] and fine_control['successful'] and
                             changes['coarse_support'] <= tol and
                             changes['fine_input'] <= config['input_refinement_fractional_tolerance'])
                confirmed = bool(valid and accepted['reference_valid'] and
                                 changes['original_reference'] is not None and changes['original_reference'] <= tol)
                record['references'].append({'pair': pair, 'primary': old['primary'],
                    'physical_support': physical, 'refined_physical_support': dense_physical,
                    'support_CQUAD': whole, 'support_CQUAD_refined': refined,
                    'fine_input_support_CQUAD': fine_control,
                    'support_reference_value': value, 'support_reference_valid': valid,
                    'original_reference_confirmed': confirmed, 'fractional_comparisons': changes,
                    'reference_tolerance': tol, 'input_refinement_tolerance': config['input_refinement_fractional_tolerance']})
                write_json(destination, record)
            record['status'] = 'completed'
        record['warnings'] = warning_records(caught)
    except Exception as error:
        record.update(status='exception', exception=repr(error), traceback=traceback.format_exc())
        raise
    finally:
        if 'caught' in locals():
            record['warnings'] = warning_records(caught)
        record['elapsed_seconds_for_resource_planning'] = time.perf_counter()-start
        write_json(destination, record)
    print(json.dumps({'case': case['id'], 'status': record['status'],
        'pairs': len(record['references']), 'confirmed': sum(r['original_reference_confirmed'] for r in record['references'])}, indent=2))


if __name__ == '__main__':
    main()
