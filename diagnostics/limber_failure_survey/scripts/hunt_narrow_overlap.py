"""Move a narrow n(z) population across quadrature nodes at fixed cosmology.

These deliberately selected roots are stress tests, never random incidence
samples. Both codes subsequently receive the saved identical distributions.
The search only selects candidates; full native and refined integrals must
still validate whether each candidate is a failure.
"""
import argparse
import copy
from pathlib import Path
import shutil
import time

import numpy as np
from scipy.optimize import brentq

from common import build_family, digest, read_manifest, write_json, source_hashes
from ccl_helpers import GSLRules, make_cosmology, make_tracers


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--width', type=float, choices=[.003, .01, .03], required=True)
    p.add_argument('--maximum-roots', type=int, default=2)
    a = p.parse_args()
    base = a.output.resolve()
    manifest, config = read_manifest(base)
    tag = f'{a.width:g}'.replace('.', 'p')
    folder = base/'targeted'/f'narrow_width_{tag}'
    if folder.exists():
        raise FileExistsError(folder)
    folder.mkdir(parents=True)
    (folder/'inputs').mkdir()
    shutil.copyfile(base/'config.json', folder/'config.json')
    parameters = dict(config['fixed_cosmology'], **config['known_control'])
    case = {'parameters': parameters,
            'camb_dark_energy_model': config['survey_camb_dark_energy_model']}
    cosmo = make_cosmology(case)
    rules = GSLRules()
    definition = dict(kind='positive_mixture', shared_fraction=.15,
                      shared_mean_z=1., shared_sigma_z=a.width,
                      z_nodes=20001, refined_z_nodes=40001)
    samples, roots, rejected = {}, [], []
    start = time.perf_counter()

    def evaluate(mean):
        key = float(mean)
        if key in samples:
            return samples[key]['fractional_rule_difference']
        if len(samples) >= 160 or time.perf_counter()-start > 180:
            raise TimeoutError('Narrow-overlap search reached its 160-evaluation/180-second budget')
        current = dict(definition, shared_mean_z=key)
        z, nz, centers = build_family(current)
        tracers, _ = make_tracers(cosmo, z, nz, centers)
        context = rules.context(cosmo, tracers, config['primary_pair_zero_based'], config['ell'])
        initial, nodes, _ = rules.initial_rules(context)
        assert not initial['callback_errors'] and initial['kronrod'] > 0
        value = (initial['kronrod']-initial['gauss'])/initial['kronrod']
        samples[key] = dict(initial, mean_z=key, fractional_rule_difference=value,
                            logk_bounds=[context['lo'], context['hi']],
                            kronrod_logk=nodes.tolist())
        write_json(folder/'search_progress.json', {
            'evaluations': list(samples.values()), 'candidate_means': roots,
            'excluded_from_incidence': True})
        return value

    status, problem = 'completed', None
    try:
        previous = None
        for mean in np.linspace(.8, 1.2, 41):
            value = evaluate(mean)
            if previous and previous[1]*value < 0:
                try:
                    root = float(brentq(evaluate, previous[0], mean, xtol=1e-10, maxiter=40))
                    evaluate(root)
                    candidate = samples[root]
                    if candidate['estimated_error'] <= 1e-4*abs(candidate['kronrod']):
                        if not roots or abs(root-roots[-1]) > 1e-7:
                            roots.append(root)
                    else:
                        rejected.append({'bracket': [previous[0], float(mean)],
                                         'reason': 'Rule crossing did not satisfy GSL initial acceptance',
                                         'candidate': candidate})
                except RuntimeError as error:
                    rejected.append({'bracket': [previous[0], float(mean)], 'reason': repr(error)})
                if len(roots) >= a.maximum_roots:
                    break
            previous = (float(mean), value)
    except Exception as error:
        status, problem = 'stopped', repr(error)

    new = copy.deepcopy(manifest)
    new.update(cases=[], pilot_case_ids=[], families={},
               sampling_description='Node-guided synthetic peak-position search at fixed cosmology; excluded from all incidence estimates.')
    for index, mean in enumerate(roots):
        name = f'narrow_{tag}_root{index}'
        current = dict(definition, name=name, shared_mean_z=mean)
        z, nz, centers = build_family(current)
        np.savez_compressed(folder/'inputs'/f'{name}.npz', z=z, nz=nz, centers=centers)
        np.savetxt(folder/'inputs'/f'{name}.txt', np.column_stack([z, *nz]))
        fz, fnz, fc = build_family(current, refined=True)
        np.savez_compressed(folder/'inputs'/f'{name}_refined.npz', z=fz, nz=fnz, centers=fc)
        new['families'][name] = {
            'definition': current, 'nodes': len(z), 'minimum_value': float(nz.min()),
            'npz_sha256': digest(folder/'inputs'/f'{name}.npz'),
            'text_sha256': digest(folder/'inputs'/f'{name}.txt'),
            'refined_npz_sha256': digest(folder/'inputs'/f'{name}_refined.npz')}
        ident = f'peak_{tag}_{index}'
        new['cases'].append(dict(case, id=ident, cosmology_id='fixed_ppf', family=name,
                                 arm='targeted_narrow_overlap'))
        new['pilot_case_ids'].append(ident)
    write_json(folder/'manifest.json', new)
    write_json(folder/'search.json', {
        'status': status, 'problem': problem, 'parameters': case,
        'peak_definition': definition, 'mean_z_range': [.8, 1.2], 'bracketing_nodes': 41,
        'maximum_roots': a.maximum_roots, 'evaluations': list(samples.values()),
        'rejected_brackets': rejected,
        'candidate_means': roots, 'excluded_from_incidence': True,
        'script_sha256': source_hashes(),
        'elapsed_seconds_for_resource_planning': time.perf_counter()-start,
        'elapsed_is_not_benchmark': True})
    print(f'{folder}: {len(roots)} candidate roots; {status}; {problem}', flush=True)
    if problem:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
