"""Archive finished paired evidence and plots, excluding large power bundles.

This collector only reads saved numerical results. It does not compute spectra,
compile code or overwrite an existing archive. Reference failures are preserved
and counted; a completed calculation is not automatically a validated reference.
"""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from common import PACKAGE, digest, read_manifest, source_hashes, write_json
from summarize import comparison_row


def read(path):
    return json.loads(path.read_text())


def verify_ladder(folder, case, manifest, cocoa):
    """Audit saved refinement values and the recorded convergence checks."""
    path = folder/'cocoa_integration_ladder'/f'{case["id"]}.json'
    if not path.exists():
        return set()
    record = read(path); arrays_path = path.with_suffix('.npz')
    if (record['status'] != 'completed' or not record.get('integration_ladder') or
            record['case'] != case or record['config_sha256'] != manifest['config_sha256']):
        raise ValueError(f'Unfinished or mismatched integration ladder: {path}')
    for key in ('shared_power_sha256', 'nz_text_sha256', 'core_commit', 'lsst_y1_commit'):
        if record[key] != cocoa[key]:
            raise ValueError(f'Ladder changed its base inputs/code ({key}): {path}')
    nodes = [96, 128, 256, 512, 1024]
    controls = record['controls']
    if (record['nodes'] != nodes or len(controls) != 5 or
            [r['integration_accuracy'] for r in controls] != list(range(5)) or
            [r['radial_nquad'] for r in controls] != nodes or
            any(r['accuracy_boost'] != 1 for r in controls) or
            len({r['nwindow'] for r in controls}) != 1 or controls[0]['nwindow'] <= 0):
        raise ValueError(f'Unexpected integration-accuracy controls: {path}')
    with np.load(arrays_path) as saved, np.load(folder/'cocoa'/f'{case["id"]}.npz') as base:
        values = saved['spectra'].copy()
        if values.shape != (5, 15) or not np.all(np.isfinite(values)) or np.any(values <= 0):
            raise ValueError(f'Invalid ladder arrays: {arrays_path}')
        for key in ('ell', 'pairs', 'a_edges'):
            if not np.array_equal(saved[key], base[key]):
                raise ValueError(f'Ladder changed estimator coordinates: {key}')
        if not np.array_equal(saved['nodes'], nodes):
            raise ValueError('Ladder NPZ node counts differ from JSON')
        for pair, column in zip(saved['pairs'], values.T):
            row = next(r for r in record['pairs'] if r['pair'] == pair.tolist())
            if not np.array_equal(column, row['values']):
                raise ValueError(f'Ladder JSON/NPZ values disagree: {pair}')
        for i, n in enumerate(base['nodes']):
            if int(n) not in nodes or not np.array_equal(values[nodes.index(int(n))], base['spectra'][i]):
                raise ValueError(f'Previously measured ladder level changed: {n}')
    files = {path, arrays_path}
    regression_path = path.parent/'regression_check.json'
    if regression_path.exists():
        check = read(regression_path)
        if (check['source_sha256'] != digest(path) or not check['successful'] or check['tests_run'] != 3 or
                check['convergence_rtol'] != 1e-4 or check['reference_rtol'] != 1e-5 or check['atol'] != 0):
            raise ValueError(f'Unexpected regression-check provenance or tolerances: {regression_path}')
        maxima = np.max(np.abs(values/values[-1]-1), axis=1)
        if not np.array_equal(maxima, check['max_fractional_change_all_pairs_by_level']):
            raise ValueError('Recorded regression changes differ from saved spectra')
        if not (maxima[0] > check['convergence_rtol'] and maxima[2] <= check['convergence_rtol'] and
                maxima[3] <= check['reference_rtol']):
            raise ValueError('Saved spectra do not reproduce the three recorded convergence outcomes')
        files.add(regression_path)
    return files


def verify_targeted_figures(folder):
    """Keep targeted figures distinct and verify their saved source hashes."""
    figures = folder/'figures'
    if not figures.exists():
        return []
    manifest, _ = read_manifest(folder)
    cases = {c['id']: c for c in manifest['cases']}
    for path in sorted(figures.glob('*.json')):
        meta = read(path)
        if path.name == 'plot_index.json':
            for source in meta.get('source_summaries', []):
                if digest(source['path']) != source['sha256']:
                    raise ValueError(f'Targeted figure has stale summary: {path}')
            continue
        if isinstance(meta.get('source_sha256'), dict):
            case = cases[meta['case']]; ident, family = case['id'], case['family']
            sources = {'summary': folder/'summary.json', 'inputs': folder/'inputs'/f'{family}.npz',
                'integrand': folder/'ccl'/f'{ident}_integrand.npz', 'cocoa': folder/'cocoa'/f'{ident}.npz',
                'integration_ladder': folder/'cocoa_integration_ladder'/f'{ident}.json'}
            for key, fingerprint in meta['source_sha256'].items():
                if key not in sources or digest(sources[key]) != fingerprint:
                    raise ValueError(f'Targeted figure source changed: {path}, {key}')
        elif path.name == 'integration_accuracy_convergence.json':
            matches = [p for p in (folder/'cocoa_integration_ladder').glob('*.json') if digest(p) == meta['source_sha256']]
            if len(matches) != 1:
                raise ValueError(f'Convergence figure source not identified uniquely: {path}')
            ladder = read(matches[0]); values = np.asarray([r['values'] for r in ladder['pairs']])
            primary = next(r for r in ladder['pairs'] if r['pair'] == [3, 4])
            if (meta['levels'] != list(range(5)) or meta['nodes'] != ladder['nodes'] or
                    not np.array_equal(meta['primary_error_percent'], 100*np.abs(np.asarray(primary['values'])/primary['values'][-1]-1)) or
                    not np.array_equal(meta['worst_pair_error_percent'], 100*np.max(np.abs(values/values[:, -1, None]-1), axis=0))):
                raise ValueError(f'Convergence figure metadata differ from saved spectra: {path}')
        else:
            raise ValueError(f'Unknown targeted figure provenance schema: {path}')
        if not all(path.with_suffix(ext).exists() for ext in ('.png', '.pdf')):
            raise ValueError(f'Missing targeted figure rendering: {path}')
    files = sorted(p for p in figures.iterdir() if p.suffix in ('.png', '.pdf', '.json'))
    for path in files:
        if path.suffix in ('.png', '.pdf') and not path.with_suffix('.json').exists() and not (figures/'plot_index.json').exists():
            raise ValueError(f'Targeted figure has no provenance record: {path}')
    return files


def verify_study(folder):
    manifest, config = read_manifest(folder)
    cases = manifest['cases']
    summary = read(folder/'summary.json') if (folder/'summary.json').exists() else None
    if summary is None:
        if cases:
            raise FileNotFoundError(f'Missing result summary: {folder}')
        summary = {'config_sha256': manifest['config_sha256'], 'rows': []}
    if summary['config_sha256'] != manifest['config_sha256']:
        raise ValueError(f'Stale summary: {folder}')
    if {r['case']['id'] for r in summary['rows']} != {c['id'] for c in cases}:
        raise ValueError(f'Summary does not contain every declared case: {folder}')
    files = {folder/'manifest.json', folder/'config.json'}
    if (folder/'summary.json').exists():
        files.add(folder/'summary.json')
    for name in ('search.json', 'search_progress.json'):
        if (folder/name).exists():
            files.add(folder/name)
    families = {c['family'] for c in cases}
    # Empty searches retain their own search record, but have no result/input claim.
    for family in families:
        record = manifest['families'][family]
        for suffix, key in [('.npz', 'npz_sha256'), ('.txt', 'text_sha256'), ('_refined.npz', 'refined_npz_sha256')]:
            if key in record:
                path = folder/'inputs'/f'{family}{suffix}'
                if digest(path) != record[key]:
                    raise ValueError(f'Changed input: {path}')
                files.add(path)
    omitted_power = {}
    reference_counts = {'primary_CCL_valid': 0, 'primary_CoCoA_valid': 0,
                        'primary_CCL_CQUAD_only_basis': 0, 'primary_CCL_successful_QAG_control_basis': 0}
    for case in cases:
        ident = case['id']; cpath = folder/'ccl'/f'{ident}.json'; apath = folder/'cocoa'/f'{ident}.json'
        ccl, cocoa = read(cpath), read(apath)
        for code, record in [('CCL', ccl), ('CoCoA', cocoa)]:
            if record['status'] != 'completed' or record['case'] != case or record['config_sha256'] != manifest['config_sha256']:
                raise ValueError(f'Unfinished or mismatched {code} case: {ident}')
        family = manifest['families'][case['family']]
        if ccl['input_sha256'] != family['npz_sha256'] or cocoa['nz_text_sha256'] != family['text_sha256']:
            raise ValueError(f'Input mismatch: {ident}')
        power = ccl['shared_power']; path = folder/power['path']
        if power['sha256'] != cocoa['shared_power_sha256']:
            raise ValueError(f'Codes received different power bundles: {ident}')
        if power['path'] not in omitted_power:
            if digest(path) != power['sha256']:
                raise ValueError(f'Changed local power bundle: {path}')
            metadata = path.with_suffix('.json')
            if read(metadata)['sha256'] != power['sha256']:
                raise ValueError(f'Power metadata mismatch: {metadata}')
            files.add(metadata)
            omitted_power[power['path']] = {'sha256': power['sha256'], 'bytes': path.stat().st_size,
                'reason': 'Large generated power/background bundle retained locally; not included in Git.',
                'regeneration': 'run_ccl.py exports this bundle from the recorded cosmology with the recorded CCL/CAMB setup.'}
        ap = folder/'cocoa'/f'{ident}.npz'; ip = folder/'ccl'/f'{ident}_integrand.npz'
        with np.load(ap) as arrays:
            for pair, values in zip(arrays['pairs'], arrays['spectra'].T):
                row = next(r for r in cocoa['pairs'] if r['pair'] == pair.tolist())
                if not np.array_equal(values, np.asarray(row['values'], dtype=float), equal_nan=True):
                    raise ValueError(f'CoCoA saved values differ between JSON and NPZ: {ident}')
        with np.load(ip) as arrays:
            if not all(np.all(np.isfinite(arrays[k])) for k in arrays.files):
                raise ValueError(f'Non-finite saved integrand: {ident}')
        supplement_path = folder/'reference_supplements'/f'{ident}.json'
        supplement = read(supplement_path) if supplement_path.exists() else None
        if supplement:
            if supplement['original_result_sha256'] != digest(cpath):
                raise ValueError(f'Supplement does not match original result: {ident}')
            files.add(supplement_path)
        support_path = folder/'support_references'/f'{ident}.json'
        support = read(support_path) if support_path.exists() else None
        if support:
            if (support['original_result_sha256'] != digest(cpath) or support['case'] != case or
                    support['config_sha256'] != manifest['config_sha256']):
                raise ValueError(f'Support reference does not match original result: {ident}')
            files.add(support_path)
        files.update(verify_ladder(folder, case, manifest, cocoa))
        primary = config['primary_pair_zero_based']
        expected = comparison_row(case, primary, ccl, cocoa, supplement, True, support)
        actual = next(r for r in summary['rows'] if r['case']['id'] == ident)
        if actual != expected:
            raise ValueError(f'Out-of-date primary summary: {ident}')
        reference_counts['primary_CCL_valid'] += expected['ccl_reference_valid']
        reference_counts['primary_CoCoA_valid'] += expected['cocoa_reference_valid']
        cquad_only = 'no successful split-QAG control' in expected['ccl_reference_basis']
        if expected['ccl_reference_valid']:
            reference_counts['primary_CCL_CQUAD_only_basis' if cquad_only else 'primary_CCL_successful_QAG_control_basis'] += 1
        pairs = {tuple(primary)} | {tuple(r['pair']) for r in ccl.get('references', [])}
        for pair in pairs:
            expected_pair = comparison_row(case, list(pair), ccl, cocoa, supplement, list(pair) == primary, support)
            actual_pair = next(r for r in summary['pair_rows'] if r['case']['id'] == ident and r['pair'] == list(pair))
            if actual_pair != expected_pair:
                raise ValueError(f'Out-of-date pair summary: {ident}, {pair}')
        files.update([cpath, apath, ap, ip])
    return files, {'cases': len(cases), 'families_with_published_inputs': sorted(families),
                   'reference_counts': reference_counts, 'omitted_large_power_bundles': omitted_power}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='Finished uniform study folder')
    parser.add_argument('--targeted-output', type=Path, action='append', default=[])
    parser.add_argument('--destination', type=Path, required=True, help='New permanent result archive; must not exist')
    parser.add_argument('--figures-destination', type=Path, required=True, help='New permanent figure folder; must not exist')
    args = parser.parse_args(); base = args.output.resolve()
    destination = args.destination.resolve(); figures_destination = args.figures_destination.resolve()
    if destination.exists() or figures_destination.exists():
        raise FileExistsError('Use fresh result and figure destinations; existing archives are immutable')
    studies = [(base, Path('.'))]
    for folder in args.targeted_output:
        folder = folder.resolve()
        studies.append((folder, Path('targeted')/folder.name))
    if len({str(relative) for _, relative in studies}) != len(studies):
        raise ValueError('Duplicate targeted destination names')
    # Finish every input/result check before creating a publication directory.
    checked = [(folder, relative, *verify_study(folder)) for folder, relative in studies]
    source_figures = base/'figures'
    if not (source_figures/'plot_index.json').exists():
        raise FileNotFoundError('Generate and visually inspect the final paired figures first')
    figure_index = read(source_figures/'plot_index.json')
    plotted = figure_index.get('source_summaries', [])
    expected_summaries = {str((folder/'summary.json').resolve()) for folder, _, _, record in checked if record['cases']}
    if {r['path'] for r in plotted} != expected_summaries:
        raise ValueError('Final paired figures do not cover exactly the declared nonempty studies')
    for record in plotted:
        if digest(record['path']) != record['sha256']:
            raise ValueError(f'Summary changed after plotting: {record["path"]}')
    figure_files = sorted(p for p in source_figures.iterdir() if p.suffix in ('.png', '.pdf', '.json'))
    figure_copies = [(p, p.name) for p in figure_files]
    for folder, relative in studies[1:]:
        figure_copies.extend((p, f'{folder.name}__{p.name}') for p in verify_targeted_figures(folder))
    if len({name for _, name in figure_copies}) != len(figure_copies):
        raise ValueError('Targeted figure prefixes would overwrite a figure')
    destination.mkdir(parents=True); figures_destination.mkdir(parents=True)
    records = []
    for folder, relative, files, record in checked:
        for src in sorted(files):
            dst = destination/relative/src.relative_to(folder)
            dst.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(src, dst)
            if digest(src) != digest(dst):
                raise ValueError(f'Copy changed file: {src}')
        records.append({'source_folder': str(folder), 'archive_relative_folder': str(relative), **record})
    for src, name in figure_copies:
        shutil.copy2(src, figures_destination/name)
        if digest(src) != digest(figures_destination/name):
            raise ValueError(f'Copied figure differs from source: {src}')
    shutil.copy2(PACKAGE/'PLAN.md', destination/'sampling_and_acceptance_plan.md')
    archive = {'schema': 1, 'studies': records, 'script_sha256_at_collection': source_hashes(),
        'included_result_files': {str(p.relative_to(destination)): digest(p) for p in sorted(destination.rglob('*')) if p.is_file()},
        'included_figure_files': {p.name: digest(p) for p in sorted(figures_destination.iterdir()) if p.is_file()},
        'scope': 'Saved native CCL spectra/reference controls and actual CoCoA covariance-module spectra. No covariance matrices. Completed does not imply a passing reference; validity flags and failed controls are retained.',
        'power_bundle_policy': 'Shared generated power/background NPZs are deliberately omitted. Their metadata and SHA256 remain archived. Numerical reruns must regenerate them; the plotted spectra, redshift inputs and reference evidence are included.'}
    write_json(destination/'publication_manifest.json', archive)
    print(json.dumps({'destination': str(destination), 'figures': str(figures_destination),
                      'result_files': len(archive['included_result_files']), 'figure_files': len(archive['included_figure_files'])}, indent=2))


if __name__ == '__main__':
    main()
