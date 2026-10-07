"""Summarize paired saved cases; preserve exclusions and targeted points."""
import argparse
import json
import math
from pathlib import Path

from common import read_manifest, write_json, digest


def load(path):
    return json.loads(path.read_text()) if path.exists() else None


def wilson(count, number):
    """Two-sided 95% Wilson score interval for a binomial proportion."""
    if not number:
        return None
    z = 1.959963984540054
    p = count / number
    scale = 1 + z * z / number
    center = (p + z * z / (2 * number)) / scale
    half = z * math.sqrt(p * (1 - p) / number + z * z / (4 * number * number)) / scale
    return [0. if count == 0 else max(0., center - half),
            1. if count == number else min(1., center + half)]


def comparison_row(case, pair, ccl, cocoa, supplement, primary, support=None):
    ccl, cocoa, supplement = ccl or {}, cocoa or {}, supplement or {}
    cr = next((r for r in ccl.get('references', []) if r['pair'] == pair), None)
    ar = next((r for r in cocoa.get('pairs', []) if r['pair'] == pair), None)
    refined = next((r for r in supplement.get('references', []) if r['pair'] == pair), None)
    accepted = refined if refined else cr
    refvalue = accepted.get('reference_value', accepted.get('cquad', {}).get('value')) if accepted else None
    extra = ((accepted or {}).get('supplementary_reference_checks')
             or (accepted if accepted and 'cquad_refined' in accepted else None))
    if extra:
        basis = ('two refined CQUAD controls plus successful split-QAG controls' if extra['successful_split_QAG_controls']
                 else 'two refined CQUAD controls; no successful split-QAG control')
    else:
        basis = 'CQUAD plus split-QAG controls' if accepted else 'no reference computed'
    support_check = next((r for r in (support or {}).get('references', []) if r['pair'] == pair), None)
    support_required = case['arm'] == 'targeted_narrow_overlap'
    support_confirmed = bool(support and support.get('status') == 'completed' and
                             support_check and support_check['original_reference_confirmed'])
    if support_required or support is not None:
        basis += '; physical-support CQUAD confirmed' if support_confirmed else '; physical-support CQUAD pending or failed'
    if support_confirmed:
        refvalue = support_check['support_reference_value']
    native = next((r['qag_quad'] for r in ccl.get('native', []) if r['pair'] == pair), {})
    replay = (cr or {}).get('native_qag_replay', {})
    silent = bool(native.get('value') is not None and not native.get('warnings') and not native.get('exception')
                  and replay.get('status') and not any(replay['status']) and not replay.get('callback_errors'))
    reported = bool(native.get('warnings') or native.get('exception') or any(replay.get('status', [])) or replay.get('callback_errors'))
    def offset(value):
        return value / refvalue - 1 if value is not None and refvalue is not None and refvalue > 0 else None
    return {'case': case, 'pair': pair, 'primary': primary,
        'ccl_status': ccl.get('status', 'not_run'), 'cocoa_status': cocoa.get('status', 'not_run'),
        'ccl_reference_valid': bool(accepted and accepted['reference_valid'] and
                                    (support_confirmed if support_required or support is not None else True)),
        'ccl_reference_basis': basis,
        'ccl_support_reference_required': support_required,
        'ccl_support_reference_checks': support_check,
        'ccl_silent_success': silent, 'ccl_reported_failure_or_warning': reported,
        'ccl_public_warnings': native.get('warnings', []), 'ccl_native_replay_status': replay.get('status'),
        'cocoa_reference_valid': bool(ar and ar['reference_valid']),
        'ccl_fractional_error': offset(native.get('value')),
        'cocoa_fractional_error': ar['default_fractional_change'] if ar else None,
        'common_CCL_CQUAD_reference_value': refvalue,
        'cocoa_default_offset_from_CCL_CQUAD': offset(ar['values'][0]) if ar else None,
        'cocoa_refined_offset_from_CCL_CQUAD': offset(ar['values'][-1]) if ar else None,
        'ccl_original_reference_checks': cr, 'ccl_reference_checks': accepted, 'cocoa_reference_checks': ar,
        'ccl_resource_seconds': ccl.get('elapsed_seconds_for_resource_planning'),
        'cocoa_resource_seconds': cocoa.get('elapsed_seconds_for_resource_planning')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    manifest, config = read_manifest(output)
    rows, pair_rows = [], []
    for case in manifest['cases']:
        ccl = load(output / 'ccl' / f'{case["id"]}.json')
        cocoa = load(output / 'cocoa' / f'{case["id"]}.json')
        if ccl is None and cocoa is None:
            continue
        supplement = load(output / 'reference_supplements' / f'{case["id"]}.json')
        if supplement:
            assert supplement['original_result_sha256'] == digest(output / 'ccl' / f'{case["id"]}.json')
        support = load(output / 'support_references' / f'{case["id"]}.json')
        if support:
            assert support['original_result_sha256'] == digest(output / 'ccl' / f'{case["id"]}.json')
        primary = config['primary_pair_zero_based']
        rows.append(comparison_row(case, primary, ccl, cocoa, supplement, True, support))
        # Secondary pairs selected by QAG/spline disagreement are diagnostics,
        # not an unbiased failure-frequency sample.
        pairs = {tuple(primary)} | {tuple(r['pair']) for r in (ccl or {}).get('references', [])}
        for pair in sorted(pairs):
            pair_rows.append(comparison_row(case, list(pair), ccl, cocoa, supplement, list(pair) == primary, support))
    frequencies = []
    for family in manifest['families']:
        planned = sum(c['arm'] == 'uniform' and c['family'] == family for c in manifest['cases'])
        if not planned:
            continue
        selected = [r for r in rows if r['case']['arm'] == 'uniform' and r['case']['family'] == family]
        for code in ('ccl', 'cocoa'):
            valid = [r for r in selected if r[f'{code}_reference_valid'] and r[f'{code}_fractional_error'] is not None]
            for threshold in config['failure_fractional_thresholds']:
                count = sum(abs(r[f'{code}_fractional_error']) > threshold and (code != 'ccl' or r['ccl_silent_success']) for r in valid)
                complete = len(valid) == planned
                frequencies.append({'family': family, 'code': code, 'fractional_error_threshold': threshold,
                    'planned_uniform_cases': planned, 'available_cases': len(selected),
                    'valid_reference_cases': len(valid), 'unresolved_or_unrun_cases': planned - len(valid),
                    'failures': count, 'fraction_among_valid': count / len(valid) if valid else None,
                    'wilson_95_interval_complete_design': wilson(count, planned) if complete else None,
                    'failure_definition': 'silent successful native/replay QAG above error threshold' if code == 'ccl' else 'default fixed-rule spectrum above own-reference error threshold',
                    'reported_failures_or_warnings': sum(r['ccl_reported_failure_or_warning'] for r in selected) if code == 'ccl' else None,
                    'silent_successes_with_valid_reference': sum(r['ccl_silent_success'] for r in valid) if code == 'ccl' else None,
                    'possible_fraction_over_complete_design': [count / planned, (count + planned - len(valid)) / planned],
                    'complete': complete})
    result = {'config_sha256': manifest['config_sha256'], 'rows': rows, 'pair_rows': pair_rows,
        'frequencies': frequencies,
        'interpretation': 'Conditional on the predeclared parameter ranges, uniform measure, fixed input families, ell and primary pair. Targeted points and selected secondary pairs are excluded. Partial results are not a completed failure-rate estimate. Wilson intervals are given only after all planned references are valid.'}
    write_json(output / 'summary.json', result)
    print(json.dumps({'cases_with_any_output': len(rows), 'pair_rows': len(pair_rows), 'frequencies': frequencies}, indent=2))


if __name__ == '__main__':
    main()
