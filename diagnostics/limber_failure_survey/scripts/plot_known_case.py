"""Show the reproduced QAG error and the measured CoCoA spectrum together."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from common import digest, write_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    name = 'known__pr_original'
    files = {k: a.output/f'{k}/{name}.json' for k in ('ccl', 'cocoa', 'reference_supplements')}
    records = {k: json.loads(f.read_text()) for k, f in files.items()}
    ref = next(r for r in records['reference_supplements']['references'] if r['primary'])
    assert ref['reference_valid']
    native = next(r['qag_quad'] for r in records['ccl']['native'] if r['pair'] == [3, 4])
    original = next(r for r in records['ccl']['references'] if r['primary'])
    cocoa = next(r for r in records['cocoa']['pairs'] if r['pair'] == [3, 4])
    assert cocoa['reference_valid'] and not native['warnings']
    assert original['native_qag_replay']['status'] == [0]
    split = original['split4']
    assert not any(split['status'])
    values = np.array([native['value'], cocoa['values'][0], split['value']])
    delta = 100*(values/ref['reference_value']-1)
    labels = ['CCL default QAG', 'CoCoA fixed quadrature', 'CCL split-4 QAG check']
    colors = ['#c44135', '#236ca4', '#438451']
    plt.rcParams.update({'font.size': 14, 'axes.labelsize': 16, 'axes.titlesize': 17})
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.9), gridspec_kw={'width_ratios': [1.2, 1]})
    fig.subplots_adjust(left=.21, right=.975, top=.70, bottom=.28, wspace=.34)
    for ax in axes:
        ax.axvline(0, color='.25', lw=1)
        ax.set_yticks(range(3), labels if ax is axes[0] else ['']*3)
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=.2)
        ax.set_xlabel('Offset from CCL CQUAD [%]')
    axes[0].barh(range(3), delta, color=colors, height=.52)
    axes[0].set_xlim(-4.65, .65)
    axes[0].set_title('Full error scale')
    axes[0].text(-3.9, 0, f'{delta[0]:.4f}%', ha='left', va='center', fontsize=15, color='white')
    axes[1].barh([1, 2], delta[1:], color=colors[1:], height=.52)
    axes[1].set_xlim(-.01, .08)
    axes[1].set_title('Zoom near agreement')
    axes[1].text(.035, 0, 'Default QAG lies off scale', ha='center', color=colors[0], fontsize=12)
    axes[1].text(.058, 1, f'+{delta[1]:.4f}%', va='center', fontsize=12)
    axes[1].text(.003, 2, f'{delta[2]:.2g}%', va='center', fontsize=12)
    for ax in axes:
        ax.set_ylim(2.4, -.4)
    fig.suptitle('Reproduced false success: 4.02% spectrum error', fontsize=23, y=.98)
    fig.text(.5, .89, r'$w_0=-1.0389272676$, $w_a=0$, $\Omega_m=0.3156$; bins 4 × 5, $\ell=355.655882$', ha='center', fontsize=15)
    fig.text(.5, .81, 'CCL returns GSL_SUCCESS after 41 evaluations; the refined references agree.', ha='center', fontsize=14)
    fig.text(.05, .14, f'CoCoA 96 → 1024 nodes per panel changes its own spectrum by {abs(100*cocoa["default_fractional_change"]):.4f}%.', fontsize=13)
    fig.text(.05, .09, 'The common-reference CoCoA offset also includes different kernel and interpolation conventions.', fontsize=12)
    fig.text(.05, .04, 'CoCoA covariance-module spectra, not ordinary data-vector integration. Deliberate control; not a failure frequency.', fontsize=12)
    folder = a.output/'figures'; folder.mkdir(exist_ok=True)
    for suffix in ('png', 'pdf'):
        fig.savefig(folder/f'known_false_success.{suffix}', dpi=190)
    plt.close(fig)
    write_json(folder/'known_false_success.json', {
        'reference': 'Refined CCL CQUAD; successful split-QAG controls verified separately',
        'reference_value': ref['reference_value'], 'labels': labels,
        'common_reference_offsets_percent': delta.tolist(),
        'cocoa_own_refinement_percent': 100*cocoa['default_fractional_change'],
        'source_sha256': {k: digest(f) for k, f in files.items()},
        'script_sha256': digest(__file__),
    })


if __name__ == '__main__':
    main()
