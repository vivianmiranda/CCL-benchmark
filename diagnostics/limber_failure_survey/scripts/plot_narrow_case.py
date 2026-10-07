"""Plot saved narrow-input measurements and actual quadrature locations."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from common import digest, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--case', required=True)
    args = parser.parse_args()
    root = args.output
    summary = json.loads((root/'summary.json').read_text())
    row = next(r for r in summary['rows'] if r['case']['id'] == args.case)
    assert row['ccl_reference_valid'] and row['cocoa_reference_valid']
    assert row['ccl_silent_success']
    family = row['case']['family']
    manifest = json.loads((root/'manifest.json').read_text())
    definition = manifest['families'][family]['definition']
    files = {'summary': root/'summary.json', 'inputs': root/'inputs'/f'{family}.npz',
             'integrand': root/'ccl'/f'{args.case}_integrand.npz',
             'cocoa': root/'cocoa'/f'{args.case}.npz'}
    inp, sample, cocoa = [np.load(files[key]) for key in ('inputs', 'integrand', 'cocoa')]
    z, nz = inp['z'], inp['nz']
    i, j = row['pair']
    product = nz[i]*nz[j]
    x, y = sample['logk'], sample['integrand']
    kn, ky, gn = sample['kronrod_logk'], sample['kronrod_integrand'], sample['gauss_logk']
    shared = np.array([np.min(abs(gn-v)) < 1e-12 for v in kn])
    ell = float(cocoa['ell'][0]); h = row['case']['parameters']['h']
    xc = np.log((ell+.5)/(cocoa['geometry_96'][1]*2997.92458/h))
    peak = x[np.argmax(y)]
    half = x[y >= y.max()/2]
    colors = dict(blue='#236ca4', red='#bd392f', purple='#765298', green='#14916b', orange='#db7b20')
    plt.rcParams.update({'font.size': 13, 'axes.titlesize': 15, 'axes.labelsize': 14})
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.4))
    fig.subplots_adjust(left=.075, right=.98, top=.74, bottom=.30, wspace=.26)
    fig.suptitle('A narrow overlap exposes false convergence', y=.97, fontsize=23, weight='bold')
    fig.text(.5, .895, f'Deliberate stress test: 15% shared population at z = {definition["shared_mean_z"]:.5f}, σz = {definition["shared_sigma_z"]:g}', ha='center', fontsize=15)
    p = row['case']['parameters']
    fig.text(.5, .845, f'w₀ = {p["w0"]:.6f}, wₐ = {p["wa"]:g}, Ωₘ = {p["Omega_m"]:g} · lens bins {i+1} × {j+1} · ℓ = {ell:.3f}', ha='center', fontsize=13)
    ax = axes[0]
    for k, color in [(i, colors['blue']), (j, colors['orange'])]:
        ax.plot(z, nz[k]/nz[k].max(), color=color, lw=1.8, label=f'Lens bin {k+1}')
    ax.fill_between(z, 0, product/product.max(), color=colors['purple'], alpha=.12)
    ax.plot(z, product/product.max(), color=colors['purple'], ls='--', lw=2, label='Product n₄(z)n₅(z)')
    ax.axvline(1, color='.5', ls=':', lw=1)
    ax.text(1.006, .37, 'CoCoA panel\nboundary at z = 1', fontsize=10, color='.35')
    ax.set(xlim=(.83, 1.13), ylim=(0, 1.4), xlabel='Redshift z', ylabel='Each curve / its maximum', title='The common population makes a sharp peak')
    ax.legend(loc='upper right', frameon=False, fontsize=11)
    ax = axes[1]
    ax.axvspan(half[0], half[-1], color=colors['purple'], alpha=.1)
    ax.plot(x, y/y.max(), color='#293f50', lw=2.1, label='CCL Limber integrand')
    ax.scatter(kn[shared], ky[shared]/y.max(), color=colors['green'], s=46, zorder=4, label='Gauss nodes (20)')
    ax.scatter(kn[~shared], ky[~shared]/y.max(), color=colors['orange'], s=46, zorder=4, label='Added Kronrod nodes (21)')
    ax.plot(xc, np.full(len(xc), -.09), '|', ms=14, color=colors['blue'], label='CoCoA fixed nodes (96 per panel)')
    ax.set(xlim=(peak-.09, peak+.09), ylim=(-.20, 1.48), xlabel=r'$\ln[k/({\rm Mpc}^{-1})]$', ylabel='Integrand / its maximum', title='The initial QAG pass misses the peak')
    ax.legend(loc='upper right', frameon=False, fontsize=10)
    for ax in axes:
        ax.grid(alpha=.12); ax.spines[['top', 'right']].set_visible(False)
    qerr = 100*row['ccl_fractional_error']
    aerr = 100*row['cocoa_fractional_error']
    fig.text(.075, .205, f'CCL QAG: {qerr:.2f}% error; GSL_SUCCESS after 41 evaluations.', fontsize=16, color=colors['red'], weight='bold')
    ladder_path = root/'cocoa_integration_ladder'/f'{args.case}.json'
    if ladder_path.exists():
        ladder = json.loads(ladder_path.read_text())
        lr = next(r for r in ladder['pairs'] if r['pair'] == row['pair'])
        e2, e3 = [100*(lr['values'][k]/lr['values'][-1]-1) for k in (2,3)]
        files['integration_ladder'] = ladder_path
        fig.text(.075, .148, f'CoCoA integration_accuracy: 0 → {aerr:.3f}%; 2 → {e2:.5f}%; 3 → {e3:.6f}% (vs level 4).', fontsize=14, color=colors['blue'])
    else:
        fig.text(.075, .148, f'CoCoA: {aerr:.3f}% default-versus-refined error; its fixed rule also needs more nodes here.', fontsize=14, color=colors['blue'])
    fig.text(.075, .093, 'Reference: refined CCL CQUAD, checked with subdivisions around the peak and twice as many input samples.', fontsize=11)
    fig.text(.075, .052, f'Common-reference CoCoA offsets: default {100*row["cocoa_default_offset_from_CCL_CQUAD"]:+.3f}%; refined {100*row["cocoa_refined_offset_from_CCL_CQUAD"]:+.4f}% (includes interpolation differences).', fontsize=11)
    fig.text(.075, .014, 'CoCoA covariance-module angular spectra, not ordinary data-vector integration. Targeted example, not a failure-frequency estimate.', fontsize=10)
    folder = root/'figures'; folder.mkdir(exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(folder/f'narrow_overlap_sampling.{ext}', dpi=190)
    plt.close(fig)
    write_json(folder/'narrow_overlap_sampling.json', {'case': args.case, 'source_sha256': {k: digest(v) for k,v in files.items()},
        'script_sha256': digest(__file__), 'ccl_error_percent': qerr, 'cocoa_own_error_percent': aerr,
        'peak_logk': float(peak), 'FWHM_logk': half.tolist(),
        'kronrod_nodes_in_FWHM': int(np.sum((kn>=half[0]) & (kn<=half[-1]))),
        'cocoa_nodes_in_FWHM': int(np.sum((xc>=half[0]) & (xc<=half[-1])))})


if __name__ == '__main__':
    main()
