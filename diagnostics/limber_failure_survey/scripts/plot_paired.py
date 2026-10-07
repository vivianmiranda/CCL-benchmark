"""Diagnostic triangles and overlap plots showing measured CCL and CoCoA."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, MaxNLocator, FuncFormatter

COLORS = ('#c95a24', '#246fa3')
LABELS = ('CCL QAG', 'CoCoA 96-node spectra')


def valid(row):
    return (row['ccl_reference_valid'] and row['cocoa_reference_valid']
            and row.get('cocoa_default_offset_from_CCL_CQUAD') is not None)


def failure(row, threshold):
    return valid(row) and ((row['ccl_silent_success'] and abs(row['ccl_fractional_error']) > threshold)
                          or abs(row['cocoa_fractional_error']) > threshold)


def save(fig, folder, name):
    fig.savefig(folder / f'{name}.png', dpi=190)
    fig.savefig(folder / f'{name}.pdf')
    plt.close(fig)


def own_table(ax, measured):
    """Explicitly separate quadrature convergence and inter-code offsets."""
    ax.axis('off')
    own_ccl = max(abs(r['ccl_fractional_error']) for r in measured) * 100
    own_cocoa = max(abs(r['cocoa_fractional_error']) for r in measured) * 100
    common_cocoa = max(abs(r['cocoa_default_offset_from_CCL_CQUAD']) for r in measured) * 100
    ax.text(0, .99, 'Different questions, different references', weight='bold', fontsize=12, va='top')
    table = ax.table(cellText=[['CCL QAG', f'{own_ccl:.5g}%', f'{own_ccl:.5g}%'],
                             ['CoCoA', f'{own_cocoa:.5g}%', f'{common_cocoa:.5g}%']],
                     colLabels=['Code', 'Own quadrature', 'Common reference'],
                     cellLoc='center', colWidths=[.23, .34, .43], bbox=[0, .43, 1, .42])
    table.auto_set_font_size(False); table.set_fontsize(10)
    ax.text(0, .34, 'Largest absolute deviations among displayed cases.\n'
            'Own quadrature: QAG vs CQUAD; CoCoA 96 vs 1024.\n'
            'Common reference: both codes vs CCL CQUAD.\n'
            'A cross-code offset is not itself quadrature error.', va='top', fontsize=10, linespacing=1.5)


def triangle(points, family, figures, threshold):
    measured = [r for r in points if valid(r)]
    if not measured:
        return
    coordinates = {tuple(r['case']['parameters'][k] for k in ('w0', 'wa', 'Omega_m')) for r in points}
    if not any(r['case']['arm'] == 'uniform' for r in points) and len(coordinates) < 3:
        return  # A tuned input at one cosmology is not a cosmology triangle.
    residual = np.array([[100*r['ccl_fractional_error'], 100*r['cocoa_default_offset_from_CCL_CQUAD']]
                         for r in measured])
    norm = Normalize(0, max(.1, float(np.max(np.abs(residual)))))
    cmap = plt.get_cmap('viridis')
    keys = ['w0', 'wa', 'Omega_m']; labels = [r'$w_0$', r'$w_a$', r'$\Omega_m$']
    fig, axes = plt.subplots(3, 3, figsize=(14.3, 11.6))
    fig.subplots_adjust(left=.08, right=.86, bottom=.18, top=.83, wspace=.35, hspace=.40)
    for i in range(3):
        for j in range(3):
            ax = axes[i, j]
            if j > i:
                ax.set_visible(False); continue
            if i == j:
                for row, (cc, co) in zip(measured, residual):
                    x = row['case']['parameters'][keys[j]]
                    ax.plot([x, x], [cc, co], color='0.75', lw=.7, zorder=1)
                    for y, color, marker in [(cc, COLORS[0], 'o'), (co, COLORS[1], 's')]:
                        ax.scatter(x, y, s=60, color=color, marker=marker, edgecolors='white', linewidths=.5, zorder=3)
                    for value, failed in [(cc, row['ccl_silent_success'] and abs(row['ccl_fractional_error']) > threshold),
                                          (co, abs(row['cocoa_fractional_error']) > threshold)]:
                        if failed:
                            ax.scatter(x, value, s=118, facecolors='none', edgecolors='#bc272d', linewidths=1.3, zorder=4)
                ax.axhline(0, color='0.4', lw=.7)
                ax.set_ylabel('Common-reference offset [%]', fontsize=10)
                if np.max(np.abs(residual)) > 1:
                    ax.set_yscale('symlog', linthresh=.08)
                    ax.yaxis.set_major_locator(FixedLocator([-100, -10, -1, -.1, 0, .1, 1, 10, 100]))
                    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:g}'))
                    ax.set_ylim(min(-.1, float(residual.min())*1.8), max(.12, float(residual.max())*1.8))
            else:
                for row in points:
                    p = row['case']['parameters']; x, y = p[keys[j]], p[keys[i]]
                    if not valid(row):
                        ax.plot(x, y, marker='x', color='0.7', ms=8, linestyle='None'); continue
                    cc = abs(100*row['ccl_fractional_error'])
                    co = abs(100*row['cocoa_default_offset_from_CCL_CQUAD'])
                    targeted = row['case']['arm'] != 'uniform'
                    ax.plot(x, y, marker='D' if targeted else 'o', ms=12, fillstyle='left',
                            markerfacecolor=cmap(norm(cc)), markerfacecoloralt=cmap(norm(co)),
                            markeredgecolor='#bc272d' if failure(row, threshold) else 'black',
                            markeredgewidth=2 if failure(row, threshold) else .6, linestyle='None')
                ax.set_ylabel(labels[i])
            ax.set_xlabel(labels[j]); ax.grid(alpha=.17); ax.tick_params(axis='both', labelsize=10)
            ax.xaxis.set_major_locator(MaxNLocator(4))
    table_ax = fig.add_axes([.47, .615, .38, .21])
    own_table(table_ax, measured)
    bar = fig.add_axes([.903, .285, .018, .31])
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=bar, label='Absolute offset from CCL CQUAD [%]')
    fig.suptitle(f'Paired quadrature checks: {family}', fontsize=21, y=.975)
    fig.text(.5, .934, 'Common reference: the same spectrum evaluated with CCL kernels and GSL CQUAD', ha='center', fontsize=12)
    fig.text(.5, .904, 'Split circles: left = CCL QAG; right = CoCoA. Diamonds = deliberate search/control.', ha='center', fontsize=11)
    handles = [Line2D([], [], ls='', marker=marker, color=color, label=label)
               for label, color, marker in zip(LABELS, COLORS, ['o', 's'])]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.48, .887), ncol=2, frameon=False)
    uniform = sum(r['case']['arm'] == 'uniform' for r in points)
    highlighted = sum(failure(r, threshold) for r in points)
    fig.text(.08, .127, f'Red outline: either code exceeds {100*threshold:g}% own-reference quadrature error ({highlighted} displayed cases).', fontsize=11)
    fig.text(.08, .096, f'{uniform} uniform cases available; {len(points)-uniform} targeted/control cases. Diagnostic coordinates, not posterior samples.', fontsize=10)
    fig.text(.08, .064, 'CoCoA: actual covariance-module all-pairs spectra. No covariance matrix; no ordinary data-vector cross-bin test.', fontsize=10)
    fig.text(.08, .033, 'Grey crosses: unresolved or not yet paired. Targeted points are excluded from frequency estimates. No universal immunity claim.', fontsize=10)
    save(fig, figures, f'{family}_triangle')


def overlap_plot(points, family, pair, input_path, figures, threshold):
    measured = [r for r in points if valid(r)]
    if not measured:
        return
    data = np.load(input_path); z, nz = data['z'], data['nz']
    i, j = pair; product = nz[i] * nz[j]
    fig, axes = plt.subplots(1, 3, figsize=(17.3, 6.0))
    fig.subplots_adjust(left=.057, right=.99, top=.76, bottom=.25, wspace=.29)
    for k, color in [(i, COLORS[0]), (j, COLORS[1])]:
        axes[0].plot(z, nz[k]/nz[k].max(), lw=1.8, label=f'Bin {k+1}', color=color)
    axes[0].plot(z, product/product.max(), lw=2.5, color='#78518e', label=f'Product: bins {i+1} × {j+1}')
    support = z[product > product.max() * .001]
    if len(support):
        width = max(.06, .15 * (support[-1]-support[0]))
        axes[0].set_xlim(max(z[0], support[0]-width), min(z[-1], support[-1]+width))
    axes[0].set(xlabel='Redshift z', ylabel='Each curve / its maximum', title='Shared input and its overlap')
    axes[0].legend(frameon=False, fontsize=10)
    x = np.arange(len(measured))
    common = np.array([[100*r['ccl_fractional_error'], 100*r['cocoa_default_offset_from_CCL_CQUAD']] for r in measured])
    own = np.array([[100*abs(r['ccl_fractional_error']), 100*abs(r['cocoa_fractional_error'])] for r in measured])
    for k, (label, color, marker) in enumerate(zip(LABELS, COLORS, ['o', 's'])):
        axes[1].scatter(x, common[:, k], label=label, color=color, marker=marker, s=55, alpha=.85)
        axes[2].scatter(x, np.maximum(own[:, k], 1e-10), label=label, color=color, marker=marker, s=55, alpha=.85)
    axes[1].axhline(0, color='0.4', lw=.7)
    if np.max(np.abs(common)) > 1:
        axes[1].set_yscale('symlog', linthresh=.08)
        axes[1].yaxis.set_major_locator(FixedLocator([-100, -10, -1, -.1, 0, .1, 1, 10, 100]))
        axes[1].yaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:g}'))
    axes[1].set(xlabel='Paired case index', ylabel='Signed offset [%]', title='Both codes vs CCL CQUAD')
    axes[2].axhline(100*threshold, color='#bc272d', lw=1.1, ls='--', label=f'{100*threshold:g}% threshold')
    axes[2].set(xlabel='Paired case index', ylabel='Absolute quadrature error [%]', yscale='log', title='Each code vs its own refinement')
    for ax in axes[1:]:
        ax.legend(frameon=False, fontsize=9)
        ax.set_xlim(-.6, len(measured)-.4)
        ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
        if len(measured) <= 3:
            ax.set_xticks(np.arange(len(measured)))
    for ax in axes:
        ax.grid(alpha=.17)
    fig.suptitle(f'{family} · lens bins {i+1} × {j+1}: overlap and actual spectra', fontsize=19, y=.97)
    fig.text(.5, .875, 'The redshift product is identical across cosmologies; changing cosmology changes its distance mapping.', ha='center', fontsize=12)
    fig.text(.5, .825, 'Common-reference offsets include kernel/interpolation differences. Own-code convergence isolates quadrature resolution.', ha='center', fontsize=11)
    fig.text(.057, .14, f'{len(measured)} paired cases; {len(points)-len(measured)} unresolved. Selected secondary pairs and targeted points are diagnostic, not incidence samples.', fontsize=10)
    fig.text(.057, .097, 'CCL reference: successful checked CQUAD refinements. CoCoA reference: 1024 nodes/panel, checked against 512.', fontsize=10)
    fig.text(.057, .054, 'CoCoA points are actual covariance-module spectra. Any CoCoA underresolution is retained. Redshift curves are drawn once per input/pair.', fontsize=10)
    save(fig, figures, f'{family}_pair_{i+1}_{j+1}_overlap_and_errors')


def overlap_atlas(groups, figures, threshold):
    """Four input/pair rows per page, each showing both measured codes."""
    groups = [g for g in groups if any(valid(r) for r in g['points'])]
    for page, first in enumerate(range(0, len(groups), 4), 1):
        selected = groups[first:first+4]
        fig, axes = plt.subplots(len(selected), 3, figsize=(16.5, 3.25*len(selected)+2.0), squeeze=False)
        fig.subplots_adjust(left=.07, right=.98, top=1-1.05/(3.25*len(selected)+2),
                            bottom=1.35/(3.25*len(selected)+2), hspace=.7, wspace=.33)
        for row_axes, group in zip(axes, selected):
            family, pair = group['family'], group['pair']
            measured = [r for r in group['points'] if valid(r)]
            data = np.load(group['input_path']); z, nz = data['z'], data['nz']
            i, j = pair; product = nz[i]*nz[j]
            for k, color in [(i, COLORS[0]), (j, COLORS[1])]:
                row_axes[0].plot(z, nz[k]/nz[k].max(), color=color, lw=1.1, alpha=.65)
            row_axes[0].plot(z, product/product.max(), color='#78518e', lw=2)
            support = z[product > product.max()*.001]
            if len(support):
                margin = max(.06, .1*(support[-1]-support[0]))
                row_axes[0].set_xlim(max(z[0], support[0]-margin), min(z[-1], support[-1]+margin))
            name = family.replace('shared_narrow_peak', 'Shared peak, z=1, σ=0.01').replace('pr_original', 'Original PR distributions')
            row_axes[0].set(title=f'{name}\nLens bins {i+1} × {j+1}', xlabel='Redshift z', ylabel='Input / its peak')
            x = np.arange(len(measured))
            common = np.array([[r['ccl_fractional_error'], r['cocoa_default_offset_from_CCL_CQUAD']] for r in measured])*100
            own = np.abs(np.array([[r['ccl_fractional_error'], r['cocoa_fractional_error']] for r in measured]))*100
            for k, (color, marker) in enumerate(zip(COLORS, ['o', 's'])):
                row_axes[1].scatter(x, common[:, k], c=color, marker=marker, s=45, alpha=.8)
                row_axes[2].scatter(x, np.maximum(own[:, k], 1e-10), c=color, marker=marker, s=45, alpha=.8)
            row_axes[1].axhline(0, color='0.4', lw=.6)
            if np.max(np.abs(common)) > 1:
                row_axes[1].set_yscale('symlog', linthresh=.08)
                row_axes[1].yaxis.set_major_locator(FixedLocator([-100, -10, -1, -.1, 0, .1, 1, 10, 100]))
                row_axes[1].yaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:g}'))
            row_axes[1].set(title='Both codes vs CCL CQUAD', xlabel=f'Paired case index ({len(measured)} cases)', ylabel='Signed offset [%]')
            row_axes[2].axhline(100*threshold, color='#bc272d', lw=.9, ls='--')
            row_axes[2].set(title='Each code vs own refinement', xlabel='Same paired case index', ylabel='Quadrature error [%]', yscale='log')
            for ax in row_axes[1:]:
                ax.set_xlim(-.6, len(measured)-.4); ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
                if len(measured) <= 3:
                    ax.set_xticks(np.arange(len(measured)))
            for ax in row_axes:
                ax.grid(alpha=.17); ax.tick_params(labelsize=10); ax.title.set_fontsize(11)
        fig.suptitle(f'Distinct redshift overlaps and actual numerical checks · page {page}', fontsize=18, y=.98)
        handles = [Line2D([], [], ls='', marker=m, color=c, label=l) for l,c,m in zip(LABELS,COLORS,['o','s'])]
        handles.append(Line2D([], [], color='#78518e', label='Normalized nᵢ(z)nⱼ(z)'))
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, 1-.42/(3.25*len(selected)+2)), ncol=3, frameon=False, fontsize=11)
        fig.text(.07, .34/(3.25*len(selected)+2), f'Each distinct input/pair appears once. Red dashed line: {100*threshold:g}%. Common-reference offsets include kernel/interpolation differences.', fontsize=10)
        fig.text(.07, .10/(3.25*len(selected)+2), 'CoCoA: actual covariance-module spectra, with 96 vs 1024-node convergence. Selected pairs and targeted cases are not incidence samples.', fontsize=10)
        save(fig, figures, f'overlap_error_atlas_{page:02d}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--targeted-summary', type=Path, action='append', default=[])
    parser.add_argument('--failure-threshold', type=float, default=.001, help='Fractional own-reference error highlighted in triangles (default 0.1%%)')
    parser.add_argument('--product-threshold', type=float, default=.0001, help='Plot every input/pair with any validated error above this fraction (default 0.01%%)')
    args = parser.parse_args(); folder = args.output.resolve()
    rows, pair_rows, input_paths = [], [], {}
    for path in [folder/'summary.json'] + args.targeted_summary:
        result = json.loads(path.read_text())
        rows.extend(result['rows']); pair_rows.extend(result.get('pair_rows', result['rows']))
        for row in result['rows']:
            input_paths[row['case']['family']] = path.parent/'inputs'/f'{row["case"]["family"]}.npz'
    figures = folder/'figures'; figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 12, 'axes.titlesize': 13, 'axes.labelsize': 12})
    outputs, groups = [], []
    for family in sorted({r['case']['family'] for r in rows}):
        triangle([r for r in rows if r['case']['family'] == family], family, figures, args.failure_threshold)
        selected = [r for r in pair_rows if r['case']['family'] == family]
        pairs = {tuple(r['pair']) for r in selected if failure(r, args.product_threshold)}
        # Always show the predeclared pair, even if no failure is found.
        pairs |= {tuple(r.get('pair', [3, 4])) for r in rows if r['case']['family'] == family}
        for pair in sorted(pairs):
            points = [r for r in selected if tuple(r['pair']) == pair]
            overlap_plot(points, family, pair, input_paths[family], figures, args.product_threshold)
            groups.append({'family': family, 'pair': pair, 'input_path': input_paths[family], 'points': points})
            outputs.append({'family': family, 'pair': list(pair), 'paired_cases': sum(valid(r) for r in points),
                            'case_ids_in_plot': [r['case']['id'] for r in points if valid(r)]})
    overlap_atlas(groups, figures, args.product_threshold)
    (figures/'plot_index.json').write_text(json.dumps({'overlap_figures': outputs,
        'source_summaries': [{'path': str(path.resolve()), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                             for path in [folder/'summary.json']+args.targeted_summary],
        'triangle_highlight_fraction': args.failure_threshold, 'overlap_selection_fraction': args.product_threshold}, indent=2)+'\n')


if __name__ == '__main__':
    main()
