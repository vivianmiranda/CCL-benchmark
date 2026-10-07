"""Explain the fixed 96-point rule using saved, actual CoCoA node positions.

This is a figure-only extension of the completed CCL diagnostic. It does
not rerun either cosmology code or modify any numerical results.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

from paths import plot_paths, INPUTS
root, figures = plot_paths()
frozen = np.load(INPUTS/'limber_outlier_tomography.npz')
sample = np.load(root/'trigger_nodes.npz')
cocoa = np.load(root/'cocoa_case1_level0_nodes.npz')
tables = np.load(root/'plot_geometry.npz')
a = cocoa['geometry'][0]
logk_cocoa = np.log(356.155882/(cocoa['geometry'][1]*2997.92458/.6727))
x, y = sample['x'], sample['y']
kn, ky, gn = sample['kronrod_x'], sample['kronrod_y'], sample['gauss_x']
shared = np.array([min(abs(gn-q))<1.e-12 for q in kn])
half = y >= y.max()/2
half_left, half_right = x[half][0], x[half][-1]
boundary_logk = np.log(356.155882/(np.interp(1.,tables['z_1D'],tables['chi'])/.6727))
panel_edges = [.4, .5, 1./1.7]
panel_masks = [(a>left)&(a<right) for left,right in zip(panel_edges[:-1],panel_edges[1:])]
canonical_nodes, canonical_weights = np.polynomial.legendre.leggauss(96)
fraction = (canonical_nodes+1)/2
normalized_weights = canonical_weights/2
for mask,left,right in zip(panel_masks,panel_edges[:-1],panel_edges[1:]):
    assert mask.sum() == 96
    np.testing.assert_allclose(np.sort(a[mask]), left+fraction*(right-left),
                               rtol=0,atol=2.e-16)

blue, teal = '#3067ad', '#237f91'
purple, green, orange = '#765298', '#179b78', '#e56b2f'
plt.rcParams.update({'font.size':11,'axes.titlesize':12.5})
fig, axes = plt.subplots(2,2,figsize=(13.4,9.8),
                         gridspec_kw={'height_ratios':[1.35,1.]})
fig.subplots_adjust(left=.075,right=.98,bottom=.17,top=.86,
                    wspace=.25,hspace=.44)
fig.suptitle('Why fixed quadrature points cluster near this peak',
             y=.977,fontsize=20,weight='bold')
fig.text(.5,.936,'The points are prescribed in advance; Gauss–Legendre nodes crowd the ends of every panel.',
         ha='center',fontsize=12.5)
fig.text(.5,.905,'Here, a predefined panel boundary at z = 1 lies close to the integrand peak at z ≃ 0.994.',
         ha='center',fontsize=12)

# A: the fixed input distributions, with the physical panel boundary marked.
ax=axes[0,0]
z=frozen['z']
for i,color in [(3,blue),(4,'#d36e23')]:
    nz=frozen[f'nz_{i}']
    ax.plot(z,nz/nz.max(),color=color,lw=2.2,label=f'Lens bin {i+1}')
overlap=frozen['nz_3']*frozen['nz_4']
overlap/=overlap.max()
ax.fill_between(z,0,overlap,color=purple,alpha=.15)
ax.plot(z,overlap,color=purple,ls='--',lw=2,label='Overlap: n₄(z)n₅(z)')
ax.axvline(1,color='0.45',ls=':',lw=1.4)
ax.text(1.016,1.13,'Panel boundary\nz = 1',fontsize=10,color='0.3')
ax.set(xlim=(.45,1.45),ylim=(0,1.38),xlabel='Redshift z',
       ylabel='Each curve / its maximum',title='A  The bin overlap makes a narrow peak')
ax.legend(frameon=False,fontsize=10,loc='upper left')

# B: actual CCL sampling and actual CoCoA positions, coloured by a-panel.
ax=axes[0,1]
ax.axvspan(half_left,half_right,color=purple,alpha=.08)
ax.plot(x,y/y.max(),lw=2.2,color='#27394b',label='CCL Limber integrand')
ax.scatter(kn[shared],ky[shared]/y.max(),s=42,color=green,zorder=4,label='Shared Gauss nodes')
ax.scatter(kn[~shared],ky[~shared]/y.max(),s=42,color=orange,zorder=4,label='Added Kronrod nodes')
for mask,color in zip(panel_masks,[blue,teal]):
    ax.plot(logk_cocoa[mask],np.full(mask.sum(),-.12),'|',color=color,ms=11,mew=.9)
ax.plot([],[],'|',color=blue,ms=11,
        label='CoCoA all-pairs spectra\n(covariance-module ticks)')
ax.axvline(boundary_logk,color='0.45',ls=':',lw=1.4)
ax.text(boundary_logk+.008,.77,'z = 1',fontsize=9,color='0.3')
ax.set(xlim=(-2.48,-2.035),ylim=(-.24,1.65),
       xlabel=r'$\ln[k/({\rm Mpc}^{-1})]$',ylabel='Integrand / its maximum',
       title='B  Two panels meet close to the peak')
ax.legend(frameon=False,fontsize=9.0,loc='upper left')

# C: exact canonical rule. Unequal weights account for unequal node spacing.
ax=axes[1,0]
ax.plot(fraction,normalized_weights,color=blue,lw=1,alpha=.6)
ax.scatter(fraction,normalized_weights,color=blue,s=10,zorder=3,
           label='96 node–weight pairs')
ax.plot(fraction,np.full(96,-.001),'|',color=blue,ms=9,mew=.9)
ax.axvline(0,color='0.6',ls=':',lw=1)
ax.axvline(1,color='0.6',ls=':',lw=1)
ax.annotate('Denser nodes,\nsmaller weights',xy=(fraction[5],normalized_weights[5]),
            xytext=(.13,.0075),fontsize=10,
            arrowprops=dict(arrowstyle='->',color=blue),color=blue)
ax.annotate('The same at\nthe other end',xy=(fraction[-6],normalized_weights[-6]),
            xytext=(.65,.0075),fontsize=10,
            arrowprops=dict(arrowstyle='->',color=blue),color=blue)
ax.set(xlim=(-.02,1.02),ylim=(-.002, .021),xlabel='Position within one panel (0 = left edge, 1 = right edge)',
       ylabel='Quadrature weight',title='C  All 96 nodes of a single fixed rule')
ax.yaxis.set_major_formatter(PercentFormatter(xmax=1,decimals=1))
ax.set_yticks([0,.005,.01,.015,.02])
ax.legend(frameon=False,fontsize=9.3,loc='upper center')
ax.text(.5,.15,'No node lies exactly on an endpoint.',transform=ax.transAxes,
        ha='center',fontsize=9,color='0.35')

# D: two adjoining panels from the saved run, with every node shown in a.
ax=axes[1,1]
for mask,color,left,right in zip(panel_masks,[blue,teal],panel_edges[:-1],panel_edges[1:]):
    ax.axvspan(left,right,color=color,alpha=.07)
    ax.hlines(.46,left,right,color=color,lw=.8)
    ax.plot(a[mask],np.full(mask.sum(),.46),'|',ms=21,mew=.85,color=color)
    ax.text((left+right)/2,.80,'96 fixed nodes',ha='center',color=color,fontsize=11)
    ax.annotate('',xy=(right,.69),xytext=(left,.69),
                arrowprops=dict(arrowstyle='<->',color=color,lw=1.2))
for edge in panel_edges:
    ax.axvline(edge,color='0.45',ls=':',lw=1.1)
ax.annotate('Nodes crowd both sides\nof the shared boundary',xy=(.5,.49),
            xytext=(.5,.20),ha='center',va='center',fontsize=10,
            arrowprops=dict(arrowstyle='->',color='0.3'))
ax.text(.45,1.03,'z = 1.5 → 1',ha='center',fontsize=10,color=blue)
ax.text((.5+1/1.7)/2,1.03,'z = 1 → 0.7',ha='center',fontsize=10,color=teal)
ax.set(xlim=(.397,.591),ylim=(0,1.16),yticks=[],
       xticks=[.4,.425,.45,.475,.5,.525,.55,.575],
       xlabel=r'Scale factor $a=1/(1+z)$',title='D  Actual nodes in the two neighbouring panels')
ax.tick_params(axis='x',labelsize=9)
ax.spines[['left']].set_visible(False)

for ax in axes.flat:
    ax.spines[['top','right']].set_visible(False)
    ax.grid(alpha=.13)
fig.text(.075,.090,
 'CCL 3.3.3 / GSL 2.7 · ℓ = 355.655882 · w₀ = −1.03892727 · bins 4 × 5 (one-based) · PR #1313 frozen inputs',fontsize=10)
fig.text(.075,.060,
 'CCL initial QAG pass: 4.0225% low. CoCoA covariance-module spectra: default versus refined differs by 0.0016% at this point.',fontsize=10)
fig.text(.075,.030,
 'CoCoA uses 7 panels × 96 nodes here. This is not its ordinary data-vector quadrature; no covariance matrix was computed.',fontsize=10)
fig.savefig(figures/'narrow_overlap_sampling.png',dpi=220)
fig.savefig(figures/'narrow_overlap_sampling.pdf')
plt.close(fig)
print('Saved annotated PNG/PDF; all displayed panel nodes match the saved run to 2e-16 in a.')
