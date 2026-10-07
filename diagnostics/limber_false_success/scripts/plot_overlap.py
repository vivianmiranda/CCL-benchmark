"""Make a shareable figure from the saved narrow-overlap investigation."""
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from paths import plot_paths, INPUTS
root, figures = plot_paths()
frozen = np.load(INPUTS/'limber_outlier_tomography.npz')
sample = np.load(root/'trigger_nodes.npz')
cocoa = np.load(root/'cocoa_case1_level0_nodes.npz')
original = np.load(root/'cocoa_case1_level0.npz')
assert np.array_equal(cocoa['covariance_spectra'], original['covariance_spectra'])
assert np.array_equal(cocoa['data_vector_autos'], original['data_vector_autos'])
# The returned covariance geometry contains chi in c/H0, so convert to Mpc
# before applying the same Limber coordinate used for the CCL illustration.
chi_mpc = cocoa['geometry'][1] * 2997.92458 / .6727
cocoa_logk = np.log(356.155882/chi_mpc)
x, y = sample['x'], sample['y']
kn, ky, gn = sample['kronrod_x'], sample['kronrod_y'], sample['gauss_x']
shared = np.array([min(abs(gn-q))<1.e-12 for q in kn])
half = y >= y.max()/2
left, right = x[half][0], x[half][-1]
counts = {name:int(np.sum((nodes>=left)&(nodes<=right))) for name,nodes in
          [('kronrod',kn),('gauss',gn),('cocoa',cocoa_logk)]}

plt.rcParams.update({'font.size':12,'axes.titlesize':13})
fig, axes = plt.subplots(1,2,figsize=(13,6.1))
fig.subplots_adjust(left=.07,right=.985,bottom=.23,top=.79,wspace=.23)
fig.suptitle('A narrow overlap can fool an adaptive error estimate',
             y=.97,fontsize=19,weight='bold')
fig.text(.5,.902,
         r'Unmodified CCL: 41-point and 20-point estimates agree, but both are 4.0225% low',
         ha='center',fontsize=12.5)

ax = axes[0]
z = frozen['z']
for index,color in [(3,'#3067ad'),(4,'#d36e23')]:
    nz=frozen[f'nz_{index}']
    ax.plot(z,nz/nz.max(),lw=2.3,color=color,label=f'Lens bin {index+1}')
overlap=frozen['nz_3']*frozen['nz_4']
overlap/=overlap.max()
ax.fill_between(z,0,overlap,color='#765298',alpha=.20)
ax.plot(z,overlap,color='#765298',lw=2,ls='--',label='Bin overlap: n₄(z)n₅(z)')
ax.set(xlim=(.45,1.45),ylim=(0,1.4),xlabel='Redshift z',
       ylabel='Each curve / its own maximum',title='Two broad bins make a narrower product')
ax.legend(frameon=False,fontsize=10.5,loc='upper left')
ax.annotate('Overlap near z ≃ 1',xy=(1.00,.86),xytext=(1.13,1.16),fontsize=10,
            ha='center',arrowprops=dict(arrowstyle='->',color='#765298',lw=1.4),color='#5f3d80')

ax=axes[1]
ax.axvspan(left,right,color='#765298',alpha=.08)
ax.plot(x,y/y.max(),color='#27394b',lw=2.3,label='CCL Limber integrand')
ax.scatter(kn[shared],ky[shared]/y.max(),color='#179b78',s=54,zorder=4,
           label='Gauss nodes (shared)')
ax.scatter(kn[~shared],ky[~shared]/y.max(),color='#e56b2f',s=54,zorder=4,
           label='Additional Kronrod nodes')
ax.plot(cocoa_logk,np.full_like(cocoa_logk,-.12),'|',ms=11,mew=1.0,
        color='#3067ad',label='CoCoA all-pairs spectra (covariance module)')
ax.set(xlim=(-2.48,-2.035),ylim=(-.23,1.47),
       xlabel=r'$\ln[k/({\rm Mpc}^{-1})]$',ylabel='Integrand / its maximum',
       title='The first CCL pass misses the top of the peak')
ax.legend(frameon=False,fontsize=9.3,loc='upper right')
for ax in axes:
    ax.spines[['top','right']].set_visible(False)
    ax.grid(alpha=.15)
fig.text(.07,.13,
    'CCL 3.3.3 / GSL 2.7 · ℓ = 355.655882 · w₀ = −1.03892727 · bins 4 × 5 (one-based) · frozen PR #1313 inputs',
    fontsize=10.2)
fig.text(.07,.085,
    'CoCoA covariance-module spectra: 7 panels × 96 nodes in a; default / refined − 1 = −0.0016%. No covariance matrix computed.',
    fontsize=10.2)
fig.text(.07,.04,
    f'Inside the shaded half-maximum width: {counts["kronrod"]} Kronrod nodes ({counts["gauss"]} shared Gauss) versus {counts["cocoa"]} configured CoCoA nodes.',
    fontsize=10.2)
fig.savefig(figures/'narrow_overlap_money.png',dpi=220)
fig.savefig(figures/'narrow_overlap_money.pdf')
plt.close(fig)
(figures/'plot_node_counts.json').write_text(json.dumps(dict(
    half_maximum_logk=[float(left),float(right)],nodes_in_width=counts,
    cocoa_spectra_reproduced_bitwise=True),indent=2)+'\n')
print(counts)
