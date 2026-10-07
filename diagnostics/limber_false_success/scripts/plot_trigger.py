"""Plot saved native CCL measurements and the quadrature acceptance trigger."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from paths import plot_paths, INPUTS
root, figures = plot_paths()
scan = json.loads((root/'trigger_scan.json').read_text())
nodes = np.load(root/'trigger_nodes.npz')
frozen = np.load(INPUTS/'limber_outlier_tomography.npz')
wroot = scan[20]['w']
delta = np.array([s['offset']*1.e5 for s in scan])
native = np.array([[s['native']['qag_quad'], s['native']['spline']] for s in scan])
error = np.array([s['full']['estimated_error']/s['full']['tolerance'] for s in scan])

plt.rcParams.update({'font.size':11, 'axes.titlesize':12})
fig, axes = plt.subplots(2, 2, figsize=(12, 8.3), layout='constrained')
ax = axes[0,0]
for i, color in [(3,'#2864a0'), (4,'#ca5a28')]:
    y = frozen[f'nz_{i}']
    ax.plot(frozen['z'],y/y.max(),label=f'Lens bin {i+1}',color=color)
product = frozen['nz_3']*frozen['nz_4']
ax.plot(frozen['z'],product/product.max(),color='#6c3d91',ls='--',label='Overlap, rescaled')
ax.set(xlim=(.1,2.),xlabel='Redshift z',ylabel='Each curve / its maximum',title='Fixed redshift distributions: narrow overlap + tails')
ax.legend(frameon=False,fontsize=10)

ax = axes[0,1]
x,y = nodes['x'],nodes['y']
kn,ky,gn = nodes['kronrod_x'],nodes['kronrod_y'],nodes['gauss_x']
shared = np.array([np.min(abs(gn-q))<1.e-12 for q in kn])
ax.plot(x,y/y.max(),color='#2864a0',label='Native CCL integrand')
ax.scatter(kn[shared],ky[shared]/y.max(),color='#169e78',s=38,zorder=3,label='20 shared Gauss nodes')
ax.scatter(kn[~shared],ky[~shared]/y.max(),color='#ed7938',s=38,zorder=3,label='21 added Kronrod nodes')
ax.set(xlim=(-2.65,-1.9),ylim=(-.02,1.13),xlabel=r'$\ln[k/({\rm Mpc}^{-1})]$',ylabel='Integrand / its maximum',title='First pass: nodes miss the top of the peak')
ax.legend(frameon=False,fontsize=9,loc='upper right')

ax = axes[1,0]
ax.plot(delta,100*(native[:,0]/native[:,1]-1),'o-',color='#ca5a28',ms=3,label='Unmodified CCL QAG / spline − 1')
ax.axhline(0,color='0.4',lw=.8)
ax.set(xlabel=r'$(w_0-w_\star)\times 10^5$',ylabel='Difference [%]',title='The accepted result jumps; the spline stays smooth')
ax.text(.03,.14,'QAG: 41 evaluations inside the dip\n287 evaluations outside',transform=ax.transAxes,fontsize=10)
ax.legend(frameon=False,fontsize=9,loc='upper right')

ax = axes[1,1]
ax.semilogy(delta,error,'o-',color='#6c3d91',ms=3)
ax.axhline(1,color='black',ls='--',lw=1,label='First-pass acceptance threshold')
ax.fill_between(delta,1.e-11,1,color='#d8eee4',alpha=.7)
ax.set(ylim=(1.e-11,1.e3),xlabel=r'$(w_0-w_\star)\times 10^5$',ylabel='GSL estimated error / requested tolerance',title='Accidental rule agreement defeats the error estimate')
ax.text(.04,.09,'At the crossing, both initial rules\nare 4.0225% below CQUAD.',transform=ax.transAxes,fontsize=10)
ax.legend(frameon=False,fontsize=9,loc='upper left')
for ax in axes.flat:
    ax.grid(alpha=.2)
fig.suptitle(r'CCL 3.3.3 · GSL 2.7 · $\ell=355.655882$ · lens bins 4 × 5 (one-based)'
             '\n'+rf'$w_\star={wroot:.10f}$; all other cosmology and survey inputs fixed',fontsize=13)
fig.savefig(figures/'qag_trigger.png',dpi=180)
fig.savefig(figures/'qag_trigger.pdf')
plt.close(fig)
