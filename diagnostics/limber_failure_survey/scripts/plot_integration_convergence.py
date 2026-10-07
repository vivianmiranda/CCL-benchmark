"""Slack-ready quadrature convergence figure from measured CoCoA results."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter
from common import digest, write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--case',required=True)
    a=p.parse_args()
    path=a.output/'cocoa_integration_ladder'/f'{a.case}.json'
    data=json.loads(path.read_text())
    check=json.loads((path.parent/'regression_check.json').read_text())
    assert check['successful'] and check['source_sha256']==digest(path)
    row=next(r for r in data['pairs'] if r['pair']==[3,4])
    errors=100*abs(np.asarray(row['values'])/row['values'][-1]-1)
    worst=100*np.max([abs(np.asarray(r['values'])/r['values'][-1]-1) for r in data['pairs']],axis=0)
    plt.rcParams.update({'font.size':14,'axes.labelsize':16,'axes.titlesize':18})
    fig,ax=plt.subplots(figsize=(11.5,8.2))
    fig.subplots_adjust(left=.13,right=.97,top=.76,bottom=.32)
    fig.suptitle('CoCoA: increasing integration_accuracy resolves the error',fontsize=20,y=.975,weight='bold')
    fig.text(.5,.915,'The same narrow-overlap stress test: σz = 0.003 at z = 0.94087',ha='center',fontsize=15)
    fig.text(.5,.868,'Quadrature error relative to integration_accuracy = 4 (1024 nodes per panel)',ha='center',fontsize=14)
    ax.plot(range(4),worst[:4],color='#9e6525',ls='--',marker='s',ms=7,lw=1.8,label='Largest error among all 15 bin pairs')
    ax.plot(range(4),errors[:4],color='#236ca4',marker='o',ms=11,lw=2.7,label='Affected pair: lens bins 4 × 5')
    for k,value in enumerate(errors[:4]):
        ax.annotate(f'{value:.6g}%',(k,value),xytext=(0,-23 if k<3 else -25),textcoords='offset points',ha='center',fontsize=14,color='#236ca4',weight='bold')
    ax.axhline(.01,color='#37804a',ls=':',lw=2,label='Convergence check: 0.01%')
    ax.axvspan(1.65,3.4,color='#37804a',alpha=.065)
    ax.text(2.4,.12,'Levels 2 and 3 pass\nfor all 15 pairs',ha='center',fontsize=16,color='#286b3b',weight='bold')
    ax.set(yscale='log',ylim=(.00007,1.8),xlim=(-.35,3.35),
           xticks=range(4),xticklabels=['0\n96 nodes','1\n128 nodes','2\n256 nodes','3\n512 nodes'],
           xlabel='integration_accuracy  ·  nodes per panel',ylabel='Absolute quadrature difference [%]')
    ax.set_yticks([1,.1,.01,.001,.0001])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y,_:f'{y:g}'))
    ax.grid(axis='y',which='major',alpha=.18);ax.spines[['right','top']].set_visible(False)
    ax.legend(loc='upper right',bbox_to_anchor=(1,1.16),frameon=False,fontsize=11)
    fig.text(.13,.17,'Default → level 2: 0.52036% → 0.00152% for the affected pair.',fontsize=16,color='#236ca4',weight='bold')
    fig.text(.13,.116,'A default-versus-refined regression check rejects level 0; level 2 passes 0.01%.',fontsize=12)
    fig.text(.13,.077,'Level 3 vs 4 also passes 0.001% for every pair. Inputs, interpolation and the seven panels stay fixed.',fontsize=11)
    fig.text(.13,.038,'This checks CoCoA covariance-module spectra, not ordinary data-vector integration or universal immunity.',fontsize=11)
    folder=a.output/'figures';folder.mkdir(exist_ok=True)
    for ext in ('png','pdf'):fig.savefig(folder/f'integration_accuracy_convergence.{ext}',dpi=200)
    plt.close(fig)
    write_json(folder/'integration_accuracy_convergence.json',{'source_sha256':digest(path),'script_sha256':digest(__file__),
        'levels':[0,1,2,3,4],'nodes':data['nodes'],'primary_error_percent':errors.tolist(),'worst_pair_error_percent':worst.tolist(),
        'quadrature_only':True})


if __name__=='__main__':main()
