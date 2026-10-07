"""Actual CoCoA covariance-module spectra, with three fixed-rule resolutions."""
import argparse
import json
import os
from pathlib import Path
import sys
import time
import subprocess

import numpy as np
from common import digest, read_manifest, selected_case, write_json, finite_float, source_hashes


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--case',required=True)
    parser.add_argument('--integration-ladder',action='store_true',help='Use the five public integration_accuracy levels, saving a separate diagnostic')
    args=parser.parse_args();output=args.output.resolve()
    manifest,config=read_manifest(output);case=selected_case(manifest,args.case)
    folder=output/('cocoa_integration_ladder' if args.integration_ladder else 'cocoa');folder.mkdir(exist_ok=True)
    destination=folder/f'{case["id"]}.json'
    if destination.exists():raise FileExistsError(destination)
    ccl=json.loads((output/'ccl'/f'{case["id"]}.json').read_text())
    assert ccl['status']=='completed' and ccl['case']==case
    shared=output/ccl['shared_power']['path'];assert digest(shared)==ccl['shared_power']['sha256']
    nz_path=output/'inputs'/f'{case["family"]}.txt'
    assert digest(nz_path)==manifest['families'][case['family']]['text_sha256']
    import cosmolike_lsst_y1_interface as ci
    sys.path.insert(0,str(Path(os.environ['ROOTDIR'])/'external_modules/code/cosmolike_core'))
    from cosmolike_notebook_utils.covariance.power import refine_power_tables
    from cosmolike_notebook_utils.covariance.accuracy import covariance_accuracy
    t0=time.perf_counter(); b=np.load(shared);tables={name:b[name] for name in b.files}
    tables['omegan2']=float(tables['omegan2']);tables=refine_power_tables(tables,refinement=8)
    ci.initial_setup();ci.init_accuracy_boost(accuracy_boost=1.,integration_accuracy=0)
    ci.init_probes(possible_probes='3x2pt');ci.init_IA(ia_model=0,ia_redshift_evolution=2,ia_code=0)
    ci.init_bias(bias_model=[0]*5);ci.init_photoz_conventions(interpolation_type=1,zmid_convention=1)
    ci.init_cosmo_runmode(is_linear=False)
    ci.init_redshift_distributions_from_files(lens_multihisto_file=str(nz_path),lens_ntomo=5,
                                             source_multihisto_file=str(nz_path),source_ntomo=5)
    p=case['parameters'];ci.set_cosmology(omegam=p['Omega_m'],omegab=p['Omega_b'],H0=100*p['h'],**tables)
    zeros=[0.]*5
    ci.set_nuisance_bias(B1=ccl['bias'],B2=zeros,B_MAG=zeros,B3nl=zeros,BK=zeros)
    ci.set_nuisance_ia(A1=zeros,A2=zeros,B_TA=zeros);ci.set_nuisance_shear_photoz(bias=zeros)
    ci.set_nuisance_clustering_photoz(bias=zeros);ci.set_nuisance_shear_calib(M=zeros)
    ell=np.array([config['ell']],dtype=float)
    edges=1/(1+np.asarray(config['a_panel_redshifts']))
    spectra={};geometry=None
    controls = ([covariance_accuracy(accuracy_boost=1,integration_accuracy=i) for i in range(5)]
                if args.integration_ladder else [{'radial_nquad':n,'nwindow':4097} for n in config['cocoa_nodes_per_panel']])
    node_counts = [c['radial_nquad'] for c in controls]
    for control in controls:
        n=control['radial_nquad']
        if args.integration_ladder:
            ci.init_accuracy_boost(accuracy_boost=1.,integration_accuracy=control['integration_accuracy'])
        result=ci.covariance.covariance_spectra(ell=ell,a_edges=edges,nquad=n,nwindow=control['nwindow'],
                    include_ia=False,include_rsd=False,linear=False,nonlimber_lmax=0)
        spectra[str(n)]=result['spectra'][0,:5,:5]
        if n==96:geometry=result['geometry']
    pairs=np.array([(i,j) for i in range(5) for j in range(i,5)])
    values=np.array([[spectra[str(n)][i,j] for i,j in pairs] for n in node_counts])
    np.savez_compressed(folder/f'{case["id"]}.npz',ell=ell,pairs=pairs,spectra=values,
        nodes=np.array(node_counts),geometry_96=geometry,a_edges=edges)
    rows=[]
    for pair,column in zip(pairs,values.T):
        valid=bool(np.all(np.isfinite(column)) and column[-1]>0)
        refine=abs(column[-2]/column[-1]-1) if valid else None
        valid = bool(valid and refine is not None and refine<=config['cocoa_reference_consistency_fractional_tolerance'])
        row={'pair':pair.tolist(),'values':[finite_float(v) for v in column],
            'reference_valid':valid,'512_vs_1024_fractional_change':finite_float(refine) if refine is not None else None,
            'default_fractional_change':finite_float(column[0]/column[-1]-1) if np.isfinite(column).all() and column[-1]>0 else None}
        cr=next((r for r in ccl['references'] if r['pair']==pair.tolist()),None)
        if cr and cr['reference_valid']:
            value=cr.get('reference_value',cr['cquad']['value'])
            row['common_CCL_CQUAD_reference_value']=value
            row['default_offset_from_CCL_CQUAD']=finite_float(column[0]/value-1)
            row['refined_offset_from_CCL_CQUAD']=finite_float(column[-1]/value-1)
        rows.append(row)
    write_json(destination,{'status':'completed','case':case,'config_sha256':manifest['config_sha256'],
        'script_sha256':source_hashes(),
        'core_commit':subprocess.check_output(['git','-C',str(Path(os.environ['ROOTDIR'])/'external_modules/code/cosmolike_core'),'rev-parse','HEAD'],text=True).strip(),
        'lsst_y1_commit':subprocess.check_output(['git','-C',str(Path(os.environ['ROOTDIR'])/'projects/lsst_y1'),'rev-parse','HEAD'],text=True).strip(),
        'shared_power_sha256':digest(shared),'nz_text_sha256':digest(nz_path),'pairs':rows,
        'nodes':node_counts,
        'integration_ladder': args.integration_ladder,
        'controls':[{k:c[k] for k in ('radial_nquad','nwindow','integration_accuracy','accuracy_boost') if k in c} for c in controls],
        'scope':'Density-only all-pairs spectra from covariance module; not ordinary data-vector integration; no covariance matrix computed.',
        'elapsed_seconds_for_resource_planning':time.perf_counter()-t0,'elapsed_is_not_benchmark':True})
    print(json.dumps({'case':case['id'],'status':'completed','primary':next(r for r in rows if r['pair']==config['primary_pair_zero_based'])},indent=2))


if __name__=='__main__':main()
