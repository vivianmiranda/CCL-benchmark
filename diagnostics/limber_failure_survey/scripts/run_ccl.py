"""One CCL correctness case; save status, references and matched CoCoA inputs."""
import argparse
import json
from pathlib import Path
import platform
import sys
import time
import traceback
import warnings

import numpy as np
import pyccl as ccl

from common import digest, read_manifest, selected_case, write_json, source_hashes
from ccl_helpers import make_cosmology, make_tracers, native_spectra, GSLRules, reference, warning_records


def export_power(cosmo, case, output):
    folder = output/'shared'
    folder.mkdir(exist_ok=True)
    filename = folder/f'{case["cosmology_id"]}.npz'
    meta = filename.with_suffix('.json')
    if filename.exists():
        old = json.loads(meta.read_text())
        assert old['parameters'] == case['parameters']
        assert old['camb_dark_energy_model'] == case['camb_dark_energy_model']
        assert old['sha256'] == digest(filename)
        return old
    h = case['parameters']['h']
    k = np.geomspace(1e-5, 1e3, 1500)
    zp = np.linspace(0,4.4,221); zd = np.linspace(0,50.1,10021); zg = np.linspace(0,4.4,881)
    linear = np.array([ccl.linear_matter_power(cosmo,k,1/(1+z)) for z in zp])*h**3
    nonlinear = np.array([ccl.nonlin_matter_power(cosmo,k,1/(1+z)) for z in zp])*h**3
    growth = ccl.growth_factor(cosmo,1/(1+zg))*(1+zg); growth /= growth[-1]
    np.savez_compressed(filename, log10k_2D=np.log10(k/h),z_2D=zp,
        lnP_linear=np.log(linear).ravel(order='F'),lnP_nonlinear=np.log(nonlinear).ravel(order='F'),
        lnP_linear_cb=np.log(linear).ravel(order='F'),G=growth,z_G=zg,z_1D=zd,
        chi=ccl.comoving_radial_distance(cosmo,1/(1+zd))*h,omegan2=np.array(0.))
    record = {'parameters': case['parameters'], 'camb_dark_energy_model': case['camb_dark_energy_model'],
              'sha256': digest(filename), 'path': str(filename.relative_to(output)),
              'seed_k_nodes':1500,'cocoa_refinement':8,'dense_k_nodes':11993}
    write_json(meta,record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--case',required=True)
    args = parser.parse_args(); output=args.output.resolve()
    manifest,config=read_manifest(output); case=selected_case(manifest,args.case)
    destination=output/'ccl'/f'{case["id"]}.json'
    if destination.exists():raise FileExistsError(destination)
    record={'case':case,'config_sha256':manifest['config_sha256'],'status':'started',
            'ccl_version':ccl.__version__,'python':sys.version,'platform':platform.platform(),
            'script_sha256':source_hashes()}
    start=time.perf_counter(); write_json(destination,record)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            inp=output/'inputs'/f'{case["family"]}.npz'
            assert digest(inp)==manifest['families'][case['family']]['npz_sha256']
            data=np.load(inp);z,nz,centers=data['z'],data['nz'],data['centers']
            cosmo=make_cosmology(case); tracers,bias=make_tracers(cosmo,z,nz,centers)
            ell=config['ell']; native=native_spectra(cosmo,tracers,ell)
            record.update(native=native,bias=bias.tolist(),ell=ell,input_sha256=digest(inp))
            rules=GSLRules(); record.update(gsl_version=rules.version,gsl_library=rules.path)
            primary=config['primary_pair_zero_based']
            candidates=[row for row in native if row['pair']==primary]
            for row in native:
                q,s=row['qag_quad']['value'],row['spline']['value']
                if row['pair']!=primary and q is not None and s is not None and s>0:
                    if abs(q/s-1)>config['screen_additional_pair_threshold']:candidates.append(row)
            checks=[]
            for row in candidates:
                context=rules.context(cosmo,tracers,row['pair'],ell)
                result=reference(rules,context,config)
                result['pair']=row['pair']; result['primary']=row['pair']==primary
                q=row['qag_quad']['value']; replay=result['native_qag_replay']['value']
                result['native_QAG_replay_matches']=bool(q is not None and replay is not None
                    and np.isclose(q,replay,rtol=1e-11,atol=0))
                result['reference_valid'] &= result['native_QAG_replay_matches']
                value=result.get('reference_value',result['cquad']['value'])
                result['native_QAG_fractional_error']=q/value-1 if q is not None and value is not None and value>0 else None
                if result['primary']:
                    initial,kn,gn=rules.initial_rules(context)
                    x=np.linspace(context['lo'],context['hi'],5001)
                    y=np.array([context['evaluate'](float(v),None) for v in x])
                    ky=np.array([context['evaluate'](float(v),None) for v in kn])
                    (output/'ccl').mkdir(exist_ok=True)
                    np.savez_compressed(output/'ccl'/f'{case["id"]}_integrand.npz',
                        logk=x,integrand=y,kronrod_logk=kn,kronrod_integrand=ky,gauss_logk=gn)
                    result['initial_rules']=initial
                    if manifest['families'][case['family']]['definition']['kind']=='positive_mixture':
                        fine_path=output/'inputs'/f'{case["family"]}_refined.npz'
                        assert digest(fine_path)==manifest['families'][case['family']]['refined_npz_sha256']
                        fine=np.load(fine_path)
                        ft,_=make_tracers(cosmo,fine['z'],fine['nz'],fine['centers'])
                        refined=rules.integrate(rules.context(cosmo,ft,primary,ell),'cquad',config['reference_epsrel'])
                        result['input_grid_refined_reference']=refined
                        good=refined['value'] is not None and not any(refined['status']) and not refined['callback_errors']
                        change=abs(refined['value']/value-1) if good and value is not None and value>0 else None
                        result['input_grid_fractional_change']=change
                        result['input_grid_valid']=bool(change is not None and change<=config['input_refinement_fractional_tolerance'])
                        result['reference_valid'] &= result['input_grid_valid']
                checks.append(result)
            record['references']=checks
            record['shared_power']=export_power(cosmo,case,output)
            record['status']='completed'
        record['warnings']=warning_records(caught)
    except Exception as error:
        if 'caught' in locals():record['warnings']=warning_records(caught)
        record.update(status='setup_or_numerical_exception',exception=repr(error),traceback=traceback.format_exc())
        raise
    finally:
        record['elapsed_seconds_for_resource_planning']=time.perf_counter()-start
        record['elapsed_is_not_benchmark']=True
        write_json(destination,record)
    print(json.dumps({'case':case['id'],'status':record['status'],
        'elapsed_seconds':record['elapsed_seconds_for_resource_planning'],
        'primary':next(r for r in record['references'] if r['primary'])},indent=2))


if __name__=='__main__':main()
