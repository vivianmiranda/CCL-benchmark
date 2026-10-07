"""Bounded initial-rule crossing search; never an incidence sample.

Run one predeclared Omega_m/wa plane per process. Candidate roots require
the same full CCL and CoCoA validation as every other plotted case.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import shutil
import time

import numpy as np
from scipy.optimize import brentq

from common import read_manifest, write_json, digest
from ccl_helpers import make_cosmology, make_tracers, GSLRules


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True,help='Prepared uniform-study folder')
    parser.add_argument('--plane',type=int,required=True)
    parser.add_argument('--family',default='pr_original')
    parser.add_argument('--maximum-roots',type=int,default=1,
                        help='Resource cap, at most the predeclared maximum (default: one per plane)')
    parser.add_argument('--saved-search',type=Path,
                        help='Read verified samples from a preserved previous search; never modify it')
    args=parser.parse_args();base=args.output.resolve();manifest,config=read_manifest(base)
    search=config['targeted_search'];planes=[(om,wa) for om in search['Omega_m'] for wa in search['wa']]
    if not 0<=args.plane<len(planes):raise ValueError(f'plane must be in 0..{len(planes)-1}')
    if args.family not in manifest['families']:raise ValueError(args.family)
    if not 1<=args.maximum_roots<=search['maximum_roots_per_plane']:
        raise ValueError('maximum-roots exceeds the predeclared bound')
    om,wa=planes[args.plane]
    grid=np.linspace(*config['uniform_ranges']['w0'],search['w0_bracket_nodes'])
    cache={};unsuccessful=[];reuse=None
    if args.saved_search:
        path=args.saved_search.resolve();old=json.loads(path.read_text())
        old_manifest,old_config=read_manifest(path.parent)
        if old_manifest['config_sha256']!=manifest['config_sha256'] or old_config!=config:
            raise ValueError('Saved search used a different configuration')
        if old['Omega_m']!=om or old['wa']!=wa or old['family']!=args.family:
            raise ValueError('Saved search belongs to a different plane or input family')
        family=manifest['families'][args.family]
        if old_manifest['families'][args.family]!=family:
            raise ValueError('Saved family metadata differ')
        for parent in (base,path.parent):
            if digest(parent/'inputs'/f'{args.family}.npz')!=family['npz_sha256']:
                raise ValueError('Saved/current family input fingerprint differs')
        for item in old['evaluations']:
            item=copy.deepcopy(item);w=item['w0']
            if not grid[0]<=w<=grid[-1] or not all(math.isfinite(item[k]) for k in ('w0','kronrod','gauss','estimated_error','signed_rule_difference')):
                raise ValueError('Invalid cached initial-rule sample')
            if item['callback_errors'] or not math.isclose(item['signed_rule_difference'],item['kronrod']-item['gauss'],rel_tol=1e-12,abs_tol=1e-30):
                raise ValueError('Cached initial-rule sample is inconsistent')
            item.setdefault('integration_logk_bounds',None)
            cache[float(w)]=item
        unsuccessful=copy.deepcopy(old.get('unsuccessful_brackets',[]))
        if old.get('status')=='stopped' and 'Failed to converge' in str(old.get('problem')) and not unsuccessful:
            # Legacy searches stopped at the first failed bracket. The final
            # solver sample identifies its enclosing predeclared grid interval.
            final=old['evaluations'][-1]['w0'];index=int(np.searchsorted(grid,final))-1
            left,right=float(grid[index]),float(grid[index+1])
            if not 0<=index<len(grid)-1 or left not in cache or right not in cache or cache[left]['signed_rule_difference']*cache[right]['signed_rule_difference']>=0:
                raise ValueError('Cannot safely infer the legacy unsuccessful bracket')
            unsuccessful.append({'bounds_w0':[left,right],'reason':'legacy_solver_did_not_converge',
                'exception':old['problem'],'inferred_from_saved_trace':True})
        reuse={'search_path':str(path),'search_sha256':digest(path),'manifest_sha256':digest(path.parent/'manifest.json'),
               'config_sha256':manifest['config_sha256'],'input_sha256':family['npz_sha256'],
               'reused_evaluations':len(cache),'skipped_saved_brackets':copy.deepcopy(unsuccessful)}
    folder=base/'targeted'/f'{args.family}_plane{args.plane:02d}'
    if folder.exists():raise FileExistsError(folder)
    folder.mkdir(parents=True);shutil.copyfile(base/'config.json',folder/'config.json')
    shutil.copytree(base/'inputs',folder/'inputs')
    data=np.load(base/'inputs'/f'{args.family}.npz')
    rules=GSLRules();start=time.perf_counter();roots=[]
    params=dict(config['fixed_cosmology'],Omega_m=om,wa=wa)
    budget=search['maximum_function_evaluations_per_plane']
    def evaluate(w):
        key=float(w)
        if key in cache:return cache[key]['signed_rule_difference']
        if len(cache)>=budget or time.perf_counter()-start>config['case_wall_seconds']:
            raise TimeoutError('Targeted plane budget exhausted')
        case={'parameters':dict(params,w0=key),'camb_dark_energy_model':config['survey_camb_dark_energy_model']}
        cosmo=make_cosmology(case);tracers,_=make_tracers(cosmo,data['z'],data['nz'],data['centers'])
        context=rules.context(cosmo,tracers,config['primary_pair_zero_based'],config['ell'])
        initial,_,_=rules.initial_rules(context)
        if initial['callback_errors'] or not all(initial.get(k) is not None and math.isfinite(initial[k]) for k in ('kronrod','gauss','estimated_error')):
            raise ValueError(initial)
        initial.update(w0=key,signed_rule_difference=initial['kronrod']-initial['gauss'],
                       integration_logk_bounds=[context['lo'],context['hi']])
        cache[key]=initial
        write_json(folder/'search_progress.json',{'Omega_m':om,'wa':wa,'family':args.family,
            'evaluations':list(cache.values()),'candidate_roots':roots,
            'unsuccessful_brackets':unsuccessful,'saved_search_reuse':reuse,'excluded_from_incidence':True})
        return initial['signed_rule_difference']
    status='completed';problem=None
    try:
        previous=None
        for w in grid:
            value=evaluate(w)
            if previous and value*previous[1]<0:
                bounds=[previous[0],float(w)]
                if any(np.allclose(row['bounds_w0'],bounds,rtol=0,atol=1e-14) for row in unsuccessful):
                    previous=(float(w),value);continue
                try:
                    root=float(brentq(evaluate,*bounds,xtol=1e-10,maxiter=30))
                    evaluate(root);item=cache[root]
                    nearby=sorted((key,row) for key,row in cache.items() if bounds[0]<=key<=bounds[1])
                    jumps=[]
                    for (x,a),(y,b) in zip(nearby[:-1],nearby[1:]):
                        if y-x<=1e-8 and x-1e-10<=root<=y+1e-10 and a['signed_rule_difference']*b['signed_rule_difference']<0:
                            fractions=[abs(v['signed_rule_difference']/v['kronrod']) if v['kronrod'] else None for v in (a,b)]
                            if all(v is not None for v in fractions) and min(fractions)>1e-3:
                                jumps.append({'bounds_w0':[x,y],'fractional_rule_differences':fractions})
                    accepted=bool(item['kronrod']>0 and item['estimated_error']<=1e-4*abs(item['kronrod'])
                                  and abs(item['signed_rule_difference'])<=1e-4*abs(item['kronrod']) and not jumps)
                    if accepted:
                        if not roots or abs(root-roots[-1])>1e-8:roots.append(root)
                    else:
                        unsuccessful.append({'bounds_w0':bounds,'reason':'initial_GSL_estimate_or_discontinuous_root_rejected',
                            'solver_candidate_w0':root,'initial_rules':item,'nearby_value_jumps':jumps})
                except (RuntimeError,ValueError) as error:
                    unsuccessful.append({'bounds_w0':bounds,'reason':'individual_bracket_solver_failed','exception':repr(error)})
                if len(roots)>=args.maximum_roots:break
            previous=(float(w),value)
    except TimeoutError as error:
        status='budget_exhausted';problem=str(error)
    except Exception as error:
        status='stopped';problem=repr(error)
    new=copy.deepcopy(manifest);new['cases']=[];new['pilot_case_ids']=[]
    for i,w in enumerate(roots):
        ident=f't{args.plane:02d}r{i:02d}__{args.family}'
        new['cases'].append({'id':ident,'cosmology_id':f't{args.plane:02d}r{i:02d}',
            'parameters':dict(params,w0=w),'family':args.family,'arm':'targeted_search',
            'camb_dark_energy_model':config['survey_camb_dark_energy_model']})
        new['pilot_case_ids'].append(ident)
    new['sampling_description']='Deliberate initial-rule crossing search. Excluded from every failure-frequency denominator.'
    write_json(folder/'manifest.json',new)
    write_json(folder/'search.json',{'status':status,'problem':problem,'Omega_m':om,'wa':wa,
        'family':args.family,'evaluations':list(cache.values()),'candidate_roots':roots,
        'maximum_roots_for_this_run':args.maximum_roots,
        'unsuccessful_brackets':unsuccessful,'saved_search_reuse':reuse,
        'candidate_acceptance':'Finite positive K41; initial GSL error and |K41-G20| each <= 1e-4 |K41|; no detected local jump. Full native/reference checks still required.',
        'excluded_from_incidence':True,'elapsed_seconds_for_resource_planning':time.perf_counter()-start})
    print(json.dumps({'folder':str(folder),'status':status,'candidates':len(roots),'problem':problem},indent=2))
    if status=='stopped':raise SystemExit(1)


if __name__=='__main__':main()
