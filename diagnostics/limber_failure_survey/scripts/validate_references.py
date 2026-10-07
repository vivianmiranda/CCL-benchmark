"""Supplement saved reference checks without rewriting original results."""
import argparse
import json
from pathlib import Path
import time
import warnings

import numpy as np
from common import read_manifest, selected_case, write_json, digest, source_hashes
from ccl_helpers import make_cosmology, make_tracers, GSLRules, supplement_reference, warning_records


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--case',required=True)
    args=parser.parse_args();output=args.output.resolve();manifest,config=read_manifest(output)
    case=selected_case(manifest,args.case)
    original_path=output/'ccl'/f'{case["id"]}.json'
    original=json.loads(original_path.read_text())
    assert original['status']=='completed' and original['case']==case
    destination=output/'reference_supplements'/f'{case["id"]}.json'
    if destination.exists():raise FileExistsError(destination)
    start=time.perf_counter();data=np.load(output/'inputs'/f'{case["family"]}.npz')
    checks=[]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        cosmo=make_cosmology(case);tracers,_=make_tracers(cosmo,data['z'],data['nz'],data['centers'])
        rules=GSLRules()
        for old in original['references']:
            context=rules.context(cosmo,tracers,old['pair'],config['ell'])
            result=supplement_reference(rules,context,config,old)
            result.update(pair=old['pair'],primary=old['primary'])
            result['reference_valid'] &= old['native_QAG_replay_matches']
            if 'input_grid_refined_reference' in old:
                fine=old['input_grid_refined_reference'];value=result['reference_value']
                success=fine['value'] is not None and not any(fine['status']) and not fine['callback_errors']
                change=abs(fine['value']/value-1) if success and value is not None and value>0 else None
                result['input_grid_fractional_change']=change
                result['input_grid_valid']=bool(change is not None and change<=config['input_refinement_fractional_tolerance'])
                result['reference_valid'] &= result['input_grid_valid']
            checks.append(result)
    record={'case':case,'status':'completed','original_result_sha256':digest(original_path),
        'config_sha256':manifest['config_sha256'],'script_sha256':source_hashes(),
        'references':checks,'warnings':warning_records(caught),
        'elapsed_seconds_for_resource_planning':time.perf_counter()-start,'elapsed_is_not_benchmark':True}
    write_json(destination,record)
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
