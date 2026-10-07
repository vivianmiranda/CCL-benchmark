"""Sequential, bounded fresh-process launcher for either runtime."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from common import read_manifest, write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--engine',choices=['ccl','cocoa'],required=True)
    choose=parser.add_mutually_exclusive_group(required=True)
    choose.add_argument('--pilot',action='store_true');choose.add_argument('--uniform',action='store_true')
    parser.add_argument('--resume',action='store_true',help='Reuse completed cases after checking their configuration hash')
    args=parser.parse_args();output=args.output.resolve();manifest,config=read_manifest(output)
    cases=manifest['pilot_case_ids'] if args.pilot else [c['id'] for c in manifest['cases'] if c['arm']=='uniform']
    env=dict(os.environ,OMP_NUM_THREADS=str(config['threads']),OPENBLAS_NUM_THREADS='1',
             MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    logs=output/'logs';logs.mkdir(exist_ok=True)
    script=Path(__file__).with_name(f'run_{args.engine}.py')
    ledger=[];ledger_path=output/f'{args.engine}_{"pilot" if args.pilot else "uniform"}_status.json'
    for case in cases:
        saved=output/args.engine/f'{case}.json'
        if saved.exists():
            old=json.loads(saved.read_text())
            if args.resume and old['status']=='completed' and old['config_sha256']==manifest['config_sha256']:
                ledger.append({'case':case,'reused_completed':True});write_json(ledger_path,ledger);continue
            raise FileExistsError(f'Preserve {saved}; inspect before resuming')
        log=logs/f'{args.engine}_{case}.log'
        if log.exists():raise FileExistsError(f'Preserve {log}; inspect before retrying')
        print(f'{args.engine}: {case}',flush=True)
        with log.open('w') as out:
            try:
                result=subprocess.run([sys.executable,str(script),'--output',str(output),'--case',case],
                    env=env,stdout=out,stderr=subprocess.STDOUT,timeout=config['case_wall_seconds'])
                code=result.returncode
            except subprocess.TimeoutExpired:code=124
        ledger.append({'case':case,'returncode':code,'log':str(log.relative_to(output))})
        write_json(ledger_path,ledger)
        if code:raise SystemExit(f'{case} failed with {code}; inspect {log}')


if __name__=='__main__':main()
