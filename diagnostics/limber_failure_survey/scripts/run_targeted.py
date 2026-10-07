"""Sequential launcher for the bounded nine-plane search and saved candidates.

Only the coordinating agent/user should launch numerical stages. Run search
and CCL with the CCL runtime, and CoCoA with its normal activated environment.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--engine', choices=['search', 'ccl', 'cocoa', 'summary'], required=True)
    parser.add_argument('--planes', nargs='+', type=int, default=list(range(9)))
    parser.add_argument('--family', default='pr_original')
    args = parser.parse_args()
    base = args.output.resolve(); here = Path(__file__).resolve().parent
    if len(args.planes) != len(set(args.planes)) or any(p not in range(9) for p in args.planes):
        raise ValueError('Planes must be distinct integers from 0 to 8')
    env = dict(os.environ, OMP_NUM_THREADS='6', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    for plane in args.planes:
        folder = base/'targeted'/f'{args.family}_plane{plane:02d}'
        if args.engine == 'search':
            if folder.exists():
                record = json.loads((folder/'search.json').read_text()) if (folder/'search.json').exists() else {}
                if record.get('status') not in ('completed','budget_exhausted'):
                    raise RuntimeError(f'Preserved incomplete search needs inspection: {folder}')
                print(f'Reuse completed search: {folder}', flush=True)
                continue
            command = [sys.executable, str(here/'hunt.py'), '--output', str(base), '--plane', str(plane),
                       '--family', args.family, '--maximum-roots', '1']
        else:
            manifest = json.loads((folder/'manifest.json').read_text())
            if not manifest['cases']:
                print(f'No candidate in plane {plane}; retain empty search.', flush=True)
                continue
            if args.engine == 'summary':
                command = [sys.executable, str(here/'summarize.py'), '--output', str(folder)]
            else:
                command = [sys.executable, str(here/'run_batch.py'), '--engine', args.engine,
                           '--pilot', '--resume', '--output', str(folder)]
        print(' '.join(command), flush=True)
        subprocess.run(command, check=True, env=env)


if __name__ == '__main__':
    main()
