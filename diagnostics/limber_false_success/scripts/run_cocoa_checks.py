"""Run only the documented CoCoA checks, sequentially in fresh processes."""
import json
import os
from pathlib import Path
import subprocess
import sys

from collect_cocoa import collect
from paths import output_dir


def main():
    root = output_dir()
    if list(root.glob('cocoa_case*.npz')):
        raise FileExistsError('CoCoA outputs exist; choose a fresh directory')
    for case in range(3):
        if not (root/f'shared_input_{case}.npz').exists():
            raise FileNotFoundError('Run export_ccl_cases.py in the CCL environment first')
    threads = int(os.environ.get('OMP_NUM_THREADS', '8'))
    if not 1 <= threads <= 8:
        raise ValueError('Set OMP_NUM_THREADS between 1 and 8')
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    script = Path(__file__).with_name('check_cocoa_case.py')
    for case in range(3):
        for level in range(5):
            command = [sys.executable, str(script), str(case), str(level)]
            if case == 1 and level == 0:
                command.append('--save-nodes')
            print(f'CoCoA case {case}, integration level {level}', flush=True)
            with (root/f'cocoa_case{case}_level{level}.log').open('w') as log:
                subprocess.run(command, env=env, stdout=log,
                               stderr=subprocess.STDOUT, check=True, timeout=600)
    (root/'cocoa_convergence.json').write_text(json.dumps(collect(root), indent=2)+'\n')


if __name__ == '__main__':
    main()
