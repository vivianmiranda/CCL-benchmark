"""Configuration and immutable sampling inputs; imports no cosmology code."""
import hashlib
import json
from pathlib import Path

import numpy as np

PACKAGE = Path(__file__).resolve().parents[1]
FROZEN = PACKAGE.parent/'limber_false_success/inputs/limber_outlier_tomography.npz'
VENDOR = PACKAGE.parent/'limber_false_success/vendor'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    return {p.name: digest(p) for p in sorted((PACKAGE/'scripts').glob('*.py'))}


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temp.replace(path)


def read_manifest(output):
    output = Path(output)
    result = json.loads((output/'manifest.json').read_text())
    assert digest(output/'config.json') == result['config_sha256']
    assert digest(FROZEN) == result['frozen_input_sha256']
    return result, json.loads((output/'config.json').read_text())


def build_family(definition, refined=False):
    """Keep broad supports while explicitly adding a shared narrow population."""
    with np.load(FROZEN) as original:
        z0 = original['z'].copy()
        nz0 = np.array([original[f'nz_{i}'] for i in range(5)])
        centers = original['bin_centers'].copy()
    if definition['kind'] == 'frozen_pr1313':
        return z0, nz0, centers
    n = definition['refined_z_nodes'] if refined else definition['z_nodes']
    z = np.linspace(z0[0], z0[-1], n)
    broad = np.array([np.interp(z, z0, row, left=0, right=0) for row in nz0])
    trapz = getattr(np, 'trapezoid', np.trapz) if hasattr(np, 'trapz') else np.trapezoid
    broad /= trapz(broad, z, axis=1)[:, None]
    peak = np.exp(-0.5*((z-definition['shared_mean_z'])/definition['shared_sigma_z'])**2)
    peak /= trapz(peak, z)
    fraction = definition['shared_fraction']
    nz = (1-fraction)*broad + fraction*peak[None, :]
    nz /= trapz(nz, z, axis=1)[:, None]
    assert np.all(np.isfinite(nz)) and np.all(nz >= 0)
    assert definition['shared_sigma_z'] / (z[1]-z[0]) >= 10
    return z, nz, centers


def selected_case(manifest, case_id):
    matches = [case for case in manifest['cases'] if case['id'] == case_id]
    if len(matches) != 1:
        raise ValueError(f'Unknown or duplicate case {case_id}')
    return matches[0]


def finite_float(value):
    value = float(value)
    return value if np.isfinite(value) else None
