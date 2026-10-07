"""Freeze a paired random design before inspecting its numerical outcomes."""
import argparse
from pathlib import Path
import shutil

import numpy as np

from common import PACKAGE, FROZEN, build_family, digest, write_json, source_hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--config', type=Path, default=PACKAGE/'config.json')
    args = parser.parse_args()
    if (args.output/'manifest.json').exists():
        raise FileExistsError('Sampling manifest already exists; preserve it')
    args.output.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.config, args.output/'config.json')
    import json
    config = json.loads(args.config.read_text())
    rng = np.random.Generator(np.random.PCG64(config['seed']))
    unit = rng.random((config['uniform_cosmologies'], 3))
    cosmologies = []
    cases = []
    for number, row in enumerate(unit):
        params = dict(config['fixed_cosmology'])
        for name, value in zip(config['parameter_order'], row):
            low, high = config['uniform_ranges'][name]
            params[name] = float(low+(high-low)*value)
        cosmologies.append({'id': f'u{number:03d}', 'parameters': params})
    control = dict(config['fixed_cosmology'], **config['known_control'])
    cosmologies.append({'id': 'known', 'parameters': control})
    families = {}
    for definition in config['families']:
        name = definition['name']
        folder = args.output/'inputs'
        folder.mkdir(exist_ok=True)
        z, nz, centers = build_family(definition)
        path = folder/f'{name}.npz'
        np.savez_compressed(path, z=z, nz=nz, centers=centers)
        np.savetxt(folder/f'{name}.txt', np.column_stack([z, *nz]))
        families[name] = {'definition': definition, 'npz_sha256': digest(path),
                          'text_sha256': digest(folder/f'{name}.txt'),
                          'nodes': len(z), 'minimum_value': float(nz.min())}
        if definition['kind'] == 'positive_mixture':
            z, nz, centers = build_family(definition, refined=True)
            path = folder/f'{name}_refined.npz'
            np.savez_compressed(path, z=z, nz=nz, centers=centers)
            families[name]['refined_npz_sha256'] = digest(path)
    for cosmology in cosmologies:
        for name in families:
            if cosmology['id'] == 'known' and name != 'pr_original':
                continue
            cases.append({'id': f'{cosmology["id"]}__{name}',
                'cosmology_id': cosmology['id'], 'parameters': cosmology['parameters'],
                'family': name, 'arm': 'uniform' if cosmology['id'] != 'known' else 'targeted_control',
                'camb_dark_energy_model': config['survey_camb_dark_energy_model']
                    if cosmology['id'] != 'known' else config['known_control_camb_dark_energy_model']})
    write_json(args.output/'manifest.json', {
        'config_sha256': digest(args.output/'config.json'),
        'frozen_input_sha256': digest(FROZEN),
        'prepared_script_sha256': source_hashes(),
        'sampling_description': 'Independent uniform cosmology draws; same draws paired across families. Targeted controls never enter incidence denominators.',
        'primary_statistic': 'Per-family fraction of valid primary-pair spectra exceeding each prespecified relative-integration-error threshold; not a posterior or a general failure probability.',
        'families': families, 'cases': cases,
        'pilot_case_ids': ['known__pr_original', 'u000__pr_original', 'u000__shared_narrow_peak']})
    print(f'Frozen {len(cases)} cases; run the three pilot cases before the full study.')


if __name__ == '__main__':
    main()
