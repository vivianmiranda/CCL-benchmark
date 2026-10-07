"""Recompute the CoCoA convergence table from saved, unaltered spectra."""
import json
import numpy as np

from paths import output_dir


def collect(root):
    records = []
    for case in range(3):
        reference = np.load(root/f'cocoa_case{case}_level4.npz')
        ell_index = np.flatnonzero(reference['ell'] == 355.655882).item()
        levels = []
        for level in range(5):
            current = np.load(root/f'cocoa_case{case}_level{level}.npz')
            np.testing.assert_array_equal(current['ell'], reference['ell'])
            for key in current.files:
                if not np.all(np.isfinite(current[key])):
                    raise ValueError(f'Non-finite {key} in case {case}, level {level}')
            gg = current['covariance_spectra'][:, :5, :5]
            gg_ref = reference['covariance_spectra'][:, :5, :5]
            relative = gg / gg_ref - 1
            levels.append({
                'level': level, 'cov_nodes': int(current['nodes']),
                'covariance_gg_max_abs_fractional_change_vs_level4':
                    float(np.max(np.abs(relative))),
                'reported_pair_max_abs_fractional_change_vs_level4':
                    float(np.max(np.abs(relative[:, 3, 4]))),
                'reported_pair_at_ell': float(gg[ell_index, 3, 4]),
                'reported_pair_at_ell_fractional_change':
                    float(relative[ell_index, 3, 4]),
                'data_vector_autos_max_abs_fractional_change_vs_level4':
                    float(np.max(np.abs(current['data_vector_autos'] /
                                        reference['data_vector_autos'] - 1))),
            })
        records.append({'case': case, 'w': float(reference['w']), 'levels': levels})
    return records


if __name__ == '__main__':
    root = output_dir()
    (root/'cocoa_convergence.json').write_text(json.dumps(collect(root), indent=2)+'\n')
