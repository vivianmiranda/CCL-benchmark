"""Check actual CoCoA fixed-rule spectra in a fresh initialized process."""
import argparse
import json
from pathlib import Path
import os
import sys
import numpy as np
import cosmolike_lsst_y1_interface as ci
sys.path.insert(0, str(Path(os.environ['ROOTDIR'])/'external_modules/code/cosmolike_core'))
from cosmolike_notebook_utils.covariance.power import refine_power_tables

parser = argparse.ArgumentParser()
parser.add_argument('case', type=int)
parser.add_argument('level', type=int, choices=range(5))
parser.add_argument('--save-nodes', action='store_true')
args = parser.parse_args()
from paths import output_dir, INPUTS, VENDOR
root = output_dir()
bundle = np.load(root/f'shared_input_{args.case}.npz')
names = ('log10k_2D', 'z_2D', 'lnP_linear', 'lnP_nonlinear',
         'lnP_linear_cb', 'G', 'z_G', 'z_1D', 'chi', 'omegan2')
tables = {name: bundle[name] for name in names}
tables['omegan2'] = float(tables['omegan2'])
tables = refine_power_tables(tables, refinement=8)
ci.initial_setup()
ci.init_accuracy_boost(accuracy_boost=1., integration_accuracy=args.level)
ci.init_probes(possible_probes='3x2pt')
ci.init_IA(ia_model=0, ia_redshift_evolution=2, ia_code=0)
ci.init_bias(bias_model=[0, 0, 0, 0, 0])
ci.init_photoz_conventions(interpolation_type=1, zmid_convention=1)
ci.init_cosmo_runmode(is_linear=False)
ci.init_redshift_distributions_from_files(
    lens_multihisto_file=str(root/'frozen_nz.txt'), lens_ntomo=5,
    source_multihisto_file=str(root/'frozen_nz.txt'), source_ntomo=5)
ci.set_cosmology(omegam=.3156, omegab=.0492, H0=67.27, **tables)
zeros = [0.]*5
ci.set_nuisance_bias(B1=bundle['bias'].tolist(), B2=zeros, B_MAG=zeros,
                     B3nl=zeros, BK=zeros)
ci.set_nuisance_ia(A1=zeros, A2=zeros, B_TA=zeros)
ci.set_nuisance_shear_photoz(bias=zeros)
ci.set_nuisance_clustering_photoz(bias=zeros)
ci.set_nuisance_shear_calib(M=zeros)
ell = np.ascontiguousarray(bundle['ell'])
# Production LSST panel boundaries, with its endpoint covering this frozen n(z).
a_edges = 1./(1.+np.array([3.5, 2., 1.5, 1., .7, .4, .2, 1.e-5]))
nodes = (96, 128, 256, 512, 1024)[args.level]
cov = ci.covariance.covariance_spectra(
    ell=ell, a_edges=a_edges, nquad=nodes, nwindow=4097,
    include_ia=False, include_rsd=False, linear=False, nonlimber_lmax=0)
# This notebook exposure calls the real data-vector batch C calculation;
# only the diagonal galaxy-bin entries are supported on this path.
dv = ci.C_gg_tomo_limber(ell)
output = root/f'cocoa_case{args.case}_level{args.level}.npz'
payload = dict(ell=ell, covariance_spectra=cov['spectra'],
         data_vector_autos=np.diagonal(dv, axis1=1, axis2=2),
         w=bundle['w'], a_edges=a_edges, nodes=np.array(nodes))
np.savez(output, **payload)
if args.save_nodes:
    np.savez(root/f'cocoa_case{args.case}_level{args.level}_nodes.npz',
             **payload, geometry=cov['geometry'])
if not np.all(np.isfinite(cov['spectra'])) or not np.all(np.isfinite(dv)):
    raise ValueError('non-finite CoCoA spectra')
ibad = np.flatnonzero(ell == 355.655882)[0]
print(json.dumps(dict(case=args.case, level=args.level,
    w=float(bundle['w']), cov_nodes=nodes,
    bad_pair_value=float(cov['spectra'][ibad,3,4]), output=str(output))), flush=True)
