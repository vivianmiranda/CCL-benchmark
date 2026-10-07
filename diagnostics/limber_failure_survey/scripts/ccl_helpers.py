"""Native CCL kernels and native GSL rules, with all statuses retained."""
import ctypes as ct
from ctypes.util import find_library
import os
from pathlib import Path
import sys
import warnings

import numpy as np
import pyccl as ccl

from common import VENDOR, finite_float
sys.path.insert(0, str(VENDOR))
from validate_limber_bad_point import limber_integrand


def make_cosmology(case):
    p = case['parameters']
    return ccl.Cosmology(Omega_c=p['Omega_m']-p['Omega_b'], Omega_b=p['Omega_b'],
        h=p['h'], n_s=p['n_s'], sigma8=p['sigma8'], w0=p['w0'], wa=p['wa'],
        m_nu=p['m_nu'], mass_split='equal', transfer_function='boltzmann_camb',
        extra_parameters={'camb': {'dark_energy_model': case['camb_dark_energy_model']}})


def make_tracers(cosmo, z, nz, centers):
    bias = 1.05 / ccl.growth_factor(cosmo, 1/(1+centers))
    tracers = [ccl.NumberCountsTracer(cosmo, has_rsd=False, dndz=(z, row),
        bias=(z, np.full_like(z, b)), mag_bias=None) for row, b in zip(nz, bias)]
    return tracers, bias


def warning_records(records):
    return [{'category': r.category.__name__, 'message': str(r.message)} for r in records]


def native_spectra(cosmo, tracers, ell):
    rows = []
    for i in range(5):
        for j in range(i, 5):
            row = {'pair': [i, j]}
            for method in ('qag_quad', 'spline'):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    try:
                        value = ccl.angular_cl(cosmo, tracers[i], tracers[j], ell,
                            l_limber=-1, limber_integration_method=method)
                        row[method] = {'value': finite_float(value)}
                    except Exception as error:
                        row[method] = {'value': None, 'exception': repr(error)}
                row[method]['warnings'] = warning_records(caught)
            rows.append(row)
    return rows


class GSLRules:
    def __init__(self):
        candidates = [Path(sys.base_prefix)/'lib'/n for n in ('libgsl.dylib', 'libgsl.so')]
        self.path = os.environ.get('GSL_LIBRARY') or next(
            (str(p) for p in candidates if p.exists()), find_library('gsl'))
        if not self.path:
            raise RuntimeError('Set GSL_LIBRARY to CCL’s GSL library')
        self.lib = ct.CDLL(self.path)
        self.version = ct.c_char_p.in_dll(self.lib, 'gsl_version').value.decode()
        self.callback = ct.CFUNCTYPE(ct.c_double, ct.c_double, ct.c_void_p)
        class Function(ct.Structure):
            _fields_ = [('function', self.callback), ('params', ct.c_void_p)]
        self.Function = Function
        ptr, dbl, size, fun = ct.c_void_p, ct.c_double, ct.c_size_t, ct.POINTER(Function)
        realptr = ct.POINTER(dbl)
        self.lib.gsl_set_error_handler_off()
        for prefix in ('gsl_integration_workspace', 'gsl_integration_cquad_workspace', 'gsl_integration_glfixed_table'):
            fn = getattr(self.lib, prefix+'_alloc'); fn.argtypes = [size]; fn.restype = ptr
            getattr(self.lib, prefix+'_free').argtypes = [ptr]
        self.lib.gsl_integration_qag.argtypes = [fun,dbl,dbl,dbl,dbl,size,ct.c_int,ptr,realptr,realptr]
        self.lib.gsl_integration_cquad.argtypes = [fun,dbl,dbl,dbl,dbl,ptr,realptr,realptr,ct.POINTER(size)]
        self.lib.gsl_integration_qk41.argtypes = [fun,dbl,dbl,realptr,realptr,realptr,realptr]
        self.lib.gsl_integration_glfixed_point.argtypes = [dbl,dbl,size,realptr,realptr,ptr]

    def context(self, cosmo, tracers, pair, ell):
        first, second = [tracers[i] for i in pair]
        power = cosmo.get_nonlin_power()
        near = max(first._trc[0].chi_min, second._trc[0].chi_min)
        far = min(first._trc[0].chi_max, second._trc[0].chi_max)
        if near <= 0:
            near = 0.5*(ell+0.5)/ccl.spline_params.K_MAX
        lo = float(np.log(max(ccl.spline_params.K_MIN, (ell+0.5)/far)))
        hi = float(np.log(min(ccl.spline_params.K_MAX, 2*(ell+0.5)/near)))
        calls, errors = [], []
        def evaluate(logk, unused):
            calls.append(logk)
            try:
                return limber_integrand(logk, ell, cosmo, first, second, power)/(ell+0.5)
            except Exception as error:
                errors.append(repr(error))
                return float('nan')
        callback = self.callback(evaluate)
        function = self.Function(callback, None)
        return {'function': function, 'callback': callback, 'evaluate': evaluate,
                'lo': lo, 'hi': hi, 'calls': calls, 'errors': errors}

    def integrate(self, context, method, epsrel, parts=1):
        c = context; c['calls'].clear(); c['errors'].clear()
        value, error = ct.c_double(), ct.c_double()
        total = estimated = 0.; statuses = []
        prefix = 'gsl_integration_cquad_workspace' if method == 'cquad' else 'gsl_integration_workspace'
        workspace = getattr(self.lib, prefix+'_alloc')(1000)
        if not workspace:
            raise MemoryError('GSL workspace allocation failed')
        try:
            edges = np.linspace(c['lo'], c['hi'], parts+1)
            for left, right in zip(edges[:-1], edges[1:]):
                if method == 'cquad':
                    count = ct.c_size_t()
                    code = self.lib.gsl_integration_cquad(ct.byref(c['function']), left,right,
                        0.,epsrel,workspace,ct.byref(value),ct.byref(error),ct.byref(count))
                else:
                    code = self.lib.gsl_integration_qag(ct.byref(c['function']),left,right,
                        0.,epsrel,1000,4,workspace,ct.byref(value),ct.byref(error))
                total += value.value; estimated += error.value; statuses.append(code)
        finally:
            getattr(self.lib, prefix+'_free')(workspace)
        return {'value': finite_float(total), 'estimated_error': finite_float(estimated),
                'status': statuses, 'calls': len(c['calls']), 'callback_errors': list(c['errors']),
                'epsrel': epsrel, 'parts': parts}

    def initial_rules(self, context):
        c = context; c['calls'].clear(); c['errors'].clear()
        v,e,absolute,asc = [ct.c_double() for _ in range(4)]
        self.lib.gsl_integration_qk41(ct.byref(c['function']),c['lo'],c['hi'],
            ct.byref(v),ct.byref(e),ct.byref(absolute),ct.byref(asc))
        nodes = np.array(c['calls'])
        table = self.lib.gsl_integration_glfixed_table_alloc(20)
        x,w=ct.c_double(),ct.c_double(); gauss=0.; gn=[]
        try:
            for i in range(20):
                self.lib.gsl_integration_glfixed_point(c['lo'],c['hi'],i,ct.byref(x),ct.byref(w),table)
                gauss += w.value*c['evaluate'](x.value,None); gn.append(x.value)
        finally:
            self.lib.gsl_integration_glfixed_table_free(table)
        return {'kronrod': finite_float(v.value), 'gauss': finite_float(gauss),
                'estimated_error': finite_float(e.value), 'callback_errors': list(c['errors'])}, nodes, np.array(gn)


def reference(rules, context, config):
    result = {'native_qag_replay': rules.integrate(context, 'qag', 1e-4),
              'cquad': rules.integrate(context, 'cquad', config['reference_epsrel']),
              'split2': rules.integrate(context, 'qag', config['reference_epsrel'], 2),
              'split4': rules.integrate(context, 'qag', config['reference_epsrel'], 4)}
    value = result['cquad']['value']
    valid = value is not None and value > 0
    for key in ('cquad', 'split2', 'split4'):
        item = result[key]
        valid &= item['value'] is not None and not any(item['status']) and not item['callback_errors']
    if valid:
        delta = max(abs(result[key]['value']/value-1) for key in ('split2','split4'))
        valid &= delta <= config['reference_consistency_fractional_tolerance']
    else:
        delta = None
    result['reference_valid'] = bool(valid)
    result['maximum_split_reference_fractional_change'] = delta
    result['reference_value'] = value
    if not valid:
        extra = supplement_reference(rules, context, config, result)
        result['supplementary_reference_checks'] = extra
        result['reference_valid'] = extra['reference_valid']
        result['reference_value'] = extra['reference_value']
    return result


def supplement_reference(rules, context, config, original):
    """Add successful refined controls; preserve every failed GSL return."""
    extra = {
        'cquad_refined': rules.integrate(context, 'cquad', 1e-8),
        'cquad_refined_split2': rules.integrate(context, 'cquad', 1e-8, 2),
        'qag_split8': rules.integrate(context, 'qag', config['reference_epsrel'], 8),
        'acceptance_protocol': 'Both refined CQUAD evaluations must succeed and agree; every available successful split-QAG control must agree. Failed controls remain recorded and do not pass.'}
    def successful(record):
        return (record.get('value') is not None and record['value']>0 and
                not any(record['status']) and not record['callback_errors'])
    value = extra['cquad_refined']['value']
    valid = successful(extra['cquad_refined']) and successful(extra['cquad_refined_split2'])
    comparisons = {};failed = [];successful_qag=[]
    if valid:
        comparisons['cquad_refined_split2'] = abs(extra['cquad_refined_split2']['value']/value-1)
        if successful(original['cquad']):
            comparisons['original_cquad'] = abs(original['cquad']['value']/value-1)
        for key,record in [('original_split2',original['split2']),('original_split4',original['split4']),('qag_split8',extra['qag_split8'])]:
            if successful(record):
                successful_qag.append(key)
                comparisons[key] = abs(record['value']/value-1)
            else:failed.append({'key':key,'status':record['status'],'callback_errors':record['callback_errors']})
        valid = max(comparisons.values())<=config['reference_consistency_fractional_tolerance']
    extra.update(reference_valid=bool(valid),reference_value=value,
        successful_split_QAG_controls=successful_qag,failed_split_QAG_controls=failed,
        fractional_reference_comparisons=comparisons,
        acceptance_tolerance=config['reference_consistency_fractional_tolerance'])
    return extra
