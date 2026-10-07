"""Locate accidental agreement of the initial 41/20 rules in the frozen case."""
import ctypes as ct
import json
import numpy as np
import pyccl as ccl
from scipy.optimize import brentq
from probe_quadrature import (gsl, real, function_pointer, real_pointer,
    callback_type, Function, diagnostic, data, z, bins, root)
gsl.gsl_integration_qk41.argtypes = [function_pointer,real,real,real_pointer,
    real_pointer,real_pointer,real_pointer]
records=[]
def evaluate(w):
    """Return the signed disagreement before GSL's error rescaling."""
    c=diagnostic.make_cosmology(w)
    t=diagnostic.make_tracers(c,z,bins,data['bin_centers'])
    ell=355.655882
    p=c.get_nonlin_power()
    lo=np.log((ell+.5)/min(t[3]._trc[0].chi_max,t[4]._trc[0].chi_max))
    hi=np.log(2*(ell+.5)/max(t[3]._trc[0].chi_min,t[4]._trc[0].chi_min))
    def integrand(x,unused):
        return diagnostic.limber_integrand(x,ell,c,t[3],t[4],p)/(ell+.5)
    callback=callback_type(integrand); f=Function(callback,None)
    value,error,resabs,resasc=real(),real(),real(),real()
    gsl.gsl_integration_qk41(ct.byref(f),lo,hi,ct.byref(value),ct.byref(error),ct.byref(resabs),ct.byref(resasc))
    rule=gsl.gsl_integration_glfixed_table_alloc(20)
    x,weight=real(),real(); gauss=0.
    for i in range(20):
        gsl.gsl_integration_glfixed_point(lo,hi,i,ct.byref(x),ct.byref(weight),rule)
        gauss+=weight.value*integrand(x.value,None)
    gsl.gsl_integration_glfixed_table_free(rule)
    record=dict(w=float(w),kronrod=value.value,gauss=gauss,
                signed_difference=value.value-gauss,estimated_error=error.value)
    records.append(record)
    print(json.dumps(record),flush=True)
    return value.value-gauss
previous=None
roots=[]
for w in np.linspace(-1.10,-.9,21):
    value=evaluate(w)
    if previous is not None and value*previous[1]<0:
        solution=brentq(evaluate,previous[0],w,xtol=1.e-10)
        roots.append(solution)
    previous=(w,value)
(root/'first_panel_roots.json').write_text(json.dumps(dict(roots=roots,records=records),indent=2)+'\n')
print('ROOTS',roots,flush=True)
