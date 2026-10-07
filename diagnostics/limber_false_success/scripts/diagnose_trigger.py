"""Separate the smooth physical change from a discontinuous QAG decision.

Frozen-power and frozen-geometry evaluations are diagnostic factor swaps,
not physical cosmologies or changes to the installed CCL implementation.
"""
import ctypes as ct
import json
import numpy as np
import pyccl as ccl
from probe_quadrature import (gsl, real, function_pointer, real_pointer,
    callback_type, Function, diagnostic, data, z, bins, root)

gsl.gsl_integration_qk41.argtypes = [function_pointer,real,real,real_pointer,
    real_pointer,real_pointer,real_pointer]
ell = 355.655882
wroot = json.loads((root/'root_spike.json').read_text())['w']

def prepare(w):
    """Prepare native CCL tracer kernels, power and integration boundaries."""
    cosmo = diagnostic.make_cosmology(w)
    tracers = diagnostic.make_tracers(cosmo,z,bins,data['bin_centers'])
    power = cosmo.get_nonlin_power()
    lo = np.log((ell+.5)/min(tracers[3]._trc[0].chi_max,tracers[4]._trc[0].chi_max))
    hi = np.log(2*(ell+.5)/max(tracers[3]._trc[0].chi_min,tracers[4]._trc[0].chi_min))
    return cosmo, tracers, power, float(lo), float(hi)

def rules(geometry, power, endpoints=None, extra=False):
    """Evaluate the embedded pair, native QAG and optionally CQUAD."""
    c,t,_,lo,hi = geometry
    if endpoints is not None:
        lo,hi = endpoints
    calls=[]
    def integrand(x,unused):
        calls.append(x)
        return diagnostic.limber_integrand(x,ell,c,t[3],t[4],power)/(ell+.5)
    callback=callback_type(integrand)
    f=Function(callback,None)
    value,error,resabs,resasc = real(),real(),real(),real()
    gsl.gsl_integration_qk41(ct.byref(f),lo,hi,ct.byref(value),ct.byref(error),ct.byref(resabs),ct.byref(resasc))
    kronrod_nodes=np.array(calls)
    rule=gsl.gsl_integration_glfixed_table_alloc(20)
    x,weight=real(),real(); gauss=0.
    gauss_nodes=[]
    for i in range(20):
        gsl.gsl_integration_glfixed_point(lo,hi,i,ct.byref(x),ct.byref(weight),rule)
        gauss+=weight.value*integrand(x.value,None)
        gauss_nodes.append(x.value)
    gsl.gsl_integration_glfixed_table_free(rule)
    record=dict(kronrod=value.value,gauss=gauss,estimated_error=error.value,
                tolerance=1.e-4*abs(value.value),bounds=[lo,hi])
    work=gsl.gsl_integration_workspace_alloc(1000)
    calls.clear()
    code=gsl.gsl_integration_qag(ct.byref(f),lo,hi,0.,1.e-4,1000,4,work,
                                ct.byref(value),ct.byref(error))
    gsl.gsl_integration_workspace_free(work)
    record.update(qag=value.value,qag_status=code,qag_calls=len(calls))
    if extra:
        work=gsl.gsl_integration_cquad_workspace_alloc(1000)
        count=ct.c_size_t()
        code=gsl.gsl_integration_cquad(ct.byref(f),lo,hi,0.,1.e-9,work,
                         ct.byref(value),ct.byref(error),ct.byref(count))
        gsl.gsl_integration_cquad_workspace_free(work)
        record.update(cquad=value.value,cquad_error=error.value,cquad_status=code)
    return record,kronrod_nodes,np.array(gauss_nodes)

anchor=prepare(wroot)
records=[]
for offset in np.linspace(-1.e-4,1.e-4,41):
    w=wroot+offset
    current=prepare(w)
    c,t,p,lo,hi=current
    native={method:float(ccl.angular_cl(c,t[3],t[4],ell,l_limber=-1,
            limber_integration_method=method)) for method in ('qag_quad','spline')}
    full,kn,gn=rules(current,p,extra=abs(offset)<1.e-12)
    record=dict(w=w,offset=offset,native=native,full=full)
    if abs(offset)<1.e-12 or abs(abs(offset)-1.e-4)<1.e-12:
        record['frozen_power']=rules(current,anchor[2])[0]
        record['frozen_geometry']=rules(anchor,p)[0]
        record['fixed_root_interval']=rules(current,p,endpoints=anchor[3:])[0]
    records.append(record)
    (root/'trigger_scan.json').write_text(json.dumps(records,indent=2)+'\n')
    print(f'{offset:+.2e}: QAG {full["qag_calls"]} calls, native residual {native["qag_quad"]/native["spline"]-1:.7g}',flush=True)
    if abs(offset)<1.e-12:
        x=np.linspace(lo,hi,10001)
        y=np.array([diagnostic.limber_integrand(xx,ell,c,t[3],t[4],p)/(ell+.5) for xx in x])
        node_y=np.array([diagnostic.limber_integrand(xx,ell,c,t[3],t[4],p)/(ell+.5) for xx in kn])
        node_z=1./ccl.scale_factor_of_chi(c,(ell+.5)/np.exp(kn))-1.
        np.savez(root/'trigger_nodes.npz',x=x,y=y,kronrod_x=kn,kronrod_y=node_y,gauss_x=gn,kronrod_z=node_z)

# Moving an otherwise unused endpoint changes node positions while preserving
# the cosmology and integrand. CQUAD checks the small tail-integral change.
endpoint=[]
for shift in (-.01, -.001, .001, .01):
    r,_,_=rules(anchor,anchor[2],endpoints=(anchor[3],anchor[4]+shift),extra=True)
    endpoint.append(dict(upper_logk_shift=shift,**r))
(root/'endpoint_controls.json').write_text(json.dumps(endpoint,indent=2)+'\n')
