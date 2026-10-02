"""Diagnostic (Section 4 of the README): what makes CAMB's linear growth
depend on k at w = -0.9. Prints D(k,z)/D(k0,z) - 1 at z = 1, with
D(k,z) = sqrt(P_lin(k,z)/P_lin(k,0)) and k0 = 5e-4 /Mpc (where CoCoA's
likelihood measures its growth factor), with and without CAMB's
dark-energy perturbations and massive neutrinos. The fiducial parameters
of the study; run in the Cocoa environment (CAMB).
python diag_de_perturbations.py"""
import numpy as np
import camb

KS = np.array([5e-4, 1e-3, 3e-3, 0.01, 0.05, 0.2])
H0, OMB, OMM = 67.32, 0.04, 0.3

def growth(w, mnu, perturbations, z=1.0):
  h = H0/100.0
  p = camb.CAMBparams()
  p.set_cosmology(H0=H0, ombh2=OMB*h**2, omch2=(OMM - OMB)*h**2 - mnu/93.14, mnu=mnu,
                  num_massive_neutrinos=1 if mnu > 0 else 0)
  p.InitPower.set_params(As=2.1e-9, ns=0.96605)
  p.set_dark_energy(w=w, dark_energy_model="fluid")
  # the ctypes field is "___no_perturbations"; setting "no_perturbations"
  # creates a Python attribute and changes nothing
  setattr(p.DarkEnergy, "___no_perturbations", not perturbations)
  p.set_matter_power(redshifts=[z, 0.0], kmax=2.0)
  p.NonLinear = camb.model.NonLinear_none
  PK = camb.get_results(p).get_matter_power_interpolator(nonlinear=False, hubble_units=False, k_hunit=False)
  D = np.sqrt(PK.P(z, KS)/PK.P(0.0, KS))
  return 100*(D/D[0] - 1)

print("D(k,z=1)/D(k0,z=1) - 1 [%%], k = %s /Mpc" % ", ".join("%g" % k for k in KS))
for w, mnu, pert in ((-0.9, 0.0, True), (-0.9, 0.0, False), (-0.9, 0.06, True), (-0.9, 0.06, False),
                     (-1.0, 0.0, True), (-1.0, 0.06, True)):
  print("w = %4.1f, mnu = %.2f eV, dark-energy perturbations %-3s: %s"
        % (w, mnu, "on" if pert else "off", " ".join("%+.3f" % x for x in growth(w, mnu, pert))))
