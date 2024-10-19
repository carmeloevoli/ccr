#ifndef CCRH_UTILITIES_H_
#define CCRH_UTILITIES_H_

#include "constants.h"

double beta_lorentz(double E_k);

double gamma_lorentz(double E_k);

double dtdz(double z);

double min_star_forming_halo(double z);

double hubble_time(double z);

double larmor_radius(double B, double E_k);

double D_Bohm(double B, double E_k);

double diffusion_time(double B, double d, double E_k);

double n_H_physical(double z);

double inelastic_time(double n_H, double E_k);

double dEdt_a(double z, double E_k);

double dEdt_i(double n_HI, double E_k);

double dEdt_C(double n_e, double E_k);

double SN_Spectrum(double Ek, double alpha);

double SN_ESpectrum_integral(double alpha);

#endif /* UTILITIES_H_ */
