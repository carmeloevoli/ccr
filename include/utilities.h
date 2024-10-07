#ifndef CCRH_UTILITIES_H_
#define CCRH_UTILITIES_H_

#include <string>

#include "constants.h"

double beta(double E_k);

double lorentz_factor(double E_k);

double dtdz(double z);

double min_star_forming_halo(double z);

double hubble_time(double z);

double larmor_radius(double B, double E_k);

double D_Bohm(double B, double E_k);

double diffusion_time(double B, double d, double E_k);

double n_H_physical(double z);

double inelastic_time(double n_H, double E_k);

double dEdz_H(double z, double E_k);

double dEdt_i(double z, double E_k);

double dEdt_C(double n_e, double E_k);

double dEdt_C_Galprop(double n_e, double E_k);

#endif /* UTILITIES_H_ */
