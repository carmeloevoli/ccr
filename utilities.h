#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "constants.h"
#include "cosmo_progs.h"

using namespace std;

double dEdt_ionization(const double& n_HI, const double& E_k);

double dEdt_coulomb(const double& n_e, const double& E_k);

double dEdt_pp(const double& n_H, const double& E_k);

double dEdt_coulomb_Galprop(const double& n_e, const double& E_k);

double dEdt_adiabatic(const double& z, const double& E_k);

double fragmentation_timescale(const double& n_H);

double n_H_physical(const double & z);

void print_timescales(string filename, const double& z);

double min_star_forming_halo(const double& z);

double free_fall_timescale(const double& z, const double& M);

double spectrum(const double& E_k);

double compute_spectrum_normalization(double E_0, double E_min, double E_max, double alpha);

double compute_initial_tau(const double& init_redshift);

double UV_mean_free_path(const double& z);

double baryon_number(const double& z);

double hubble_time(const double& z);

#endif /* UTILITIES_H_ */
