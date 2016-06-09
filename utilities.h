#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "constants.h"
#include "cosmo_progs.h"

using namespace std;

double dEdt_ionization(const double& n_H, const double& E_k);

double dEdt_coulomb(const double& n_e, const double& E_k);

double dEdt_coulomb_Galprop(const double& n_e, const double& E_k);

double dEdt_adiabatic(const double& z, const double& E_k);

double fragmentation_timescale(const double& n_H);

double n_H_physical(const double & z);

void print_timescales(string filename, const double& z);

double min_star_forming_halo(const double& z);

double free_fall_timescale(const double& z, const double& M);

#endif /* UTILITIES_H_ */
