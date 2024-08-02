#ifndef CCRH_UTILITIES_H_
#define CCRH_UTILITIES_H_

#include <string>

#include "constants.h"
// #include "cosmo_progs.h"

// double dEdt_ionization(const double& n_HI, const double& E_k);

// double dEdt_coulomb(const double& n_e, const double& E_k);

// double dEdt_pp(const double& n_H, const double& E_k);

// double dEdt_coulomb_Galprop(const double& n_e, const double& E_k);

// double dEdt_adiabatic(const double& z, const double& E_k);

// double fragmentation_timescale(const double& n_H);

// double n_H_physical(const double & z);

// void print_timescales(std::string filename, const double& z);

// double free_fall_timescale(const double& z, const double& M);

// double circular_velocity_at_rvir(const double& M, const double& z);

// double halo_mass_given_vvir(const double& v, const double& z);

// double spectrum(const double& E_k, const double& SN_slope);

// double compute_spectrum_normalization(double E_0, double E_min, double E_max, double alpha);

// double compute_initial_tau(const double& init_redshift);

// double UV_mean_free_path(const double& z);

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
