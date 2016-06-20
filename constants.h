#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <cmath>
#include <iostream>

using namespace std;

// pows
#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define pow4(A) ((A)*(A)*(A)*(A))

// CGS units
static const double centimeter = 1;
static const double cm = centimeter;
static const double gram = 1;
static const double second = 1;
static const double s = second;
static const double dyne = 1;
static const double erg = 1;
static const double esu = 1;
static const double gauss = 1;
static const double kelvin = 1;
static const double K = kelvin;

// derived units
static const double meter = 1e2 * centimeter;
static const double kilometer = 1e3 * meter;
static const double km = kilometer;
static const double cm2 = cm * cm;
static const double cm3 = cm * cm * cm;
static const double kilogram = 1e3 * gram;
static const double newton = 1e5 * dyne;
static const double joule = 1e7 * erg;
static const double coulomb = 2997924580. * esu;
static const double tesla = 1e4 * gauss;
static const double microgauss = 1e-6 * gauss;
static const double nanogauss = 1e-9 * gauss;
static const double muG = microgauss;
static const double nG = nanogauss;

// electron volt
static const double electronvolt = 1.60217657e-12 * erg;
static const double kiloelectronvolt = 1e3 * electronvolt;
static const double megaelectronvolt = 1e6 * electronvolt;
static const double gigaelectronvolt = 1e9 * electronvolt;
static const double teraelectronvolt = 1e12 * electronvolt;
static const double petaelectronvolt = 1e15 * electronvolt;
static const double exaelectronvolt = 1e18 * electronvolt;
static const double eV = electronvolt;
static const double keV = kiloelectronvolt;
static const double MeV = megaelectronvolt;
static const double GeV = gigaelectronvolt;
static const double TeV = teraelectronvolt;
static const double PeV = petaelectronvolt;
static const double EeV = exaelectronvolt;

// time
static const double year = 3.15569e7 * second;
static const double kiloyear = 1e3 * year;
static const double megayear = 1e6 * year;
static const double gigayear = 1e9 * year;
static const double kyr = kiloyear;
static const double Myr = megayear;
static const double Gyr = gigayear;

// parsec
static const double parsec = 3.0856775807e18 * centimeter;
static const double kiloparsec = 1e3 * parsec;
static const double megaparsec = 1e6 * parsec;
static const double gigaparsec = 1e9 * parsec;
static const double pc = parsec;
static const double kpc = kiloparsec;
static const double Mpc = megaparsec;
static const double Gpc = gigaparsec;

// physical constants
static const double c_light = 2.99792458e10 * centimeter / second;
static const double c_light_squared = c_light * c_light;
static const double mass_proton = 1.67262158e-24 * gram;
static const double mass_proton_c2 = mass_proton * c_light_squared;
static const double mass_neutron = 1.67492735e-24 * gram;
static const double mass_electron = 9.10938291e-28 * gram;
static const double mass_electron_c2 = mass_electron * c_light_squared;
static const double electron_radius = 2.8179409238e-13 * cm; /* e^2 / mc^2 classical electron radius. */
static const double mass_sun = 1.989e30 * kilogram;
static const double h_planck = 6.62606957e-27 * erg * second;
static const double h_bar_planck = h_planck / 2. / M_PI;
static const double k_boltzmann = 1.3806488e-16 * erg / kelvin;
static const double G = 6.67259e-8 * pow3(cm) / gram / pow2(second);
static const double pi_re2_me_c2_c = M_PI * pow2(electron_radius) * c_light * mass_electron_c2; /* Pi * e^4 / mc */
static const double pi_re_h_bar2_c2 = M_PI * electron_radius * pow2(h_bar_planck) * pow2(c_light); // erg^2* cm^3

static const double A_H = 1;
static const double Z_H = 1;
static const double ionization_potential_H = 13.6 * eV;
static const double A_He = 4;
static const double Z_He = 2;
static const double ionization_potential_He = 24.6 * eV;

// PLANCK Cosmological constants
static const double hlittle = 0.6711;
static const double Omega_m = 0.3175;
static const double Omega_l = 1.0 - Omega_m;
static const double Omega_b = (0.022068 / hlittle) / hlittle;
static const double Omega_n = 0.0;
static const double Omega_k = 0.0;
static const double Omega_r = 8.6e-5;
static const double Omega_tot = 1.0;
static const double Y_He = 0.247695;
static const double T_cmb = 2.728 * K;

// PLANCK derived
static const double H_0 = hlittle * 3.2407e-18 / s;
static const double rho_crit = 3.0 * H_0 * H_0 / (8.0 * M_PI * G);  /* at z = 0 */
static const double n_H_0 = rho_crit * Omega_b * (1.0 - Y_He) / mass_proton; /*  current hydrogen number density estimate ~1.92e-7 */
static const double n_He_0 = rho_crit * Omega_b * Y_He / (4.0 * mass_proton); /*  current helium number density estimate */
static const double f_H = n_H_0 / (n_H_0 + n_He_0); /* hydrogen number fraction */
static const double f_He = n_He_0 / (n_H_0 + n_He_0); /* helium number fraction */

// source constants
static const double reference_energy = 1. * GeV;
static const double SN_slope = 2.2;
static const double SN_efficiency = 0.1;
static const double SN_kinetic_energy = 1e51 * erg;
static const double SN_fraction = 0.01 / mass_sun;
static const double SN_E_min = 0.1 * GeV;
static const double SN_E_max = 1e6 * GeV;
static const double initial_redshift = 30.;

#include "cosmology_params.h"
#include "anal_params.h"
#include "heat_params.h"

#endif /* CONSTANTS_H_ */
