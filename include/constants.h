#ifndef CCRH_CONSTANTS_H_
#define CCRH_CONSTANTS_H_

#include <cmath>

// Replacing macros with constexpr functions for type safety and better performance
constexpr double pow2(double A) { return A * A; }
constexpr double pow3(double A) { return A * A * A; }
constexpr double pow4(double A) { return A * A * A * A; }

namespace cgs {

// CGS units
constexpr double centimeter = 1;
constexpr double cm = centimeter;
constexpr double gram = 1;
constexpr double second = 1;
constexpr double s = second;
constexpr double erg = 1;
constexpr double statC = 1;
constexpr double Gauss = 1;
constexpr double Kelvin = 1;

// Derived units
constexpr double K = Kelvin;
constexpr double meter = 1e2 * centimeter;
constexpr double kilometer = 1e3 * meter;
constexpr double km = kilometer;
constexpr double cm2 = cm * cm;
constexpr double cm3 = cm * cm * cm;
constexpr double kilogram = 1e3 * gram;
constexpr double joule = 1e7 * erg;
constexpr double tesla = 1e4 * Gauss;
constexpr double microgauss = 1e-6 * Gauss;
constexpr double nanogauss = 1e-9 * Gauss;
constexpr double muG = microgauss;
constexpr double nG = nanogauss;
constexpr double barn = 1e-24 * pow2(centimeter);
constexpr double mbarn = 1e-27 * pow2(centimeter);

// Electron volt
constexpr double electronvolt = 1.60217657e-12 * erg;
constexpr double kiloelectronvolt = 1e3 * electronvolt;
constexpr double megaelectronvolt = 1e6 * electronvolt;
constexpr double gigaelectronvolt = 1e9 * electronvolt;
constexpr double teraelectronvolt = 1e12 * electronvolt;
constexpr double petaelectronvolt = 1e15 * electronvolt;
constexpr double exaelectronvolt = 1e18 * electronvolt;
constexpr double eV = electronvolt;
constexpr double keV = kiloelectronvolt;
constexpr double MeV = megaelectronvolt;
constexpr double GeV = gigaelectronvolt;
constexpr double TeV = teraelectronvolt;
constexpr double PeV = petaelectronvolt;
constexpr double EeV = exaelectronvolt;

// Time
constexpr double year = 3.15569e7 * second;
constexpr double kiloyear = 1e3 * year;
constexpr double megayear = 1e6 * year;
constexpr double gigayear = 1e9 * year;
constexpr double kyr = kiloyear;
constexpr double Myr = megayear;
constexpr double Gyr = gigayear;

// Parsec
constexpr double parsec = 3.0856775807e18 * centimeter;
constexpr double kiloparsec = 1e3 * parsec;
constexpr double megaparsec = 1e6 * parsec;
constexpr double gigaparsec = 1e9 * parsec;
constexpr double pc = parsec;
constexpr double kpc = kiloparsec;
constexpr double Mpc = megaparsec;
constexpr double Gpc = gigaparsec;

// Physical constants
constexpr double c_light = 2.99792458e10 * centimeter / second;
constexpr double c_light_squared = c_light * c_light;
constexpr double mass_proton = 1.67262158e-24 * gram;
constexpr double mass_proton_c2 = mass_proton * c_light_squared;
constexpr double mass_neutron = 1.67492735e-24 * gram;
constexpr double mass_electron = 9.10938291e-28 * gram;
constexpr double mass_electron_c2 = mass_electron * c_light_squared;
constexpr double electron_radius = 2.8179409238e-13 * cm;
constexpr double electron_charge = 4.80320425e-10 * statC;
constexpr double mass_sun = 1.989e30 * kilogram;
constexpr double h_planck = 6.62606957e-27 * erg * second;
constexpr double h_bar_planck = h_planck / (2. * M_PI);
constexpr double k_boltzmann = 1.3806488e-16 * erg / Kelvin;
constexpr double G_N = 6.67259e-8 * pow3(cm) / gram / pow2(second);
constexpr double sigma_th = 6.6524e-25 * pow2(cm);
constexpr double alpha = 0.0072973525643;

constexpr double ionization_potential_H = 13.6 * eV;
constexpr double ionization_potential_He = 24.6 * eV;

// PLANCK Cosmological constants
constexpr double hlittle = 0.6711;
constexpr double Omega_m = 0.3175;
constexpr double Omega_l = 1.0 - Omega_m;
constexpr double Omega_b = (0.022068 / hlittle) / hlittle;
constexpr double Omega_n = 0.0;
constexpr double Omega_k = 0.0;
constexpr double Omega_r = 8.6e-5;
constexpr double Omega_tot = 1.0;
constexpr double Y_He = 0.247695;
constexpr double T_cmb = 2.728 * K;

// PLANCK derived constants
constexpr double H_0 = hlittle * 3.2407e-18 / s;
constexpr double rho_c = 3.0 * pow2(H_0) / (8.0 * M_PI * G_N);
constexpr double n_H_0 = rho_c * Omega_b * (1.0 - Y_He) / mass_proton;
constexpr double n_He_0 = rho_c * Omega_b * Y_He / (4.0 * mass_proton);
constexpr double f_H = n_H_0 / (n_H_0 + n_He_0);
constexpr double f_He = n_He_0 / (n_H_0 + n_He_0);
constexpr double IGM_Temperature = 1e4 * cgs::K;

// Source constants
// constexpr double reference_energy = 1. * GeV;
// constexpr double SN_efficiency = 0.1;
// constexpr double SN_kinetic_energy = 1e51 * erg;
// constexpr double SN_fraction = 0.01 / mass_sun;
// constexpr double initial_redshift = 20.;
// constexpr double UV_photoionization_cs = 6.3e-18 * pow2(cm);
// constexpr double PopII_spectrum_slope = 5.;
// constexpr double PopII_dNdM = 8e60 / mass_sun;
// constexpr double clumping_factor = 2.;

}  // namespace cgs

#endif /* CCRH_CONSTANTS_H_ */
