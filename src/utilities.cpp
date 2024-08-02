#include "utilities.h"

#include <iostream>

double lorentz_factor(double E_k) { return 1. + E_k / cgs::mass_proton_c2; }

double beta(double E_k) { return std::sqrt(1.0 - 1.0 / pow2(1.0 + E_k / cgs::mass_proton_c2)); }

double larmor_radius(double B, double E_k) {
    const auto p = E_k * sqrt(1. + 2. * cgs::mass_proton_c2 / E_k);
    const auto Z = 1.;
    return p / Z / cgs::electron_charge / B;
}

double D_Bohm(double B, double E_k) {
    const auto rL = larmor_radius(B, E_k);
    return rL * cgs::c_light / 3.;
}

double diffusion_time(double B, double d, double E_k) {
    const auto DB = D_Bohm(B, E_k);
    return d * d / DB;
}

double hubble_time(double z) {
    return 2. / 3. / cgs::H_0 / std::sqrt(cgs::Omega_m) * std::pow(1. + z, -1.5);
}

double n_H_physical(double z) { return cgs::n_H_0 * std::pow(1. + z, 3.); }

double min_star_forming_halo(double z) {
    return 1e8 * cgs::mass_sun * std::pow(10. / (1. + z), 1.5);
}

double sigma_pp(double E_k) {  // Kafexhiu et al. PRD 90, 2014
    const auto E_threshold = 0.2797 * cgs::GeV;
    const auto x = E_k / E_threshold;
    double value = 0;
    if (x > 1) {
        value = 30.7 - 0.96 * log(x) + 0.18 * pow2(log(x));
        value *= pow3(1 - pow(x, -1.9));
    }
    return value * cgs::mbarn;
}

double inelastic_time(double n_H, double E_k) {
    const auto sigma = sigma_pp(E_k);
    return 1. / beta(E_k) / cgs::c_light / n_H / sigma;
}

double dEdz_H(double z, double E_k) {
    auto value = E_k / (1. + z);
    value *= 2. - 1. / (1. + cgs::mass_proton_c2 / E_k);
    return value;
}

/* function DTDZ returns the value of dt/dz at the redshift parameter z. */
double dtdz(double z) {
    using std::pow;
    using std::sqrt;

    const auto x = sqrt(cgs::Omega_l / cgs::Omega_m) * pow(1 + z, -3.0 / 2.0);
    const auto dxdz = sqrt(cgs::Omega_l / cgs::Omega_m) * pow(1 + z, -5.0 / 2.0) * (-3.0 / 2.0);
    const auto const1 = 2 * sqrt(1 + cgs::Omega_m / cgs::Omega_l) / (3.0 * cgs::H_0);

    const auto numer = dxdz * (1 + x * pow(pow(x, 2) + 1, -0.5));
    const auto denom = x + sqrt(pow(x, 2) + 1);
    return (const1 * numer / denom);
}

double dEdt_i(double n_HI, double E_k) {
    const auto p2 = pow2(lorentz_factor(E_k)) - 1.;
    const auto beta2 = pow2(beta(E_k));
    const auto factor = 4. * cgs::pi_re2_me_c2_c / beta(E_k);
    const auto A_H = cgs::Z_H * n_HI *
                     (log(2. * cgs::mass_electron_c2 * p2 / cgs::ionization_potential_H) - beta2);
    const auto A_He = cgs::Z_He * cgs::f_He * n_HI *
                      (log(2. * cgs::mass_electron_c2 * p2 / cgs::ionization_potential_He) - beta2);
    return factor * (A_H + A_He);
}

double dEdt_C(double n_e, double E_k) {
    const auto p = sqrt(pow2(lorentz_factor(E_k)) - 1.);
    const auto w_pl = 5.64e4 * sqrt(n_e) / cgs::s;
    const auto Coulomb_log =
        log(2. * cgs::mass_electron_c2 * beta(E_k) / cgs::h_bar_planck / w_pl * p);
    return 4. * cgs::pi_re2_me_c2_c / beta(E_k) * n_e * (Coulomb_log - pow2(beta(E_k)) / 2.);
}

double dEdt_C_Galprop(double n_e, double E_k) {
    const auto beta_2 = pow2(beta(E_k));
    const auto beta_3 = pow3(beta(E_k));
    const auto temperature_electron = 1e4 * cgs::K;
    const auto x_m = pow(3. * sqrt(M_PI) / 4., 1. / 3.) *
                     sqrt(2. * cgs::k_boltzmann * temperature_electron / cgs::mass_electron_c2);
    const auto x_m_3 = pow3(x_m);
    const auto log_1 = pow2(cgs::mass_electron_c2) / cgs::pi_re_h_bar2_c2 / n_e;
    const auto log_2 = cgs::mass_proton_c2 * pow2(lorentz_factor(E_k)) * pow4(beta(E_k)) /
                       (cgs::mass_proton_c2 + 2. * lorentz_factor(E_k) * cgs::mass_electron_c2);
    const auto Coulomb_log = 0.5 * log(log_1 * log_2);
    return 4. * cgs::pi_re2_me_c2_c * n_e * Coulomb_log * (beta_2 / (x_m_3 + beta_3));
}

// double fragmentation_timescale(const double& n_H) {
// 	return 6e7 * year / n_H; // arXiv:1208.4979
// }

// double free_fall_timescale(const double& z, const double& halo_mass) {
// 	double virial_radius_comoving = fast::MtoRvir(z, halo_mass / mass_sun) * Mpc; // M in M_sun,
// Rvir in comoving Mpc 	double virial_radius_physical = virial_radius_comoving / (1. + z);
// double virial_volume_phyisical = 4./3. * M_PI * pow3(virial_radius_physical); 	double
// halo_dark_density = halo_mass / virial_volume_phyisical; 	return sqrt(3. * M_PI / 32. / G /
// halo_dark_density);
// }

// double UV_mean_free_path(const double& z) {
// 	//double measurement = 3.3;
// 	//double lambda_ll = c_light / fast::hubble(z) / (1. + z) / measurement;
// 	//return lambda_ll / sqrt(M_PI);
// 	return 3.9 * Mpc * pow((1. + z) / 4., -5.0); // 4.5
// }

// double circular_velocity_at_rvir(const double& M, const double& z) {
// 	return 24.0 * km / s * pow(M / (1e8 * mass_sun), 1./3.) * pow((1. + z) / 10., 0.5);
// }

// double halo_mass_given_vvir(const double& v, const double& z) {
//     double v_ = v / (24.0 * km / s);
//     double z_ = pow((1. + z) / 10., -0.5);
// 	return 1e8 * mass_sun * pow(v_ * z_, 3.);
// }

// void print_timescales(string filename, const double& z) {
// 	cout << "Print timescales in " << filename << "\n";

// 	const double n_H = n_H_physical(z);
// 	const double n_He = f_He * n_H;
// 	const double ionization_fraction = 1e-4;
// 	const double n_e = ionization_fraction * n_H + 2. * ionization_fraction * n_He;
// 	const double n_HI = (1. - ionization_fraction) * n_H;

// 	cout << scientific << z << "\t" << n_H << "\t" << n_HI << "\t" << n_He << "\t" <<
// ionization_fraction << "\t" << n_e << "\n";

// 	ofstream outfile;
// 	outfile.open(filename.c_str());
// 	outfile << "#E_k [GeV] - t_Hubble - t_ionization - t_Coulomb - t_adiabatic - t_fragmentation
// " << endl; 	outfile << scientific << setprecision(3); 	for (double E_k = MeV; E_k < 200. *
// GeV; E_k
// *= 1.2) { 		outfile << E_k / GeV << "\t" << fast::t_hubble(z) << "\t";
// 		//outfile << spectrum(E_k) << "\t";
// 		outfile << E_k / dEdt_ionization(n_HI, E_k) << "\t";
// 		outfile << E_k / dEdt_coulomb_Galprop(n_e, E_k) << "\t" << E_k / dEdt_coulomb(n_e,
// E_k) << "\t"; 		outfile << E_k / dEdt_adiabatic(z, E_k) << "\t" ;
// outfile << fragmentation_timescale(n_H) << "\n";
// 	}
// 	outfile.close();
// }

// double spectrum(const double& E_k, const double& SN_slope) {
//     double p = sqrt(E_k * E_k + 2. * mass_proton_c2 * E_k);
//     double beta = p / (E_k + mass_proton_c2);
//     double E_k_0 = reference_energy;
//     double p_0 = sqrt(E_k_0 * E_k_0 + 2. * mass_proton_c2 * E_k_0);

//     return 1. / beta * pow(p / p_0, -SN_slope);
// }

// double I_spectrum(double x, void * params) {
// 	double alpha = *(double *) params;
// 	double E_0 = *((double *)params + 1);
//     double E = x * E_0;
// 	//double f = x * pow((E_0 * x * x + 2. * mass_proton_c2 * x) / (E_0 + 2. * mass_proton_c2),
// - .5 * alpha); 	return x * spectrum(E, alpha);
// }

// double compute_spectrum_normalization(double E_0, double E_min, double E_max, double alpha) {
//     gsl_integration_workspace * w
// 	= gsl_integration_workspace_alloc (10000);

// 	double result, error;
// 	double params[2] = {alpha, E_0};
// 	gsl_function F;
// 	F.function = &I_spectrum;
// 	F.params = params;

// 	gsl_integration_qag (&F, E_min / E_0, E_max / E_0, 0, 1e-5, 10000, 3, w, &result, &error);

// 	gsl_integration_workspace_free (w);

// 	return result;
// }

// double compute_initial_tau(const double& init_redshift){
// 	double result = 0;
// 	double dz = 0.1;
// 	for (double z = 1000; z > init_redshift; z -= dz) {
// 		double dt = -fast::dtdz(z) * dz;
// 		double n_e = 1e-4 * n_H_physical(z);
// 		result += SIGMAT * n_e * c_light * dt;
// 	}
// 	return result;
// }
