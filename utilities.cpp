#include "utilities.h"

#define gamma(E) (1.0 + E / mass_proton_c2)
#define beta(E) (sqrt(1.0 - 1.0 / pow2(1.0 + E / mass_proton_c2)))

double dEdt_ionization(const double& n_HI, const double& E_k) {
	double p_squared = pow2(gamma(E_k)) - 1.;
	double beta_squared = pow2(beta(E_k));
	double factor = 4. * pi_re2_me_c2_c / beta(E_k);
	double A_H = Z_H * n_HI * (log(2. * mass_electron_c2 * p_squared / ionization_potential_H) - beta_squared);
	double A_He = Z_He * f_He * n_HI * (log(2. * mass_electron_c2 * p_squared / ionization_potential_He) - beta_squared);
	return factor * (A_H + A_He);
}

double dEdt_coulomb(const double& n_e, const double& E_k) {
	double p = sqrt(pow2(gamma(E_k)) - 1.);
	double w_pl = 5.64e4 * sqrt(n_e) / s;
	double Coulomb_log = log(2. * mass_electron_c2 * beta(E_k) / h_bar_planck / w_pl * p);
	return 4. * pi_re2_me_c2_c / beta(E_k) * n_e * (Coulomb_log - pow2(beta(E_k)) / 2.);
}

double dEdt_coulomb_Galprop(const double& n_e, const double& E_k) {
	double beta_2 = pow2(beta(E_k));
	double beta_3 = pow3(beta(E_k));
	double temperature_electron = 1e4 * K;
	double x_m = pow(3. * sqrt(M_PI) / 4., 1./3.) * sqrt(2. * k_boltzmann * temperature_electron / mass_electron_c2);
	double x_m_3 = pow3(x_m);
	double log_1 = pow2(mass_electron_c2) / pi_re_h_bar2_c2 / n_e;
	double log_2 = mass_proton_c2 * pow2(gamma(E_k)) * pow4(beta(E_k)) / (mass_proton_c2 + 2. * gamma(E_k) * mass_electron_c2);
	double Coulomb_log = 0.5 * log(log_1 * log_2);
	return 4. * pi_re2_me_c2_c * n_e * Coulomb_log * (beta_2 / (x_m_3 + beta_3));
}

double dEdt_pp(const double& n_H, const double& E_k) {
	double n = n_H * (1. + f_He);
	return 3.85e-16 * GeV / s * (n * pow3(cm)) * pow(E_k / GeV, 1.28) * pow(E_k / GeV + 200., -0.2);
}

double dEdt_adiabatic(const double& z, const double& E_k) {
	const double dzdt = 1. / fast::dtdz(z);
	return -dzdt / (1. + z) * E_k;
}

double fragmentation_timescale(const double& n_H) {
	return 6e7 * year / n_H; // arXiv:1208.4979
}

double hubble_time(const double& z) {
	return 2. / 3. / H_0 / sqrt(Omega_m) * pow(1. + z, -1.5);
}

double free_fall_timescale(const double& z, const double& halo_mass) {
	double virial_radius_comoving = fast::MtoRvir(z, halo_mass / mass_sun) * Mpc; // M in M_sun, Rvir in comoving Mpc
	double virial_radius_physical = virial_radius_comoving / (1. + z);
	double virial_volume_phyisical = 4./3. * M_PI * pow3(virial_radius_physical);
	double halo_dark_density = halo_mass / virial_volume_phyisical;
	return sqrt(3. * M_PI / 32. / G / halo_dark_density);
}

/* returns the proper mean baryonic density at z in cm^-3 */
double n_H_physical(const double & z) {
	return n_H_0 * pow(1. + z, 3.);
}

double UV_mean_free_path(const double& z) {
	//double measurement = 3.3;
	//double lambda_ll = c_light / fast::hubble(z) / (1. + z) / measurement;
	//return lambda_ll / sqrt(M_PI);
	return 3.9 * Mpc * pow((1. + z) / 4., -5.0); // 4.5
}

double spectrum(const double& E_k) {
	double p = sqrt(E_k * E_k + 2. * mass_proton_c2 * E_k);
    double beta = p / (E_k + mass_proton_c2);
	double E_k_0 = reference_energy;
	double p_0 = sqrt(E_k_0 * E_k_0 + 2. * mass_proton_c2 * E_k_0);
	return 1. / beta * pow(p / p_0, -SN_slope);
}

/*double compute_source_normalization() {
	return 0;
}*/

double min_star_forming_halo(const double& z) {
	return 1e8 * mass_sun * pow(10. / (1.+z), 1.5);
}

double circular_velocity_at_rvir(const double& M, const double& z) {
	return 24.0 * km / s * pow(M / (1e8 * mass_sun), 1./3.) * pow((1. + z) / 10., 0.5);
}

double halo_mass_given_vvir(const double& v, const double& z) {
    double v_ = v / (24.0 * km / s);
    double z_ = pow((1. + z) / 10., -0.5);
	return 1e8 * mass_sun * pow(v_ * z_, 3.);
}

void print_timescales(string filename, const double& z) {
	cout << "Print timescales in " << filename << "\n";

	const double n_H = n_H_physical(z);
	const double n_He = f_He * n_H;
	const double ionization_fraction = 1e-4;
	const double n_e = ionization_fraction * n_H + 2. * ionization_fraction * n_He;
	const double n_HI = (1. - ionization_fraction) * n_H;

	cout << scientific << z << "\t" << n_H << "\t" << n_HI << "\t" << n_He << "\t" << ionization_fraction << "\t" << n_e << "\n";

	ofstream outfile;
	outfile.open(filename.c_str());
	outfile << "#E_k [GeV] - t_Hubble - t_ionization - t_Coulomb - t_adiabatic - t_fragmentation " << endl;
	outfile << scientific << setprecision(3);
	for (double E_k = MeV; E_k < 200. * GeV; E_k *= 1.2) {
		outfile << E_k / GeV << "\t" << fast::t_hubble(z) << "\t";
		outfile << spectrum(E_k) << "\t";
		outfile << E_k / dEdt_ionization(n_HI, E_k) << "\t";
		outfile << E_k / dEdt_coulomb_Galprop(n_e, E_k) << "\t" << E_k / dEdt_coulomb(n_e, E_k) << "\t";
		outfile << E_k / dEdt_adiabatic(z, E_k) << "\t" ;
		outfile << fragmentation_timescale(n_H) << "\n";
	}
	outfile.close();
}

double I_spectrum(double x, void * params) {
	double alpha = *(double *) params;
	double E_0 = *((double *)params + 1);
	double f = x * pow((E_0 * x * x + 2. * mass_proton_c2 * x) / (E_0 + 2. * mass_proton_c2), - .5 * alpha);
	return f;
}

double compute_spectrum_normalization(double E_0, double E_min, double E_max, double alpha) {

	gsl_integration_workspace * w
	= gsl_integration_workspace_alloc (10000);

	double result, error;
	double params[2] = {alpha, E_0};
	gsl_function F;
	F.function = &I_spectrum;
	F.params = params;

	gsl_integration_qag (&F, E_min / E_0, E_max / E_0, 0, 1e-5, 10000, 3, w, &result, &error);

	gsl_integration_workspace_free (w);

	return result;
}

double compute_initial_tau(const double& init_redshift){
	double result = 0;
	double dz = 0.1;
	for (double z = 1000; z > init_redshift; z -= dz) {
		double dt = -fast::dtdz(z) * dz;
		double n_e = 1e-4 * n_H_physical(z);
		result += SIGMAT * n_e * c_light * dt;
	}
	return result;
}
