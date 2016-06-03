#include "utilities.h"

double dEdt_ionization(const double& n_H, const double& E_k) {
	const double gamma = 1.0 + E_k / mass_proton_c2;
	const double beta = sqrt(1.0 - 1.0 / pow2(gamma));
	const double p_squared = (pow2(gamma) - 1.);

	const double factor = 4. * pi_re2_me_c2_c / beta;

	const double A_H = Z_H * n_H * (log(2. * mass_electron_c2 * p_squared / ionization_potential_H) - pow2(beta));

	const double A_He = Z_He * f_He * n_H * (log(2. * mass_electron_c2 * p_squared / ionization_potential_He) - pow2(beta));

	return factor * (A_H + A_He);
}

double dEdt_coulomb(const double& n_e, const double& E_k) {
	const double gamma = 1.0 + E_k / mass_proton_c2;
	const double beta = sqrt(1.0 - 1.0 / pow2(gamma));
	const double p = sqrt(pow2(gamma) - 1.);

	const double factor = 4. * pi_re2_me_c2_c / beta;

	const double w_pl = 5.64e4 * sqrt(n_e) / s;

	const double Coulomb_log = log(2. * mass_electron_c2 * beta / h_bar_planck / w_pl * p);

	return factor * n_e * (Coulomb_log - pow2(beta) / 2.);
}

double dEdt_adiabatic(const double& z, const double& E_k) {
	const double dzdt = 1. / dtdz(z);

	return -dzdt / (1. + z) * E_k;
}

double fragmentation_timescale(const double& n_H) {
	return 6e7 * year / n_H; // arXiv:1208.4979
}

/* returns the proper mean baryonic density at z in cm^-3 */
double n_H_physical(const double & z) {
	return n_H_0 * pow(1. + z, 3.);
}

/* returns the hubble "constant" (in 1/sec) at z */
double hubble(const double& z){
	return H_0 * sqrt(Omega_m * pow3(1. + z) + Omega_r * pow4(1. + z) + Omega_l);
}

/* returns hubble time (in sec), t_h = 1/H */
double t_hubble(const double& z){
	return 1.0 / hubble(z);
}

/* function DTDZ returns the value of dt/dz at the redshift parameter z. */
double dtdz(const double& z){
  const double x = sqrt(Omega_l / Omega_m) * pow(1. + z, -1.5);
  const double dxdz = sqrt(Omega_l / Omega_m) * pow(1. + z, -2.5) * (-1.5);
  const double const1 = 2. * sqrt(1. + Omega_m / Omega_l) / (3.0 * H_0) ;

  const double numer = dxdz * (1. + x * pow(pow(x, 2) + 1, -0.5));
  const double denom = x + sqrt(pow(x, 2) + 1);

  return (const1 * numer / denom);
}

void print_timescales(string filename, const double& z) {
    cout << "Print timescales in " << filename << "\n";

    const double n_H = n_H_physical(z);
    const double n_He = f_He * n_H;
    const double ionization_fraction = 1e-4;
    const double n_e = ionization_fraction * n_H + 2. * ionization_fraction * n_He;

    cout << scientific << z << "\t" << n_H << "\t" << n_He << "\t" << ionization_fraction << "\t" << n_e << "\n";

    ofstream outfile;
    outfile.open(filename.c_str());
    outfile << "#E_k [GeV] - t_Hubble - t_ionization - t_Coulomb - t_adiabatic - t_fragmentation " << endl;
    outfile << scientific << setprecision(3);
    for (double E_k = MeV; E_k < 200. * GeV; E_k *= 1.2) {
    	double p =
        outfile << E_k / GeV << "\t" << t_hubble(z) << "\t" << E_k / dEdt_ionization(n_H, E_k) << "\t" << E_k / dEdt_coulomb(n_e, E_k) << "\t" <<  E_k / dEdt_adiabatic(z, E_k) << "\t" << fragmentation_timescale(n_H) << "\n";
    }
    outfile.close();
}
