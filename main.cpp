#include <iostream>

#include "reionization.h"
#include "SFR.h"
#include "utilities.h"

using namespace std;

double sigma_norm, R, theta_cmb, omhh, z_equality, y_d, sound_horizon, alpha_nu, f_nu, f_baryon, beta_c, d2fact, R_CUTOFF, DEL_CURR, SIG_CURR;

#define v(A) (r*r+1.)

void simple_main();

int main() {

	bool doSFR = false;
	bool doReionization = true;

	fast::init_ps();
    
    simple_main();
    
	if (doSFR) {
		SFR* S = new SFR("SFR_new_60.txt");

        S->set_vvir_cut(60. * km / s);
        
		//S->print_hmf(10, 1e6, 1e15);

		S->evolve();

		delete S;
	}

	if (doReionization) {
		Reionization* R = new Reionization("test_new_2.5");

		R->read_SFR("SFR_new_100.txt");

        R->set_dz(1e-7);

		R->set_f_sfr(0.02);

		R->set_f_esc(1e-2);

        R->set_SN_slope(2.5);

        R->init_grids();
        
        R->init_reionization();
        
        //R->plot_source_function(10);

        R->evolve(true);

		delete R;
	}

	return 0;
}

void simple_main() {

	//cout << Omega_b / Omega_m << "\n";

	//double z = 20.;

	//double M = min_star_forming_halo(z);

	//cout << z << "\t" << M / mass_sun << "\t" << free_fall_timescale(z,M) << "\n";
	
	double alpha = 2.5;

	double E_0 = 1.0 * GeV;

	cout << alpha << "\t";

	cout << compute_spectrum_normalization(E_0, 10. * keV, .01 * GeV, alpha) << "\t";

    cout << compute_spectrum_normalization(E_0, 10. * keV, 1e6 * GeV, alpha) << "\t";

	cout << compute_spectrum_normalization(E_0, 10. * keV, .01 * GeV, alpha) / compute_spectrum_normalization(E_0, 10. * keV, 1e6 * GeV, alpha) * 100.;

	cout << "\n";

	//cout << compute_spectrum_normalization(1. * GeV, 0.1 * GeV, 1e4 * GeV, 1e51 * erg, 2.2) << endl;

	//cout << compute_spectrum_normalization(1. * GeV, 0.1 * GeV, 1e5 * GeV, 1e51 * erg, 2.2) << endl;

	//cout << compute_spectrum_normalization(1. * GeV, 0.1 * GeV, 1e6 * GeV, 1e51 * erg, 2.2) << endl;

	//cout << compute_spectrum_normalization(1. * GeV, 0.01 * GeV, 1e6 * GeV, 1e51 * erg, 2.2) << endl;

	//cout << compute_spectrum_normalization(1. * GeV, 0.001 * GeV, 1e6 * GeV, 1e51 * erg, 2.2) << endl;

	//print_timescales("output/timescales_at_z10.txt", 10);

	//print_timescales("output/timescales_at_z20.txt", 20);

	//for (double z = 0; z < 30; z += 1) {
	//double l = UV_mean_free_path(z);
	//cout << z << "\t" << l / Mpc << "\t" << c_light / sqrt(M_PI) / l / fast::hubble(z) / (1+z) << "\t" << 39. * pow((1. + z) / 4., -4.5) << "\n";
	//}

	//cout << 1. / H_0 / sqrt(Omega_m) / n_H_0 << "\n";

    exit(1);

	return;
}
