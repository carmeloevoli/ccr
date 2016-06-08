#include <iostream>

// #include "TEvolve.h"

#include "reionization.h"
//#include "utilities.h"

using namespace std;

double sigma_norm, R, theta_cmb, omhh, z_equality, y_d, sound_horizon, alpha_nu, f_nu, f_baryon, beta_c, d2fact, R_CUTOFF, DEL_CURR, SIG_CURR;

int main() {

    
	//print_timescales("output/timescales_at_z10.txt", 10);

	//print_timescales("output/timescales_at_z20.txt", 20);

    Reionization* R = new Reionization();
    
    //R->print_hmf(30, 1e6, 1e12);
    
    //cout << R->integrate_hmf(30, min_star_forming_halo(30), 1e12) << "\n";
    
    R->init_reionization();
    
    R->set_dz(0.5);
    
    R->evolve();
    
    delete R;
    
/*	TEvolve* E = new TEvolve(400,"var");

	E->init_energy_grid();
	E->init_equation_vectors();
	E->init_time_scales_vectors();
	E->print_time_scales();

	E->evolve();

	delete E;
*/
	return 0;
}
