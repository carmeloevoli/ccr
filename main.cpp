#include <iostream>

// #include "TEvolve.h"

#include "reionization.h"
//#include "utilities.h"

using namespace std;

double sigma_norm, R, theta_cmb, omhh, z_equality, y_d, sound_horizon, alpha_nu, f_nu, f_baryon, beta_c, d2fact, R_CUTOFF, DEL_CURR, SIG_CURR;

int main() {
    
//    cout << compute_spectrum_normalization(1. * GeV, 0.1 * GeV, 1e4 * GeV, 1e51 * erg, 2.2) << endl;

//    cout << compute_spectrum_normalization(1. * GeV, 0.1 * GeV, 1e5 * GeV, 1e51 * erg, 2.2) << endl;

//    cout << compute_spectrum_normalization(1. * GeV, 0.1 * GeV, 1e6 * GeV, 1e51 * erg, 2.2) << endl;

  //  cout << compute_spectrum_normalization(1. * GeV, 0.01 * GeV, 1e6 * GeV, 1e51 * erg, 2.2) << endl;

    //cout << compute_spectrum_normalization(1. * GeV, 0.001 * GeV, 1e6 * GeV, 1e51 * erg, 2.2) << endl;


    //print_timescales("output/timescales_at_z10.txt", 10);
    
    //print_timescales("output/timescales_at_z20.txt", 20);
    
    Reionization* R = new Reionization();
    
    //R->print_hmf(30, 1e6, 1e12);
    
    //cout << R->integrate_hmf(30, min_star_forming_halo(30), 1e12) << "\n";
    
    R->init_reionization();
    
    R->set_dz(0.01);
    
    R->set_f_sfr(0.04);
    
    R->set_f_lya(3e63 / mass_sun);
    
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
