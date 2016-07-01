#include <iostream>

// #include "TEvolve.h"

#include "reionization.h"
#include "SFR.h"
#include "utilities.h"

using namespace std;

double sigma_norm, R, theta_cmb, omhh, z_equality, y_d, sound_horizon, alpha_nu, f_nu, f_baryon, beta_c, d2fact, R_CUTOFF, DEL_CURR, SIG_CURR;

#define v(A) (r*r+1.)

int main() {
    
    const double dt = 0.1;

    double r = 0;
    
    double t = 0;

    while (r < 10) {
    
    double dr = v(r) * dt;
 
        r += dr;
        
        t += dt;
        
        cout << t << "\t" << r << "\n";
    }
    
    
    
//    cout << compute_spectrum_normalization(1. * GeV, 0.1 * GeV, 1e4 * GeV, 1e51 * erg, 2.2) << endl;

//    cout << compute_spectrum_normalization(1. * GeV, 0.1 * GeV, 1e5 * GeV, 1e51 * erg, 2.2) << endl;

//    cout << compute_spectrum_normalization(1. * GeV, 0.1 * GeV, 1e6 * GeV, 1e51 * erg, 2.2) << endl;

  //  cout << compute_spectrum_normalization(1. * GeV, 0.01 * GeV, 1e6 * GeV, 1e51 * erg, 2.2) << endl;

    //cout << compute_spectrum_normalization(1. * GeV, 0.001 * GeV, 1e6 * GeV, 1e51 * erg, 2.2) << endl;

    //print_timescales("output/timescales_at_z10.txt", 10);
    
    //print_timescales("output/timescales_at_z20.txt", 20);
    
    fast::init_ps();
    
    //for (double z = 0; z < 30; z += 1) {
    //    double l = UV_mean_free_path(z);
        //cout << z << "\t" << l / Mpc << "\t" << c_light / sqrt(M_PI) / l / fast::hubble(z) / (1+z) << "\t" << 39. * pow((1. + z) / 4., -4.5) << "\n";
    //}

    //SFR* S = new SFR("SFR.txt");
    
    //S->print_hmf(30, 1e6, 1e12);
    
    //S->evolve();
    
    //delete S;
    
    Reionization* R = new Reionization("test_no_CR");
    
    R->read_SFR("SFR.txt");
    
    R->init_reionization();
    
    R->set_dz(1e-6);
    
    R->set_f_sfr(0.04);
    
    R->set_f_esc(2e-4);
        
    //R->evolve(false);
    
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
