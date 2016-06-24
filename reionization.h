#ifndef REIONIZATION_H_
#define REIONIZATION_H_

#include "ps.h"
#include "utilities.h"

class Reionization {
public:
    Reionization();
    ~Reionization();
    
    void init_reionization();
    void evolve();
    void evolve_CR(const double& dt);
    void evolve_IGM(const double& dt);
    void print_status(bool doTitle);
    void read_SFR(const string& filename);
    double hmf_integral_interpolate(const double& z);

    inline void set_dz(const double& dz) {
        this->dz = dz;
    }
    
    inline void set_xII(const double& x_II) {
        this->x_II = x_II;
    }
    
    inline void set_Tk(const double& T_k) {
        this->T_k = T_k;
    }
    
    inline void set_f_sfr(const double& f_sfr) {
        this->f_sfr = f_sfr;
    }
    
    inline void set_f_esc(const double& f_esc) {
        this->f_esc = f_esc;
    }
    
private:
    double optical_depth_PLANCK;
    double z;
    double dz;
    double x_II;
    double n_H;
    double n_e;
    double n_n;
    double T_k;
    double f_sfr;
    double f_esc;
    double ionization_rate;
    double recombination_rate;
    double heating_rate;
    double star_formation_rate_comoving;
    double star_formation_rate_physical;
    double optical_depth;
    double sn_energy_rate;
    double cz;
    double normalization_integral;

    vector<double> hmf_z;
    vector<double> hmf_integral;
};

#endif