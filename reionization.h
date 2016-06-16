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
    void print_status(bool doTitle);
    
    void print_hmf(const double& z, const double& M_min, const double& M_max);
    double integrate_hmf(double z, const double& M_min, const double& M_max);
    double Galaxy_ionization_rate(const double& z);
    double Galaxy_heating_rate(const double& z);
    
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
    
    inline void set_f_lya(const double& f_lya) {
        this->f_lya = f_lya;
    }
    
private:
    double z;
    double dz;
    double x_II;
    double T_k;
    double f_sfr;
    double f_lya;
    double clumping_factor;
    double min_sfr_halo;
    double hmf_integral;
    double ionization_rate;
    double recombination_rate;
    double star_formation_rate_comoving;
    double star_formation_rate_physical;
    double optical_depth;
    double sn_energy_rate;
    double cz;
};

#endif