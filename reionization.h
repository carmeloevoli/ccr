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
    void print_status();
    
    void print_hmf(const double& z, const double& M_min, const double& M_max);
    double integrate_hmf(double z, const double& M_min, const double& M_max);
    double Galaxy_ionization_rate(const double& z);
    double Galaxy_heating_rate(const double& z);

    inline void set_dz(const double& dz) {
        this->dz = dz;
    }
    
    inline void set_xe(const double& x_e) {
        this->x_e = x_e;
    }
    
    inline void set_Tk(const double& T_k) {
        this->T_k = T_k;
    }
    
    
private:
    double z;
    double dz;
    double x_e;
    double T_k;
    double f_sfr;
    double f_alpha;
    double clumping_factor;
    double min_sfr_halo;
    double mass_in_halos;
    double ionization_rate;
    double recombination_rate;
    double star_formation_rate;
};

#endif