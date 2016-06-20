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
    
    inline void set_f_UV(const double& f_UV) {
        this->f_UV = f_UV;
    }
    
private:
    double optical_depth_PLANCK;
    double z;
    double dz;
    double x_II;
    double T_k;
    double f_sfr;
    double f_UV;
    double clumping_factor;
    double ionization_rate;
    double recombination_rate;
    double heating_rate;
    double star_formation_rate_comoving;
    double star_formation_rate_physical;
    double optical_depth;
    double sn_energy_rate;
    double cz;

    vector<double> hmf_z;
    vector<double> hmf_integral;
    
};

#endif