#ifndef REIONIZATION_H_
#define REIONIZATION_H_

#include <iostream>
#include <sstream>

#include "ps.h"
#include "tridiag.h"
#include "utilities.h"

class Reionization {
public:
    Reionization(const string& init_filename_);
    ~Reionization();
    
    void init_reionization();
    void init_grids();
    void evolve(const bool& doCR);
    void evolve_CR(const double& dt);
    void evolve_IGM(const double& dt);
    void build_losses();

    double compute_ionization_rate_CR();
    double compute_heating_rate_CR();
    
    void print_status(bool doTitle);
    void dump_N(const double& z);
    void open_output_files();
    void close_output_files();
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
    
    inline void set_E_min(const double& E_min) {
        this->SN_E_min = E_min;
    }
    
    inline void set_E_max(const double& E_max) {
        this->SN_E_max = E_max;
    }
    
    inline void set_SN_slope(const double& slope) {
        this->SN_slope = slope;
    }
    
    inline void set_E_size(const int& E_size) {
        this->SN_E_size = E_size;
    }
    
private:
    double optical_depth_PLANCK;
    double z;
    double dz;
    double x_II;
    double n_H;
    double n_HI;
    double n_e;
    double T_k;
    double f_sfr;
    double f_esc;
    double ionization_rate;
    double ionization_rate_CR;
    double recombination_rate;
    double heating_rate;
    double heating_rate_CR;
    double star_formation_rate_comoving;
    double star_formation_rate_physical;
    double optical_depth;
    double sn_energy_rate;
    double cz;
    double normalization_integral;
    double SN_E_min;
    double SN_E_max;
    double SN_slope;
    
    size_t SN_E_size;
    
    vector<double> hmf_z;
    vector<double> hmf_integral;
    
    vector<double> E_k;
    vector<double> Q_sn;
    vector<double> b_losses;
    vector<double> N_cr;
    
    vector<double> knownTerm;
    vector<double> diagonal;
    vector<double> upperDiagonal;
    vector<double> lowerDiagonal;
    vector<double> N_next;
    
    string init_filename;
    
    ofstream fout_losses;
    ofstream fout_igm;
};

#endif