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
    void open_spectrum_file(const double& z);
    void close_output_files();
    void close_spectrum_file();
    void read_SFR(const string& filename);
    double hmf_integral_interpolate(const double& z);
    void plot_source_function(const double& z);

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
    ofstream fout_spectrum;
};

#endif
