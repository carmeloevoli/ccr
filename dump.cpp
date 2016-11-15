#include "reionization.h"

void Reionization::open_output_files() {
    string losses_filename = "output/";
    losses_filename += init_filename;
    losses_filename += "_losses.txt";
    fout_losses.open(losses_filename.c_str());
    
    string igm_filename = "output/";
    igm_filename += init_filename;
    igm_filename += "_igm.txt";
    fout_igm.open(igm_filename.c_str());
}

void Reionization::open_spectrum_file(const double& z) {
    string filename = "output/";
    filename += init_filename;
    filename += "_spectrum.txt";
    fout_spectrum.open(filename.c_str());
}

void Reionization::close_output_files() {
    fout_losses.close();
    fout_igm.close();
}

void Reionization::close_spectrum_file() {
    fout_spectrum.close();
}

void Reionization::dump_N(const double& z) {
    cout << "Dump spectrum at z = " << z << "\n";
    stringstream sstream;
    string filename;
    sstream << "output/" << init_filename << "_spectra_" << E_size << "_at_" << z << "_fesc_" << f_esc << "_fsfr_" << f_sfr << ".txt";
    filename = sstream.str();
    ofstream outfile;
    outfile.open(filename.c_str());
    outfile << "#E - N " << endl;
    outfile << scientific << setprecision(2);
    for (size_t i = 0; i < E_k.size(); ++i) {
        outfile << E_k.at(i) << "\t" << N_cr.at(i) << "\t" << cz * Q_sn.at(i) << "\t" << n_H << "\n";
    }
    outfile.close();
}

void Reionization::print_status(bool doTitle) {
    if (doTitle) {
        fout_losses << "#z - x_II - t_I[E_0] - t_C[E_0] - t_a[E_0] - t_I[E_1] - t_C[E_1] - t_a[E_1] - t_H" << "\n";
        fout_igm << "#z - x_II - T_k[K] - optical_depth - SFR [M_sun Mpc^-3 yr^-1] - Ion Rate [Myr^-1] - Rec Rate [Myr^-1]" << "\n";
    }
    else if (z > 0) {
        fout_losses << scientific << setprecision(3);
        fout_losses << z << "\t";
        fout_losses << x_II << "\t";
        double E_k = 1. * MeV;
        fout_losses << (E_k / dEdt_ionization(n_HI, E_k)) / Gyr << "\t";
        fout_losses << (E_k / dEdt_coulomb(n_e, E_k)) / Gyr << "\t";
        fout_losses << (E_k / dEdt_pp(n_H, E_k)) / Gyr << "\t";
        fout_losses << (E_k / dEdt_adiabatic(z, E_k)) / Gyr << "\t";
        E_k = 10. * MeV;
        fout_losses << (E_k / dEdt_ionization(n_HI, E_k)) / Gyr << "\t";
        fout_losses << (E_k / dEdt_coulomb(n_e, E_k)) / Gyr << "\t";
        fout_losses << (E_k / dEdt_pp(n_H, E_k)) / Gyr << "\t";
        fout_losses << (E_k / dEdt_adiabatic(z, E_k)) / Gyr << "\t";
        fout_losses << hubble_time(z) / Gyr << "\t";
        //cout << fast::t_hubble(z) / Gyr << "\t";
        fout_losses << "\n";
        
        fout_igm << scientific << setprecision(3);
        fout_igm << z << "\t";
        fout_igm << x_II << "\t";
        fout_igm << T_k << "\t";
        fout_igm << optical_depth << "\t";
        //cout << optical_depth_PLANCK << "\t";
        //fout_igm << min_sfr_halo / mass_sun << "\t";
        //fout_igm << hmf_integral / (mass_sun / pow3(Mpc) / year) << "\t";
        fout_igm << star_formation_rate_comoving / (mass_sun / pow3(Mpc) / year) << "\t";
        fout_igm << sn_energy_rate / (erg / pow3(cm) / s) << "\t";
        //cout << cz / (1. / erg / pow3(cm) / s) << "\t";
        fout_igm << ionization_rate / (1. / Myr) << "\t";
        fout_igm << ionization_rate_CR / (1. / Myr) << "\t";
        fout_igm << recombination_rate / (1. / Myr) << "\t";
        fout_igm << heating_rate / (1. / Myr) << "\t";
        fout_igm << heating_rate_CR / (1. / Myr) << "\t";
        fout_igm << "\n";
    }
}

