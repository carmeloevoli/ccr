#include "reionization.h"

Reionization::Reionization(const string& init_filename_) {
    init_filename = init_filename_;
    f_sfr = 4e-2;
    f_esc = 0.1;
    SN_E_min = 100. * keV;
    SN_E_max = 1. * TeV;
    SN_slope = 2.2;
    SN_E_size = 7 * 32;
    open_output_files();
    
    //double z = 20;
    //double M = min_star_forming_halo(z);
    //double SFR = f_sfr * Omega_b / (Omega_m - Omega_b) * M / free_fall_timescale(z, M);
    //cout << z << "\t" << M / mass_sun << "\t" << SN_efficiency * SN_kinetic_energy * SN_fraction * SFR / (erg / s) << "\n";
}

Reionization::~Reionization() {
    close_output_files();
}

void Reionization::init_grids() {
    cout << "Init energ grids ...\n";
    double deltaLogE = exp(log(SN_E_max / SN_E_min) / (SN_E_size - 1));
    for (size_t iE = 0; iE < SN_E_size; ++iE) {
        double E_k_ = exp(log(SN_E_min) + iE * log(deltaLogE));
        E_k.push_back(E_k_);
        Q_sn.push_back(spectrum(E_k_, SN_slope));
    }
    b_losses.resize(SN_E_size);
    N_cr.resize(SN_E_size);
    knownTerm.resize(SN_E_size - 1);
    diagonal.resize(SN_E_size - 1);
    upperDiagonal.resize(SN_E_size - 2);
    lowerDiagonal.resize(SN_E_size - 2);
    N_next.resize(SN_E_size - 1);
}

void Reionization::init_reionization() {
    cout << "Init reionization parameters ...\n";
    x_II = 1e-4;
    T_k = 1e1 * K;
    z = initial_redshift;
    optical_depth = 0;
    heating_rate = 0;
    heating_rate_CR = 0;
    ionization_rate = 0;
    ionization_rate_CR = 0;
    recombination_rate = 0;
    optical_depth_PLANCK = 0.055 + 3. * 0.009;
    normalization_integral = compute_spectrum_normalization(reference_energy, SN_E_min, SN_E_max, SN_slope);
    //double initial_tau = compute_initial_tau(initial_redshift);
    cout << "... normalization integral is " << normalization_integral << "\n";
}

void Reionization::read_SFR(const string& filename) {
    cout << "Reading SFR from file " << filename << " ...\n";
    ifstream file_to_read(filename.c_str());
    const int num_of_header_lines = 1;
    
    if (file_to_read) {
        for (int i = 0; i < num_of_header_lines; ++i)
            file_to_read.ignore(512, '\n');
        
        while(!file_to_read.eof()) {
            double x, y1, y2;
            file_to_read >> x >> y1 >> y2;
            if (!file_to_read.eof()) {
                hmf_z.push_back(x);
                hmf_integral.push_back(y2);
            }
        }
        file_to_read.close();
    }
    else {
        cout << "File " << filename << " does not exist!\n";
        exit(1);
    }
    //cout << "HMF vector size = " << hmf_z.size() << endl;
}

double Reionization::hmf_integral_interpolate(const double& x) {
    if (x > 49) {
        cout << "z too large in interpolation!" << "\n";
        exit(1);
    }
    
    size_t i = 0;
    bool found = false;
    vector<double>::iterator it = hmf_z.begin();
    while (it != hmf_z.end() && !found) {
        it++;
        if (*it > x) {
            found = true;
            i = it - hmf_z.begin() - 1;
        }
    }
    
    double x_0 = hmf_z.at(i);
    double x_1 = hmf_z.at(i + 1);
    double y_0 = hmf_integral.at(i);
    double y_1 = hmf_integral.at(i + 1);
    double y = y_0 + (y_1 - y_0) * (x - x_0) / (x_1 - x_0);
    
    return y;
}


void Reionization::build_losses() {
    for (size_t iE = 0; iE < SN_E_size; ++iE) {
        double E_k_iE = E_k.at(iE);
        double dEdt = dEdt_adiabatic(z, E_k_iE);
        dEdt += dEdt_ionization(n_HI, E_k_iE);
        dEdt += dEdt_coulomb(n_e, E_k_iE);
        dEdt += dEdt_pp(n_H, E_k_iE);
        b_losses.at(iE) = -dEdt;
    }
}

void Reionization::evolve_IGM(const double& dt) {
    double A = (PopII_spectrum_slope - 1.) / (PopII_spectrum_slope + 0.5);
    double B = (PopII_spectrum_slope - 1.) / (pow2(PopII_spectrum_slope) - .25);
    
    n_H = n_H_physical(z);
    n_e =  x_II * (1. + f_He) * n_H;
    n_HI = (1. - x_II) * n_H;
    
    star_formation_rate_comoving = f_sfr * Omega_b / (Omega_m - Omega_b) * hmf_integral_interpolate(z); // M V^-1 T^-1
    star_formation_rate_physical = star_formation_rate_comoving * pow3(1. + z); // M V^-1 T^-1
    ionization_rate = A * (1. - x_II) * UV_photoionization_cs * UV_mean_free_path(z) * f_esc *  PopII_dNdM * star_formation_rate_physical; // T^-1
    recombination_rate = fast::alpha_A(T_k) * clumping_factor * (1. + f_He) * n_H * pow2(x_II); // L^3 T^-1 L^-3
    heating_rate = B * (13.6 * eV) * UV_photoionization_cs * UV_mean_free_path(z) * f_esc *  PopII_dNdM * star_formation_rate_physical * (1. - x_II); // E / T
    
    double dx_II = dt * (ionization_rate - recombination_rate); // cm^3 s^-1 cm^-3;
    double dT_k_dz = 2. * T_k / (1. + z) - T_k / (1. + x_II) * (dx_II / dz) + 2. / 3. / k_boltzmann / (1. + x_II) * fast::dtdz(z) * heating_rate ;
    
    x_II -= dx_II;
    T_k -= dT_k_dz * dz;
    
    optical_depth += SIGMAT * n_e * c_light * -dt;
}

double Reionization::compute_ionization_rate_CR() {
    double sum = 0;
    double delta_lnE = log(E_k.at(1) / E_k.at(0));
    for (size_t iE = 0; iE < E_k.size(); iE++) {
        double E_k_iE = E_k.at(iE);
        sum += dEdt_ionization(n_HI, E_k_iE) * N_cr.at(iE) * E_k_iE;
    }
    
    return delta_lnE * sum / W_H;
}

double Reionization::compute_heating_rate_CR() {
    double sum_C = 0;
    double sum_I = 0;
    double delta_lnE = log(E_k.at(1) / E_k.at(0));
    for (size_t iE = 0; iE < E_k.size(); iE++) {
        double E_k_iE = E_k.at(iE);
        sum_C += dEdt_coulomb(n_e, E_k_iE) * N_cr.at(iE) * E_k_iE;
        sum_I += dEdt_ionization(n_HI, E_k_iE) * N_cr.at(iE) * E_k_iE;
    }
    
    return delta_lnE * (sum_C +  sum_I);
}

void Reionization::evolve_CR(const double& dt) {
    sn_energy_rate = SN_efficiency * SN_kinetic_energy * SN_fraction * star_formation_rate_physical; // E V^-1 T^-1
    cz = sn_energy_rate / pow2(reference_energy) / normalization_integral;
    n_H = n_H_physical(z);
    
    build_losses();
    
    for (size_t iE = 0; iE < SN_E_size - 1; ++iE) {
        double alpha2 = - fabs(dt) * b_losses.at(iE) / (E_k.at(iE+1) - E_k.at(iE));
        double alpha3 = - fabs(dt) * b_losses.at(iE+1) / (E_k.at(iE+1) - E_k.at(iE));
        diagonal.at(iE) = .5 * (1. + alpha2);
        if (iE != SN_E_size - 2)
            upperDiagonal.at(iE) = .5 * (1. - alpha3);
        if (iE != 0)
            lowerDiagonal.at(iE - 1) = 0;
        knownTerm.at(iE) = .5 * (1. - alpha2) * N_cr.at(iE);
        knownTerm.at(iE) += .5 * (1. + alpha3) * N_cr.at(iE + 1);
        knownTerm.at(iE) += .5 * fabs(dt) * cz * (Q_sn.at(iE + 1) + Q_sn.at(iE)) / n_H;
    }
    
    gsl_linalg_solve_tridiag(diagonal, upperDiagonal, lowerDiagonal, knownTerm, N_next);
    
    for (size_t iE = 0; iE < SN_E_size - 1; ++iE) {
        N_cr.at(iE) = max(N_next.at(iE), 0.);
    }
    
    ionization_rate_CR = compute_ionization_rate_CR();
    heating_rate_CR = compute_heating_rate_CR();
    
    //cout << dt << "\t" << *min_element(N_cr.begin(),N_cr.end()) << "\t" << *max_element(N_cr.begin(),N_cr.end()) << "\n";
    
    return;
}

void Reionization::plot_source_function(const double& z) {
    cout << "Plotting source function ...\n";
    
    star_formation_rate_comoving = f_sfr * Omega_b / (Omega_m - Omega_b) * hmf_integral_interpolate(z); // M V^-1 T^-1
    star_formation_rate_physical = star_formation_rate_comoving * pow3(1. + z); // M V^-1 T^-1
    sn_energy_rate = SN_efficiency * SN_kinetic_energy * SN_fraction * star_formation_rate_physical; // E V^-1 T^-1
    cz = sn_energy_rate / pow2(reference_energy) / normalization_integral; // E^-1 V^-1 T^-1
    
    open_spectrum_file(z);
    fout_spectrum << "#E [erg] - Q [] - Cz [erg-1 cm-3 s-1] \n";
    for (size_t iE = 0; iE < SN_E_size - 1; ++iE) {
        fout_spectrum << setprecision(5) << scientific << E_k.at(iE) << "\t" << Q_sn.at(iE) << "\t" << cz << "\n";
    }
    close_spectrum_file();
}

void Reionization::evolve(const bool& doCR) {
    
    cout << "Begin solver ...\n";
    
    size_t counter = 0;
    
    print_status(true);
    
    while (z > 0) {
        
        counter++;

        double dt = fast::dtdz(z) * dz;
        
        evolve_IGM(dt);
        
        if (doCR)
            evolve_CR(dt);
        
        z -= dz;
        
        if (counter % (int)(0.1/dz) == 0) {
            print_status(false);
            if (doCR)
                dump_N(z);
        }
    }
    cout << "... end in " << 1 << " s.\n";
}
