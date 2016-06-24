#include "reionization.h"

Reionization::Reionization() {
    init_reionization();
    init_grids();
    f_sfr = 4e-2;
    f_esc = 0.1;
}

Reionization::~Reionization() {}

void Reionization::init_grids() {
    double deltaLogE = exp(log(SN_E_max / SN_E_min) / (E_size - 1));
    for (size_t iE = 0; iE < E_size; ++iE) {
        E_k.push_back(exp(log(SN_E_min) + iE * log(deltaLogE)));
    }
    Q_SN.resize(E_size, 0.);
    b_losses.resize(E_size, 0.);
}


void Reionization::init_reionization() {
    x_II = 1e-4;
    T_k = 1e1 * K;
    z = initial_redshift;
    optical_depth = 0;
    heating_rate = 0;
    ionization_rate = 0;
    recombination_rate = 0;
    normalization_integral = compute_spectrum_normalization(reference_energy, SN_E_min, SN_E_max, SN_slope);
    optical_depth_PLANCK = 0.055 + 3. * 0.009;
    //double initial_tau = compute_initial_tau(initial_redshift);
}

void Reionization::read_SFR(const string& filename) {
    ifstream file_to_read(filename.c_str());
    const int num_of_header_lines = 1;
    
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
    
    return;
}

void Reionization::evolve_CR(const double& dt) {
    sn_energy_rate = SN_efficiency * SN_kinetic_energy * SN_fraction * star_formation_rate_physical; // E V^-1 T^-1
    cz = sn_energy_rate / pow2(reference_energy) / normalization_integral;
    
    
    
    return;
}

void Reionization::evolve() {
    
    size_t counter = 0;
    
    print_status(true);
    
    while (z > 0) {
        
        counter++;
        
        double dt = fast::dtdz(z) * dz;
        
        evolve_IGM(dt);
        
        evolve_CR(dt);
        
        z -= dz;
        
        if (counter % 100000 == 0)
            print_status(false);
    }
}

void Reionization::print_status(bool doTitle) {
    if (doTitle) {
        cout << "#z - x_e - T_k [K] - optical_depth - min_sfr_halo [M_sun] - halo_integral [M_sun Mpc^-3 yr^-1] - SFR [M_sun Mpc^-3 yr^-1] - Ion Rate [Myr^-1] - Rec Rate [Myr^-1]" << "\n";
    }
    else if (z > 0) {
        cout << scientific << setprecision(3);
        cout << z << "\t";
        cout << x_II << "\t";
        cout << T_k << "\t";
        cout << optical_depth << "\t";
        //cout << optical_depth_PLANCK << "\t";
        double E_k = 1. * MeV;
        cout << (E_k / dEdt_ionization(n_HI, E_k)) / Gyr << "\t";
        cout << (E_k / dEdt_coulomb(n_e, E_k)) / Gyr << "\t";
        E_k = 10. * MeV;
        cout << (E_k / dEdt_ionization(n_HI, E_k)) / Gyr << "\t";
        cout << (E_k / dEdt_coulomb(n_e, E_k)) / Gyr << "\t";
        cout << hubble_time(z) / Gyr << "\t";
        //cout << fast::t_hubble(z) / Gyr << "\t";
        //cout << min_sfr_halo / mass_sun << "\t";
        //cout << hmf_integral / (mass_sun / pow3(Mpc) / year) << "\t";
        cout << star_formation_rate_comoving / (mass_sun / pow3(Mpc) / year) << "\t";
        cout << sn_energy_rate / (erg / pow3(cm) / s) << "\t";
        cout << cz / (1. / erg / pow3(cm) / s) << "\t";
        cout << ionization_rate / (1. / Myr) << "\t";
        cout << recombination_rate / (1. / Myr) << "\t";
        cout << "\n";
    }
}
