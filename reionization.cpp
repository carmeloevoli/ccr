#include "reionization.h"

Reionization::Reionization() {
    init_reionization();
    clumping_factor = 2.;
    f_sfr = 4e-2;
    f_esc = 0.1;
    optical_depth_PLANCK = 0.055 + 3. * 0.009;
}

Reionization::~Reionization() {}

void Reionization::init_reionization() {
    x_II = 1e-4;
    T_k = 1e1 * K;
    z = initial_redshift;
    optical_depth = 0;
    heating_rate = 0;
    ionization_rate = 0;
    recombination_rate = 0;
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

void Reionization::evolve() {
    double dt = 0;
    double dx_II = 0;
    double dT_k_dz = 0;
    double n_e = 0;
    
    size_t counter = 0;
    
    print_status(true);
    
    const double UV_photoionization_cs = 6.3e-18 * pow2(cm);
    const double PopII_spectrum_slope = 5.;
    const double PopII_dNdM = 8e60 / mass_sun;
    const double A = (PopII_spectrum_slope - 1.) / (PopII_spectrum_slope + 0.5);
    const double B = (PopII_spectrum_slope - 1.) / (pow2(PopII_spectrum_slope) - .25);
    
    double normalization_integral = compute_spectrum_normalization(reference_energy, SN_E_min, SN_E_max, SN_slope);
    
    //double initial_tau = compute_initial_tau(initial_redshift);
    
    while (z > 0) {
        
        counter++;
        
        dt = fast::dtdz(z) * dz;
        
        star_formation_rate_comoving = f_sfr * Omega_b / (Omega_m - Omega_b) * hmf_integral_interpolate(z); // M V^-1 T^-1
        star_formation_rate_physical = star_formation_rate_comoving * pow3(1. + z); // M V^-1 T^-1
        //ionization_rate = c_light * UV_photoionization_cs * f_UV * (1. - x_II) * star_formation_rate_physical * fabs(dt); // T^-1
        ionization_rate = A * (1. - x_II) * UV_photoionization_cs * UV_mean_free_path(z) * f_esc *  PopII_dNdM * star_formation_rate_physical; // T^-1
        recombination_rate = fast::alpha_A(T_k) * clumping_factor * (1. + f_He) * n_H_physical(z) * pow2(x_II); // L^3 T^-1 L^-3
        heating_rate = B * (13.6 * eV) * UV_photoionization_cs * UV_mean_free_path(z) * f_esc *  PopII_dNdM * star_formation_rate_physical * (1. - x_II); // E / T
        
        sn_energy_rate = SN_efficiency * SN_kinetic_energy * SN_fraction * star_formation_rate_physical; // E V^-1 T^-1
        cz = sn_energy_rate / pow2(reference_energy) / normalization_integral;
        
        dx_II = dt * (ionization_rate - recombination_rate); // cm^3 s^-1 cm^-3;
        
        dT_k_dz = 2. * T_k / (1. + z) - T_k / (1. + x_II) * (dx_II / dz) + 2. / 3. / k_boltzmann / (1. + x_II) * fast::dtdz(z) * heating_rate ;
        
        x_II -= dx_II;
        T_k -= dT_k_dz * dz;
        z -= dz;
        
        n_e = (1. + f_He) * x_II * n_H_physical(z);
        
        optical_depth += SIGMAT * n_e * c_light * -dt;
        
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
        cout << optical_depth_PLANCK << "\t";
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
