#include <fstream>
#include <iomanip>
#include <iostream>

#include "io.h"
#include "utilities.h"

// Constants for z
const double Z_VALUE = 12.0;

// Function declarations
void print_energy_losses();
void print_diffusion_time();
void print_spectrum();

int main() {
    std::cout << "Starting calculations...\n";

    // Call the functions in sequence
    print_energy_losses();
    std::cout << "Energy losses calculated and saved to 'losses_timescales_output.txt'.\n";

    print_diffusion_time();
    std::cout << "Diffusion time calculated and saved to 'diffusion_timescales_output.txt'.\n";

    print_spectrum();
    std::cout << "Energy spectrum calculated and saved to 'energy_spectrum_output.txt'.\n";

    std::cout << "All calculations completed.\n";
    return 0;
}

// Energy loss calculations
void print_energy_losses() {
    try {
        std::ofstream outfile = open_output_file("losses_timescales_output.txt");
        const auto n_H = n_H_physical(Z_VALUE);
        for (double E_k = 1e-3 * cgs::MeV; E_k < 1e6 * cgs::MeV; E_k *= 1.1) {
            outfile << std::scientific << std::setprecision(6);
            outfile << E_k / cgs::MeV << "\t";  // Energy in MeV
            outfile << beta_lorentz(E_k) * cgs::c_light / (cgs::Mpc / cgs::year)
                    << "\t";                                        // Lorentz velocity
            outfile << E_k / dEdt_C(n_H, E_k) / cgs::year << "\t";  // Coulomb losses (n_H)
            outfile << E_k / dEdt_C(0.0014 * n_H, E_k) / cgs::year
                    << "\t";                                        // Coulomb losses (0.0014 * n_H)
            outfile << E_k / dEdt_i(n_H, E_k) / cgs::year << "\t";  // Ionization losses
            outfile << E_k / dEdt_a(Z_VALUE, E_k) / cgs::year << "\t";  // Adiabatic losses
            outfile << "\n";
        }
    } catch (const std::ios_base::failure &e) {
        std::cerr << "File I/O error in energy losses: " << e.what() << std::endl;
    }
}

// Diffusion time calculations
void print_diffusion_time() {
    try {
        std::ofstream outfile = open_output_file("diffusion_timescales_output.txt");
        const auto z = Z_VALUE;
        const auto B_IGM = 1e-16 * cgs::Gauss;
        const auto d = cgs::Mpc / (1. + z);
        for (double E_k = 1e-2 * cgs::MeV; E_k < 1e6 * cgs::MeV; E_k *= 1.1) {
            outfile << std::scientific << std::setprecision(6);
            outfile << E_k / cgs::MeV << "\t";             // Energy in MeV
            outfile << hubble_time(z) / cgs::Gyr << "\t";  // Hubble time
            const auto v = beta_lorentz(E_k) * cgs::c_light;
            outfile << (d / v) / cgs::Gyr << "\t";                        // Particle crossing time
            outfile << diffusion_time(B_IGM, d, E_k) / cgs::Gyr << "\t";  // Diffusion time
            outfile << "\n";
        }
        outfile.close();
    } catch (const std::ios_base::failure &e) {
        std::cerr << "File I/O error in diffusion time: " << e.what() << std::endl;
    }
}

// Spectrum calculations
void print_spectrum() {
    try {
        std::ofstream outfile = open_output_file("energy_spectrum_output.txt");
        const auto E_0 = cgs::mass_proton_c2;
        const auto L = 1e-33 * cgs::erg / cgs::cm3 / cgs::s;  // Energy source term

        for (double E_k = 1e-3 * cgs::MeV; E_k < 1e9 * cgs::MeV; E_k *= 1.1) {
            outfile << std::scientific << std::setprecision(6);
            outfile << E_k / cgs::MeV << "\t";  // Energy in MeV

            // Normalized energy spectra for different slopes
            auto Q_norm = L / pow2(E_0) / SN_ESpectrum_integral(2.0);
            outfile << pow2(E_k) * Q_norm * SN_Spectrum(E_k, 2.0) / cgs::erg << "\t";

            Q_norm = L / pow2(E_0) / SN_ESpectrum_integral(2.4);
            outfile << pow2(E_k) * Q_norm * SN_Spectrum(E_k, 2.4) / cgs::erg << "\t";

            Q_norm = L / pow2(E_0) / SN_ESpectrum_integral(2.7);
            outfile << pow2(E_k) * Q_norm * SN_Spectrum(E_k, 2.7) / cgs::erg << "\t";
            outfile << "\n";
        }
        outfile.close();
    } catch (const std::ios_base::failure &e) {
        std::cerr << "File I/O error in spectrum: " << e.what() << std::endl;
    }
}
