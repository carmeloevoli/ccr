#include <fstream>
#include <iomanip>
#include <iostream>

#include "ccrh.h"
#include "constants.h"
#include "utilities.h"

// Constants for z
const double Z_VALUE = 12.0;

// Function to open file and return output stream
std::ofstream open_output_file(const std::string &filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        throw std::ios_base::failure("Failed to open output file.");
    }
    return outfile;
}

void print_energy_losses() {
    try {
        std::ofstream outfile = open_output_file("losses_timescales_output.txt");
        const auto n_H = n_H_physical(Z_VALUE);
        for (double E_k = 1e-3 * cgs::MeV; E_k < 1e6 * cgs::MeV; E_k *= 1.1) {
            outfile << std::scientific << std::setprecision(6);
            outfile << E_k / cgs::MeV << "\t";
            outfile << beta_lorentz(E_k) * cgs::c_light / (cgs::Mpc / cgs::year) << "\t";
            outfile << E_k / dEdt_C(n_H, E_k) / cgs::year << "\t";
            outfile << E_k / dEdt_C(0.0014 * n_H, E_k) / cgs::year << "\t";
            outfile << E_k / dEdt_i(n_H, E_k) / cgs::year << "\t";
            outfile << E_k / dEdt_a(Z_VALUE, E_k) / cgs::year << "\t";
            outfile << "\n";
        }
    } catch (const std::ios_base::failure &e) {
        std::cerr << "File I/O error: " << e.what() << std::endl;
    }
}

// Main function to print diffusion time
void print_diffusion_time() {
    try {
        std::ofstream outfile = open_output_file("diffusion_timescales_output.txt");
        const auto z = Z_VALUE;
        const auto B_IGM = 1e-16 * cgs::Gauss;
        const auto d = cgs::Mpc / (1. + z);
        for (double E_k = 1e-2 * cgs::MeV; E_k < 1e6 * cgs::MeV; E_k *= 1.1) {
            outfile << std::scientific << std::setprecision(6);
            outfile << E_k / cgs::MeV << "\t";             // 0
            outfile << hubble_time(z) / cgs::Gyr << "\t";  // 1
            const auto v = beta_lorentz(E_k) * cgs::c_light;
            outfile << (d / v) / cgs::Gyr << "\t";                        // 3
            outfile << diffusion_time(B_IGM, d, E_k) / cgs::Gyr << "\t";  // 4
            outfile << "\n";
        }
        outfile.close();  // Close the file
    } catch (const std::ios_base::failure &e) {
        std::cerr << "File I/O error: " << e.what() << std::endl;
    }
}

void print_spectrum() {
    try {
        std::cout << SN_ESpectrum_integral(2.4) << "\n";
        const auto E_0 = cgs::mass_proton_c2;
        const auto L = 1e-33 * cgs::erg / cgs::cm3 / cgs::s;

        std::ofstream outfile = open_output_file("energy_spectrum_output.txt");
        for (double E_k = 1e-3 * cgs::MeV; E_k < 1e9 * cgs::MeV; E_k *= 1.1) {
            outfile << std::scientific << std::setprecision(6);
            outfile << E_k / cgs::MeV << "\t";
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
        std::cerr << "File I/O error: " << e.what() << std::endl;
    }
}

void ccrh() {
    CCRH ccrh;
    // ccrh.init_vectors();
    // ccrh.evolve();
}

int main() {
    print_energy_losses();
    // print_diffusion_time();
    print_spectrum();
    // ccrh();
    return 0;
}
