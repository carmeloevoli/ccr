#include <fstream>
#include <iomanip>  // For std::scientific and std::setprecision
#include <iostream>

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
        for (double E_k = 1e-2 * cgs::MeV; E_k < 1e6 * cgs::MeV; E_k *= 1.1) {
            outfile << std::scientific << std::setprecision(6);
            outfile << E_k / cgs::MeV << "\t";
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
        std::ofstream outfile = open_output_file("diffusion_time_output.txt");

        const auto z = Z_VALUE;
        const auto B_IGM = 1e-16 * cgs::Gauss;
        const auto n_H = n_H_physical(z);
        const auto d = cgs::Mpc / (1. + z);
        const auto j = dtdz(z);

        // Loop over energy values and write data to file
        for (double E_k = 0.1 * cgs::MeV; E_k < 10. * cgs::GeV; E_k *= 1.1) {
            outfile << std::scientific << std::setprecision(6);
            outfile << E_k / cgs::MeV << "\t";             // 0
            outfile << hubble_time(z) / cgs::Gyr << "\t";  // 1
            const auto v = beta_lorentz(E_k) * cgs::c_light;
            outfile << v / (cgs::Mpc / cgs::Gyr) << "\t";                 // 2
            outfile << (d / v) / cgs::Gyr << "\t";                        // 3
            outfile << diffusion_time(B_IGM, d, E_k) / cgs::Gyr << "\t";  // 4
            outfile << inelastic_time(n_H, E_k) / cgs::Gyr << "\t";       // 5
            // outfile << -E_k / dEdz_H(z, E_k) * j / cgs::Gyr << "\t";      // 6
            outfile << E_k / dEdt_i(n_H, E_k) / cgs::Gyr << "\t";
            outfile << E_k / dEdt_C(n_H, E_k) / cgs::Gyr << "\t";
            outfile << E_k / dEdt_C(1e-3 * n_H, E_k) / cgs::Gyr << "\t";
            outfile << larmor_radius(B_IGM, E_k) / cgs::Mpc << "\t";
            outfile << "\n";
        }

        outfile.close();  // Close the file
    } catch (const std::ios_base::failure &e) {
        std::cerr << "File I/O error: " << e.what() << std::endl;
    }
}

int main() {
    print_energy_losses();
    // print_diffusion_time();
    return 0;
}
