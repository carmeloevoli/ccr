#include <iostream>

#include "constants.h"
#include "utilities.h"

void print_diffusion_time() {
    const auto z = 10.;
    const auto B_IGM = 1e-16 * cgs::Gauss;
    const auto n_H = n_H_physical(z);
    const auto d = cgs::Mpc;
    const auto j = dtdz(z);

    for (double E_k = 0.1 * cgs::MeV; E_k < 10. * cgs::GeV; E_k *= 1.1) {
        std::cout << std::scientific;
        std::cout << E_k / cgs::MeV << "\t";
        std::cout << hubble_time(z) / cgs::Gyr << "\t";
        const auto v = beta(E_k) * cgs::c_light;
        std::cout << (d / v) / cgs::Gyr << "\t";
        std::cout << diffusion_time(B_IGM, d, E_k) / cgs::Gyr << "\t";
        std::cout << inelastic_time(n_H, E_k) / cgs::Gyr << "\t";
        std::cout << -E_k / dEdz_H(z, E_k) * j / cgs::Gyr << "\t";
        std::cout << E_k / dEdt_i(n_H, E_k) / cgs::Gyr << "\t";
        std::cout << E_k / dEdt_C(n_H, E_k) / cgs::Gyr << "\t";
        std::cout << E_k / dEdt_C(1e-3 * n_H, E_k) / cgs::Gyr << "\t";
        std::cout << larmor_radius(B_IGM, E_k) / cgs::Mpc << "\t";
        std::cout << "\n";
    }
}

int main() {
    // std::cout << "Hello world!\n";

    print_diffusion_time();

    return 0;
}