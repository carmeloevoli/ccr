#include "utilities.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_erf.h>

#include <cmath>

double gamma_lorentz(double E_k) { return 1.0 + E_k / cgs::mass_proton_c2; }

double beta_lorentz(double E_k) {
    const double gamma = gamma_lorentz(E_k);
    const double gamma2 = pow2(gamma);
    return std::sqrt(1.0 - 1.0 / gamma2);
}

double larmor_radius(double B, double E_k) {
    const double p = E_k * std::sqrt(1.0 + 2.0 * cgs::mass_proton_c2 / E_k);
    constexpr double Z = 1.0;  // Assumed charge number
    return p / (Z * cgs::electron_charge * B);
}

double D_Bohm(double B, double E_k) { return larmor_radius(B, E_k) * cgs::c_light / 3.0; }

double diffusion_time(double B, double d, double E_k) { return pow2(d) / D_Bohm(B, E_k); }

double hubble_time(double z) {
    return 2.0 / (3.0 * cgs::H_0 * std::sqrt(cgs::Omega_m) * std::pow(1.0 + z, 1.5));
}

double n_H_physical(double z) { return cgs::n_H_0 * std::pow(1.0 + z, 3.0); }

double min_star_forming_halo(double z) {
    return 1e8 * cgs::mass_sun * std::pow(10.0 / (1.0 + z), 1.5);
}

double sigma_pp(double E_k) {
    constexpr double E_threshold = 0.2797 * cgs::GeV;
    const double x = E_k / E_threshold;
    if (x <= 1.0) {
        return 0.0;
    }
    double value = 30.7 - 0.96 * std::log(x) + 0.18 * pow2(std::log(x));
    return value * pow3(1.0 - std::pow(x, -1.9)) * cgs::mbarn;
}

double inelastic_time(double n_H, double E_k) {
    return 1.0 / (beta_lorentz(E_k) * cgs::c_light * n_H * sigma_pp(E_k));
}

double dtdz(double z) {
    const double x = std::sqrt(cgs::Omega_l / cgs::Omega_m) * std::pow(1.0 + z, -1.5);
    const double dxdz = -1.5 * x / (1.0 + z);
    const double numer = dxdz * (1.0 + x / std::sqrt(pow2(x) + 1.0));
    const double denom = x + std::sqrt(pow2(x) + 1.0);
    return (2.0 / (3.0 * cgs::H_0)) * numer / denom;
}

double dEdt_a(double z, double E_k) {
    double dEdz = E_k / (1.0 + z);
    dEdz *= 2.0 - 1.0 / (1.0 + cgs::mass_proton_c2 / E_k);
    return -dEdz / dtdz(z);
}

double Ionization_B(double IonPotential, double E_k) {
    const double beta = beta_lorentz(E_k);
    const double gamma = gamma_lorentz(E_k);
    const double mass_ratio = cgs::mass_electron / cgs::mass_proton;
    const double factor = 4.0 * pow4(gamma) * pow4(beta) * pow2(cgs::mass_electron_c2) /
                          (pow2(IonPotential) * (1.0 + 2.0 * gamma * mass_ratio));
    return std::log(factor - 2.0 * pow2(beta));
}

double dEdt_i(double n_HI, double E_k) {
    const double beta = beta_lorentz(E_k);
    const double factor = 3.0 * cgs::c_light * cgs::sigma_th * cgs::mass_electron_c2 / (4.0 * beta);
    const double B_H = Ionization_B(cgs::ionization_potential_H, E_k);
    const double B_He = Ionization_B(cgs::ionization_potential_He, E_k);
    const double B = B_H + cgs::f_He * B_He;
    return (B > 0.0) ? factor * n_HI * B : 0.0;
}

double Coulomb_ln(double n_e, double E_k) {
    const double factor = pow2(cgs::alpha) / (M_PI * pow3(cgs::electron_radius) * n_e);
    const double gamma = gamma_lorentz(E_k);
    return 0.5 * std::log(factor * pow2(gamma) * pow4(beta_lorentz(E_k)) /
                          (1.0 + 2.0 * gamma * (cgs::mass_electron / cgs::mass_proton)));
}

double Coulomb_W(double x) {
    return gsl_sf_erf(x) - (2.0 / std::sqrt(M_PI)) * x * std::exp(-pow2(x));
}

double dEdt_C(double n_e, double E_k) {
    const double factor = 1.5 * cgs::c_light * cgs::mass_electron_c2 * cgs::sigma_th;
    const double beta = beta_lorentz(E_k);
    const double beta_e = 2.0 * cgs::k_boltzmann * cgs::IGM_Temperature / cgs::mass_electron_c2;
    const double ln_C = Coulomb_ln(n_e, E_k);
    const double W_e = Coulomb_W(beta / beta_e);
    return factor * ln_C / beta * n_e * W_e;
}

double SN_Spectrum(double E_k, double alpha) {
    const double beta_0 = beta_lorentz(cgs::mass_proton_c2);
    const double beta = beta_lorentz(E_k);
    const double p = std::sqrt(pow2(E_k) + 2.0 * cgs::mass_proton_c2 * E_k);
    const double p_0 =
        std::sqrt(pow2(cgs::mass_proton_c2) + 2.0 * cgs::mass_proton_c2 * cgs::mass_proton_c2);
    return (beta_0 / beta) * std::pow(p / p_0, -alpha);
}

double SN_ESpectrum_integrand(double x, void* params) {
    const double alpha = *(double*)params;
    const double E = std::exp(x);
    return E * E * SN_Spectrum(E, alpha);
}

double SN_ESpectrum_integral(double alpha) {
    constexpr size_t N = 10000;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(N);
    double result = 0.0, error = 0.0;

    gsl_function F;
    F.function = &SN_ESpectrum_integrand;
    F.params = &alpha;

    const double Ek_min = 1e-6 * cgs::MeV;
    const double Ek_max = 1e6 * cgs::MeV;
    gsl_integration_qag(&F, std::log(Ek_min), std::log(Ek_max), 0, 1e-5, N, 3, w, &result, &error);

    gsl_integration_workspace_free(w);
    return result / pow2(cgs::mass_proton_c2);
}
