#include "utilities.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_erf.h>

#include <iostream>

double gamma_lorentz(double E_k) { return 1. + E_k / cgs::mass_proton_c2; }

double beta_lorentz(double E_k) {
    const auto gamma2 = pow2(gamma_lorentz(E_k));
    return std::sqrt(1.0 - 1.0 / gamma2);
}

double larmor_radius(double B, double E_k) {
    const auto p = E_k * sqrt(1. + 2. * cgs::mass_proton_c2 / E_k);
    const auto Z = 1.;
    return p / Z / cgs::electron_charge / B;
}

double D_Bohm(double B, double E_k) {
    const auto rL = larmor_radius(B, E_k);
    return rL * cgs::c_light / 3.;
}

double diffusion_time(double B, double d, double E_k) {
    const auto DB = D_Bohm(B, E_k);
    return d * d / DB;
}

double hubble_time(double z) {
    return 2. / 3. / cgs::H_0 / std::sqrt(cgs::Omega_m) * std::pow(1. + z, -1.5);
}

double n_H_physical(double z) { return cgs::n_H_0 * std::pow(1. + z, 3.); }

double min_star_forming_halo(double z) {
    return 1e8 * cgs::mass_sun * std::pow(10. / (1. + z), 1.5);
}

double sigma_pp(double E_k) {  // Kafexhiu et al. PRD 90, 2014
    const auto E_threshold = 0.2797 * cgs::GeV;
    const auto x = E_k / E_threshold;
    double value = 0;
    if (x > 1) {
        value = 30.7 - 0.96 * log(x) + 0.18 * pow2(log(x));
        value *= pow3(1 - pow(x, -1.9));
    }
    return value * cgs::mbarn;
}

double inelastic_time(double n_H, double E_k) {
    const auto sigma = sigma_pp(E_k);
    return 1. / beta_lorentz(E_k) / cgs::c_light / n_H / sigma;
}

/* function DTDZ returns the value of dt/dz at the redshift parameter z. */
double dtdz(double z) {
    using std::pow;
    using std::sqrt;

    const auto x = sqrt(cgs::Omega_l / cgs::Omega_m) * pow(1 + z, -3.0 / 2.0);
    const auto dxdz = sqrt(cgs::Omega_l / cgs::Omega_m) * pow(1 + z, -5.0 / 2.0) * (-3.0 / 2.0);
    const auto const1 = 2 * sqrt(1 + cgs::Omega_m / cgs::Omega_l) / (3.0 * cgs::H_0);

    const auto numer = dxdz * (1 + x * pow(pow(x, 2) + 1, -0.5));
    const auto denom = x + sqrt(pow(x, 2) + 1);
    return (const1 * numer / denom);
}

double dEdt_a(double z, double E_k) {
    auto dEdz = E_k / (1. + z);
    dEdz *= 2. - 1. / (1. + cgs::mass_proton_c2 / E_k);
    return -dEdz / dtdz(z);
}

double Ionization_B(double IonPotential, double E_k) {
    const auto beta = beta_lorentz(E_k);
    const auto beta2 = pow2(beta);
    const auto beta4 = pow4(beta);
    const auto gamma = gamma_lorentz(E_k);
    const auto gamma4 = pow4(gamma);
    const auto mass_ratio = cgs::mass_electron / cgs::mass_proton;
    const auto me2 = pow2(cgs::mass_electron_c2);
    const auto I2 = pow2(IonPotential);
    const auto factor = 4. * me2 * beta4 * gamma4 / I2 / (1. + 2 * gamma * mass_ratio);
    return log(factor - 2. * beta2);
}

double dEdt_i(double n_HI, double E_k) {
    const auto beta = beta_lorentz(E_k);
    const auto factor = 3. * cgs::c_light * cgs::sigma_th * cgs::mass_electron_c2 / 4. / beta;
    const auto B_H = Ionization_B(cgs::ionization_potential_H, E_k);
    const auto B_He = Ionization_B(cgs::ionization_potential_He, E_k);
    const auto B = B_H + cgs::f_He * B_He;
    return (B > 0.) ? factor * n_HI * B : 0.;
}

double Coulomb_ln(double n_e, double E_k) {
    const auto factor = pow2(cgs::alpha) / M_PI / pow3(cgs::electron_radius) / n_e;
    const auto beta = beta_lorentz(E_k);
    const auto beta4 = pow4(beta);
    const auto gamma = gamma_lorentz(E_k);
    const auto gamma2 = pow2(gamma);
    const auto mass_ratio = cgs::mass_electron / cgs::mass_proton;
    return .5 * log(factor * gamma2 * beta4 / (1. + 2. * gamma * mass_ratio));
}

double Coulomb_W(double x) {
    const auto value_a = gsl_sf_erf(x);
    const auto value_b = 2. / sqrt(M_PI) * x * exp(-pow2(x));
    return value_a - value_b;
}

double dEdt_C(double n_e, double E_k) {
    const auto factor = 3. * cgs::c_light * cgs::mass_electron_c2 * cgs::sigma_th / 2.;
    const auto beta = beta_lorentz(E_k);
    const auto beta_e = 2. * cgs::k_boltzmann * cgs::IGM_Temperature / cgs::mass_electron_c2;
    const auto ln_C = Coulomb_ln(n_e, E_k);
    const auto W_e = Coulomb_W(beta / beta_e);
    return factor * ln_C / beta * n_e * W_e;
}

double SN_Spectrum(double Ek, double alpha) {
    const auto E0 = cgs::mass_proton_c2;
    const auto beta = beta_lorentz(Ek);
    const auto beta_0 = beta_lorentz(E0);
    const auto p = sqrt(pow2(Ek) + 2. * cgs::mass_proton_c2 * Ek);
    const auto p_0 = sqrt(pow2(E0) + 2. * cgs::mass_proton_c2 * E0);
    return beta_0 / beta * pow(p / p_0, -alpha);
}

double SN_ESpectrum_integrand(double x, void* params) {
    double alpha = *(double*)params;
    double E = exp(x);
    return E * E * SN_Spectrum(E, alpha);
}

double SN_ESpectrum_integral(double alpha) {
    size_t N = 10000;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(N);

    double result, error;
    gsl_function F;
    F.function = &SN_ESpectrum_integrand;
    F.params = &alpha;

    const auto Ek_min = 1e-6 * cgs::MeV;
    const auto Ek_max = 1e6 * cgs::MeV;
    const auto E0 = cgs::mass_proton_c2;

    gsl_integration_qag(&F, log(Ek_min), log(Ek_max), 0, 1e-5, N, 3, w, &result, &error);

    gsl_integration_workspace_free(w);

    return result / pow2(E0);
}