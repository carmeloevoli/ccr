#include "ccrh.h"

// Initializes energy vector `Ek` and source spectrum `Q`
void CCRH::init_vectors() {
    const double deltaLogE = std::exp(std::log(EMax / EMin) / (ESize - 1));

    // Use std::generate to initialize `Ek` and `Q`
    std::generate(Ek.begin(), Ek.end(), [&, i = 0]() mutable {
        double Ek_i = std::exp(std::log(EMin) + i * std::log(deltaLogE));
        ++i;
        return Ek_i;
    });

    std::transform(Ek.begin(), Ek.end(), Q.begin(),
                   [](double E) { return SN_Spectrum(E, SN_slope); });
}

// Builds the energy loss rates for a given redshift `z` and ionization fraction `x_e` in-place
void CCRH::build_energy_losses(std::vector<double>& energy_losses, double z, double x_e) const {
    const double n_H = n_H_physical(z);
    const double n_e = x_e * n_H;
    const double n_HI = (1.0 - x_e) * n_H;

    // Reuse the existing vector and compute in-place to avoid unnecessary allocations
    std::transform(Ek.begin(), Ek.end(), energy_losses.begin(),
                   [&](double E) { return dEdt_a(z, E) + dEdt_C(n_e, E) + dEdt_i(n_HI, E); });
}

// Computes the source normalization factor based on redshift `z`
double CCRH::build_source_normalization(double z) const {
    constexpr double E_0 = cgs::mass_proton_c2;
    const double L_z = SN_L * std::pow(1.0 + z, SN_evolution);
    return L_z / pow2(E_0) / SN_ESpectrum_integral(SN_slope);
}

// Evolves the cosmic ray population over a timestep `dt` at redshift `z`
void CCRH::CR_evolutor(double z, double dt) {
    if (dt <= 0.0) throw std::invalid_argument("Time step must be greater than zero.");

    std::vector<double> energy_losses(ESize);
    build_energy_losses(energy_losses, z, x_e);  // Build energy losses in-place
    const double Q0 = build_source_normalization(z);

    std::vector<double> diagonal(ESize - 1);
    std::vector<double> upperDiagonal(ESize - 2);
    std::vector<double> lowerDiagonal(ESize - 2);
    std::vector<double> knownTerm(ESize - 1);
    std::vector<double> N_next(ESize - 1);

    for (size_t i = 0; i < ESize - 1; ++i) {
        const auto alpha = 0.5 * dt / (Ek[i + 1] - Ek[i]);
        diagonal[i] = 1.0 + alpha * energy_losses[i];
        if (i != ESize - 2) {
            upperDiagonal[i] = -alpha * energy_losses[i + 1];
        }
        if (i != 0) {
            lowerDiagonal[i - 1] = 0.;
        }
        knownTerm[i] = (1.0 - alpha * energy_losses[i]) * N[i];
        knownTerm[i] += alpha * energy_losses[i + 1] * N[i + 1];
        knownTerm[i] += dt * Q0 * Q[i];
    }

    // Solve the tridiagonal system for the next step of N
    gsl_linalg_solve_tridiag(diagonal, upperDiagonal, lowerDiagonal, knownTerm, N_next);

    // Ensure non-negative particle numbers
    std::transform(N_next.begin(), N_next.end(), N.begin(),
                   [](double val) { return std::max(val, 0.0); });
}

// Evolves the cosmic ray population over the redshift range
void CCRH::evolve() {
    double z = z_max;
    while (z > z_min) {
        const double dt = -dtdz(z) * dz;
        CR_evolutor(z, dt);
        z -= dz;
        std::cout << "z: " << z << ", CR Energy: " << CR_energy() / cgs::erg << " erg\n";
    }
}

// Calculates the total cosmic ray energy
double CCRH::CR_energy() const {
    const double dlnE = std::log(Ek[1] / Ek[0]);

    double sum = 0;
    for (size_t i = 0; i < ESize; ++i) {
        const double Ek_i = Ek[i];
        sum += pow2(Ek_i) * N[i];
    }

    return dlnE * sum;
}

// Calculates the cosmic ray ionization rate
std::pair<double, double> CCRH::compute_rates(size_t iMin, double z) const {
    const double dlnE = std::log(Ek[1] / Ek[0]);
    const double n_H = n_H_physical(z);
    const double n_HI = (1.0 - x_e) * n_H;
    const double n_e = x_e * n_H;
    const double W_H = 36.3 * cgs::eV;

    double sum_i = 0;
    double sum_C = 0;
    for (size_t i = iMin; i < ESize; ++i) {
        const double Ek_i = Ek[i];
        sum_i += dEdt_i(n_HI, Ek_i) * Ek_i * N[i];
        sum_C += dEdt_C(n_e, Ek_i) * Ek_i * N[i];
    }

    return {dlnE * sum_i / W_H / cgs::n_H_0, dlnE * sum_C / cgs::n_H_0};
}

// Dumps the current state of the system to a file
void CCRH::dump() {
    const double t_H = hubble_time(z_min);
    const double Q_0 = build_source_normalization(z_min);

    std::ofstream outfile = open_output_file("ccrh_solution.txt");
    for (size_t i = 0; i < ESize; ++i) {
        outfile << std::scientific << std::setprecision(6);
        outfile << Ek[i] / cgs::MeV << "\t";
        outfile << N[i] << "\t";
        outfile << Q_0 * Q[i] * t_H << "\t";
        auto rates = compute_rates(i, z_min);
        outfile << rates.first << "\t";
        outfile << rates.second << "\t";
        outfile << "\n";
    }
    outfile.close();
}