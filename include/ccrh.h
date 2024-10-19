#ifndef CCRH_H_
#define CCRH_H_

#include <algorithm>
#include <cassert>
#include <numeric>
#include <vector>

#include "constants.h"
#include "io.h"
#include "tridiag.h"
#include "utilities.h"

class CCRH {
   public:
    CCRH() : Ek(ESize), Q(ESize), N(ESize, 0.0) {
        init_vectors();  // Initialize vectors in the constructor
    }

    ~CCRH() = default;

    void evolve();
    void dump();

   private:
    void build_energy_losses(std::vector<double>& energy_losses, double z, double x_e) const;
    double build_source_normalization(double z) const;
    std::pair<double, double> compute_rates(size_t iMin, double z) const;
    void CR_evolutor(double z, double dt);
    double CR_energy() const;

   private:
    std::vector<double> Ek;  // Energies
    std::vector<double> Q;   // Source spectrum
    std::vector<double> N;   // Particle number density

    static constexpr size_t ESize = 1000;
    static constexpr double EMin = 0.1 * cgs::MeV;
    static constexpr double EMax = 1e6 * cgs::MeV;

    static constexpr double z_min = 10.0;
    static constexpr double z_max = 20.0;
    static constexpr double dz = 0.001;

    static constexpr double SN_slope = 2.4;
    static constexpr double SN_evolution = 0.0;
    static constexpr double SN_L = 1e-33 * cgs::erg / cgs::cm3 / cgs::second;
    static constexpr double x_e = 1e-3;

    void init_vectors();  // Initializes energy vectors
};

#endif