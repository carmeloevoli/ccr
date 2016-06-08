#include "reionization.h"

Reionization::Reionization() {
    fast::init_ps();
    init_reionization();
    clumping_factor = 1;
    f_sfr = 0.1 / Gyr;
    f_alpha = 1 / mass_sun / pow3(Mpc);
}

Reionization::~Reionization() {}

void Reionization::init_reionization() {
    x_e = 1e-4;
    T_k = 1e4 * K;
    z = 30.;
}

void Reionization::evolve() {
    double dt = 0;
    double dx_e = 0;
    double dT_k = 0;
    
    while (z > 0) {
    
        dt = fast::dtdz(z) * dz;
        min_sfr_halo = min_star_forming_halo(z); // M
        mass_in_halos = integrate_hmf(z, min_sfr_halo, 1e12 * mass_sun); // M V^-1
        star_formation_rate = f_sfr * mass_in_halos; // M V^-1 T^-1
        ionization_rate = f_alpha * star_formation_rate; // T^-1
        recombination_rate = fast::alpha_A(T_k) * n_H_physical(z) * pow2(x_e) * clumping_factor; // L^3 T^-1 L^-3
        
        //dx_e = -dt * (Galaxy_ionization_rate(z) - fast::alpha_A(T_k) * n_H_physical(z) * pow2(x_e) * clumping_factor); // cm^3 s^-1 cm^-3
;
        dT_k = 0;
        x_e += dx_e;
        T_k += dT_k;
        z -= dz;
        
        print_status();
    }
}

void Reionization::print_status() {
    if (z > 0) {
        cout << scientific << setprecision(3);
        cout << z << "\t";
        cout << x_e << "\t";
        cout << T_k << "\t";
        cout << min_sfr_halo / mass_sun << "\t";
        cout << mass_in_halos / (mass_sun / pow3(Mpc)) << "\t";
        cout << star_formation_rate / (mass_sun / pow3(Mpc) / year) << "\t";
        cout << ionization_rate / (1. / Myr) << "\t";
        cout << recombination_rate / (1. / Myr) << "\t";
        cout << "\n";
    }
}

double Reionization::Galaxy_ionization_rate(const double& z) {
    //double mass_formed = integrate_hmf(z, min_star_forming_halo(z), 1e12 * mass_sun);
    
    //cout << z << "\t" << mass_formed / mass_sun << "\n";
    
    return 0;
}

double Reionization::Galaxy_heating_rate(const double& z) {
    return 0;
}

void Reionization::print_hmf(const double& z, const double& M_min, const double& M_max) {
    /*
     FUNCTION dNdM(z, M)
     Computes the Press_schechter mass function with Sheth-Torman correction for ellipsoidal collapse at
     redshift z, and dark matter halo mass M (in solar masses).
     
     The return value is the number density per unit mass of halos in the mass range M to M+dM in units of:
     comoving Mpc^-3 Msun^-1
     
     Reference: Sheth, Mo, Torman 2001
     */
    
    double dNdM = -1;
    
    printf ("#Redshift - Mass [Msun] - dNdM [Mpc^-3 Msun^-1]\n");
    
    for (double M = M_min; M <= M_max; M *= 1.1){
        dNdM = fast::dNdM_st(z,M);
        printf ("%5.2e %5.2e %5.2e\n", z, M, dNdM);
    }
}

double hmf_function (double M, void *params) {
    double z = *(double *) params;
    double dNdM = fast::dNdM_st(z, M / mass_sun) / pow3(Mpc) / mass_sun;
    //printf ("%5.2e %5.2e %5.2e\n", z, M, dNdM);
    return M * dNdM;
}

double Reionization::integrate_hmf(double z, const double& M_min, const double& M_max) {
    int KEYINTEGRAL = 2;
    int LIMIT = 10000;
    
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (LIMIT);
    
    double result, error;
    
    gsl_function F;
    F.function = &hmf_function;
    F.params = &z;
    
    gsl_integration_qag (&F, M_min, M_max, 0, 1e-4, LIMIT, KEYINTEGRAL,
                         w, &result, &error);
    
    //printf ("result          = % .8e\n", result);
    //printf ("estimated error = % .8e\n", error);
    
    gsl_integration_workspace_free (w);
    
    return result;
}