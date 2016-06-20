#include "SFR.h"

SFR::SFR(const string& filename) {
    this->filename = filename;
}

SFR::~SFR() {
}

void SFR::evolve() {
    outfile.open(filename.c_str());
    outfile << "#z - Min SFR halo mass - HMF integral" << "\n";
    for (double z = 0; z < 51; z += 0.1) {
        double min_sfr_halo = min_star_forming_halo(z); // M
        double hmf_integral = integrate_hmf(z, min_sfr_halo, 1e12 * mass_sun); // M V^-1 T^-1
        outfile << scientific << setprecision(15) << z << "\t" << min_sfr_halo << "\t" << hmf_integral << "\n";
        cout << z << "\n";
    }
    outfile.close();
}

double hmf_function (double M, void *params) {
    double z = *(double *) params;
    double dNdM = fast::dNdM_st(z, M / mass_sun) / pow3(Mpc) / mass_sun;
    //printf ("%5.2e %5.2e %5.2e\n", z, M, dNdM);
    return M * dNdM / free_fall_timescale(z, M);
}

double SFR::integrate_hmf(double z, const double& M_min, const double& M_max) {
    int KEYINTEGRAL = 3;
    int LIMIT = 10000;
    
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (LIMIT);
    
    double result, error;
    
    gsl_function F;
    F.function = &hmf_function;
    F.params = &z;
    
    gsl_integration_qag (&F, M_min, M_max, 0, 1e-5, LIMIT, KEYINTEGRAL,
                         w, &result, &error);
    
    //printf ("result          = % .8e\n", result);
    //printf ("estimated error = % .8e\n", error);
    
    gsl_integration_workspace_free (w);
    
    return result;
}

void SFR::print_hmf(const double& z, const double& M_min, const double& M_max) {
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
