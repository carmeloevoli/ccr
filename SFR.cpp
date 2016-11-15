#include "SFR.h"

SFR::SFR(const string& filename) {
	this->filename = filename;
	outfile.open(filename.c_str());
	vvir_cut = 24. * km / s;
	efficiency = 1.0;
}

SFR::~SFR() {
	outfile.close();
}

void SFR::evolve() {
	outfile << "#z - Min SFR halo mass - HMF integral" << "\n";
    int i = 0;
	for (double z = 0; z < 51; z += 0.1) {
		double min_sfr_halo = min_star_forming_halo(z); // M
		double min_filtering_halo = halo_mass_given_vvir(vvir_cut, z);
		double hmf_integral = integrate_hmf(z, max(min_sfr_halo, min_filtering_halo), 1e13 * mass_sun); // M V^-1 T^-1
		outfile << scientific << setprecision(15) << z << "\t" << min_sfr_halo << "\t" << hmf_integral << "\n";
        if (i % 10 == 0) cout << z << "\n";
        i++;
	}
}

double hmf_function (double logM, void *params) {
	double M = exp(logM);
	double z = *(double *) params;
	double dNdM = fast::dNdM_st(z, M / mass_sun) / pow3(Mpc) / mass_sun;
	double tff = free_fall_timescale(z, M);
	//printf ("%5.2e %5.2e %5.2e\n", z, M, dNdM);
	return M * M * dNdM / tff;
}

double SFR::integrate_hmf(double z, const double& M_min, const double& M_max) {
	int KEYINTEGRAL = 4;
	int LIMIT = 30000;

	gsl_integration_workspace * w
	= gsl_integration_workspace_alloc (LIMIT);

	double result, error;

	gsl_function F;
	F.function = &hmf_function;
	F.params = &z;

	gsl_integration_qag (&F, log(M_min), log(M_max), 0, 1e-4, LIMIT, KEYINTEGRAL,
			w, &result, &error);

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
	double ff = -1;

	printf ("#Redshift - Mass [Msun] - dNdM [Mpc^-3 Msun^-1] - t_ff [s] \n");

	for (double M = M_min; M <= M_max; M *= 1.1){
		dNdM = fast::dNdM_st(z, M);
		ff = free_fall_timescale(z, M);
		printf ("%5.2e %5.2e %5.2e %5.2e\n", z, M, dNdM, ff);
	}
}
