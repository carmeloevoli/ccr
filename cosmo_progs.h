/*
 Usefull cosmo progs.
 CAUTION:  many of these assume standard cosmology (k=0, etc.), so check them before doing wierd things to cosmology.
 */

#ifndef _COSMO_PROGS_
#define _COSMO_PROGS_

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include <gsl/gsl_integration.h>

#include "constants.h"

namespace fast {
    
    double hubble(float z);  /* returns the hubble "constant" at z */
    double t_hubble(float z);  /* returns hubble time, t_h = 1/H */
    float t_dynamical(float z); /* dynamical time at z in seconds */
    double M_jeans_neutralIGM(float z); /* returns the cosmological jeans mass at z; in neutral IGM */
    double M_jeans_general(float z, float Del, float T, float mu); /* returns the cosmological jeans mass (in solar masses) at z, non-linear overdensity Del=rho/<rho>, temperature T, and \mu */
    float z_gasCMBdecoupling(); /* the redshift at which the gas temperature diverges from the CMB temperature */
    double neutral_fraction(double density, double gamma_ion); /* neutral fraction given H density (cm^-3) and gamma (1/s) */
    double nftogamma(double nf, float z); /* ionization rate (1/s) from neutral fraction (assuming mean density at z) */
    double alpha_A(double T); //case A hydrogen recombination coefficient (Abel et al. 1997)
    float TtoM(float z, float T, float mu);
    float MtoVcir(float z, double M); // M in M_sun, Vcir in proper km/s
    float MtoRvir(float z, double M); // M in M_sun, Rvir in comoving Mpc
    float VcirtoT(float v, float mu);
    double dVdz(float z);
    double arcmintoMpc(float z, float arcmin); /* arcmin to comoving Mpc */
    double d_A(float z); /* angular diameter distance (proper Mpc) */
    double d_L(float z); /* luminosity distance (proper Mpc) */
    double gettime(double z);
    double timesince(double z, double zsource);
    double rtoz (double R, double zsource); /* in comoving Mpc */
    double properdistance(double z, double zsource); /* in cm */
    double comovingdistance(double z, double zsource); /* in cm */
    double ttoz(double t, double zsource);
    double invsinh(double x);
    double dtdz(float z);
    double drdz(float z); /* comoving distance, (1+z)*C*dtdz(in cm) per unit z */
    double Lya_crosssec(double nu, double T); /* the differential Lya absorption crosssection (in cm) at frequency nu, and gas temperature T */
    double Lya_dopwidth(double T);
    double voigt1(double x, double a );
    double HI_ion_crosssec(double nu);
    double HeII_ion_crosssec(double nu);
    double HeI_ion_crosssec(double nu);
    double omega_mz(float z); /* Omega_m at redshift z*/
    double rho_critz(float z); /* critical density at redshift z in  Msun Mpc^-3 */
    double Deltac_nonlinear(float z);
    double anal_taudamp(float nuobs, float zstart, float zend); /* Lya damping wing optical depth contributed by a neutral IGM from zstart to zend (zstart > zend) at observed frequency nuobs */
    double mean_rho(float z); /* returns the proper mean baryonic density at z in g/cm^3 */
    double T_gas_adiabatic(float z); /* returns the adiabatically-cooled gas temperature at z */
    double tau_e(float zstart, float zend, float *zarry, float *xiarry, int len); /* returns the thompson scattering optical depth from zstart to zend.  The ionization history is specified by the last three arguements.  Set these to NULL or 0 if fully ionized IGM is assumed. */
}

#endif  /* end _COSMO_PROGS_ */

