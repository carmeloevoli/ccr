#ifndef _PS_
#define _PS_

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "constants.h"
#include "cosmo_progs.h"
#include "misc.h"

/* New in v1.1 */
#define ERFC_NPTS (int) 75
#define ERFC_PARAM_DELTA (float) 0.1
static double log_erfc_table[ERFC_NPTS], erfc_params[ERFC_NPTS];
static gsl_interp_accel *erfc_acc;
static gsl_spline *erfc_spline;

extern double sigma_norm, R, theta_cmb, omhh, z_equality, y_d, sound_horizon, alpha_nu, f_nu, f_baryon, beta_c, d2fact, R_CUTOFF, DEL_CURR, SIG_CURR;

namespace fast {
    double init_ps(); /* initialize global variables, MUST CALL THIS FIRST!!! returns R_CUTOFF */
    void free_ps(); /* deallocates the gsl structures from init_ps */
    double splined_erfc(double); /* returns erfc for x>=0, using cubic spline in logy-x space */
    double deltolindel(float del, float z); /* converts a non-linear overdensity, del, at z to a linear overdensity at z=0 */
    double lindeltodel(float lindel, float z); /* converts a linear overdensity, del, at z=0 to a non-linear overdensity at redshift z */
    double power_in_k(double k); /* Returns the value of the linear power spectrum density (i.e. <|delta_k|^2>/V) at a given k mode at z=0 */
    double RtoM(double); /* R in Mpc, M in Msun */
    double MtoR(double); /* R in Mpc, M in Msun */
    double M_J_WDM(); /* returns the "effective Jeans mass" corresponding to the gas analog of WDM ; eq. 10 in BHO 2001 */
    double sheth_delc(double del, double sig);
    double dNdM_st(double z, double M);
    double dNdM(double z, double M);
    double dnbiasdM(double M, float z, double M_o, float del_o); /* dnbiasdM */
    double FgtrM(double z, double M);  //calculates the fraction of mass contained in haloes with mass > M at redshift z
    double FgtrM_st(double z, double M);  //calculates the fraction of mass contained in haloes with mass > M at redshift z, with Sheth-Tormen correction
    double FgtrM_bias(double z, double M, double del_bias, double sig_bias);  //calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias
    double sigmaparam_FgtrM_bias(float z, float sigsmallR, float del_bias, float sig_bias);/* Uses sigma parameters instead of Mass for scale */
    double FgtrM_bias_BL08(double z, double M, double del_bias, double sig_bias); // as above, but this version uses the hybrid perscription of Barkana & Loeb 2004 (specifically the separate integral version of eq. 2 in Barkana & Loeb 2008)
    double dicke(double z); //calculates the dicke growth function at redshift z
    double ddickedz(double z); /* Redshift derivative of the growth function at z */
    double ddickedt(double z); /* Time derivative of the growth function at z */
    double sigma_z0(double M); //calculates sigma at z=0 (no dicke)
    double dsigmasqdm_z0(double M); //calculates d(sigma^2)/dm at z=0 (i.e. does not include dicke growth)
    double TFmdm(double k); //Eisenstien & Hu power spectrum transfer function
    void TFset_parameters();
    float get_R_c();  // returns R_CUTOFF
    double get_M_min_ion(float z);
}
#endif