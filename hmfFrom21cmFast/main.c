#include <stdio.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#include "Parameter_files/INIT_PARAMS.H"
#include "Parameter_files/ANAL_PARAMS.H"
#include "Parameter_files/HEAT_PARAMS.H"

#define POW2(A) ((A)*(A))
#define POW3(A) ((A)*(A)*(A))

int main(int argc, char ** argv)
{   
    init_ps();
    
    /*
     FUNCTION dNdM(z, M)
     Computes the Press_schechter mass function with Sheth-Torman correction for ellipsoidal collapse at
     redshift z, and dark matter halo mass M (in solar masses).
     
     The return value is the number density per unit mass of halos in the mass range M to M+dM in units of:
     comoving Mpc^-3 Msun^-1
     
     Reference: Sheth, Mo, Torman 2001
     */
        
    printf ("#Redshift - Mass [Msun] - dNdM [Mpc^-3 Msun^-1]\n");
    
    for (double M = 1e6; M < 1e12; M *= 1.05){
        double dNdM_0 = dNdM_st(10., M);
        double dNdM_6 = dNdM_st(6., M);
        double dNdM_3 = dNdM_st(3., M);
        double dNdM_1 = dNdM_st(1., M);
        printf ("%6.3e %6.3e %6.3e %6.3e %6.3e\n", M, dNdM_0, dNdM_6, dNdM_3, dNdM_1);
    }
    
    return 0;
}
