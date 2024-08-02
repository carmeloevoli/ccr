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
    
    for (double z = 20; z > 0; z -= 0.1){
        double M_min = 1e8 * pow(10. / (1. + z), 1.5);
        double dNdM_min = dNdM_st(z, M_min);
        printf ("%5.2e %5.2e %5.2e\n", z, M_min, dNdM);
    }
    
    return 0;
}

