#include <stdio.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#include "Parameter_files/INIT_PARAMS.H"
#include "Parameter_files/ANAL_PARAMS.H"
#include "Parameter_files/HEAT_PARAMS.H"

#define POW2(A) ((A)*(A))
#define POW3(A) ((A)*(A)*(A))

#define KEYINTEGRAL 2
#define LIMIT 10000

double hmf (double M, void * params) {
    double z = *(double *) params;
    double dNdM = dNdM_st(z,M);
    //printf ("%5.2e %5.2e %5.2e\n", z, M, dNdM);
    return dNdM;
}

double hmf_integral (double z)
{
    init_ps();

    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (LIMIT);
    
    double result, error;
    double M_min = 1e8 * pow(10. / (1.+z), 1.5);
    double M_max = 1e13;
    
    gsl_function F;
    F.function = &hmf;
    F.params = &z;
    
    gsl_integration_qag (&F, M_min, M_max, 0, 1e-4, LIMIT, KEYINTEGRAL,
                          w, &result, &error);
    
    //printf ("result          = % .8e\n", result);
    //printf ("estimated error = % .8e\n", error);
    
    gsl_integration_workspace_free (w);
    
    return result;
}

int main (int argc, char ** argv) {
    
    for (double z = 0.1; z < 20; z *= 1.1) {
        double result = hmf_integral(z);
        printf ("%5.2e %5.2e\n", z, result);
    }
        
    return 0;
}
