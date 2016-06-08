#ifndef _MISC_
#define _MISC_

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

/*** Some usefull math macros ***/
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//static double mnarg1,mnarg2;
//#define FMAX(a,b) (mnarg1=(a),mnarg2=(b),(mnarg1) > (mnarg2) ?\
//(mnarg1) : (mnarg2))

//static double mnarg1,mnarg2;
//#define FMIN(a,b) (mnarg1=(a),mnarg2=(b),(mnarg1) < (mnarg2) ?\
//(mnarg1) : (mnarg2))

/*** Wrapper functions for the std library functions fwrite and fread
 which should improve stability on certain 64-bit operating systems
 when dealing with large (>4GB) files ***/
size_t mod_fwrite (const void *, unsigned long long, unsigned long long, FILE *);
size_t mod_fread(void *, unsigned long long, unsigned long long, FILE *);

/* generic function to compare floats */
int compare_floats(const void *, const void *);

/* fast complimentary error function approximation */
float erfcc(float);

double FMAX(double a, double b);

#endif