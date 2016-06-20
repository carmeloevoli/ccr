#ifndef __SFR_H
#define __SFR_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "ps.h"
#include "utilities.h"

using namespace std;

class SFR {
public:
    SFR(const string& filename);
    ~SFR();
    
    void evolve();
    
private:
    double integrate_hmf(double z, const double& M_min, const double& M_max);
    void print_hmf(const double& z, const double& M_min, const double& M_max);

    string filename;
    ofstream outfile;
};

#endif