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
    void print_hmf(const double& z, const double& M_min, const double& M_max);

	void set_efficiency(double efficiency) {
		this->efficiency = efficiency;
	}

	void set_filename(const string& filename) {
		this->filename = filename;
	}

	void set_vvir_cut(double vvir_cut) {
		this->vvir_cut = vvir_cut;
	}

private:
    double integrate_hmf(double z, const double& M_min, const double& M_max);
    double efficiency;
    double vvir_cut;
    string filename;
    ofstream outfile;
};

#endif
