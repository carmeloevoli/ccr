#ifndef TRIDIAG_H_
#define TRIDIAG_H_

#include <iostream>
#include <vector>

using namespace std;

#define GSL_SUCCESS 0
#define GSL_ENOMEM 1
#define GSL_EZERODIV 2
#define GSL_EBADLEN 3

int gsl_linalg_solve_tridiag(const vector<double> & diag,
		const vector<double> & abovediag,
		const vector<double> & belowdiag,
		const vector<double> & rhs,
		vector<double> & solution);

int solve_tridiag_nonsym(const vector<double> & diag,
		const vector<double> & abovediag,
		const vector<double> & belowdiag,
		const vector<double> & rhs,
		vector<double> & x,
		size_t N);

void solveTridiagonalSystem(vector<double>& a,
		vector<double>& b,
		vector<double>& c,
		vector<double>& r,
		vector<double>& u,
		int n);

#endif /* TRIDIAG_H_ */
