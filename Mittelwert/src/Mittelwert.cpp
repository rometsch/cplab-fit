//============================================================================
// Name        : Mittelwert.cpp
// Author      : Thomas Rometsch
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
using namespace std;

int main() {
	vector<double> vel;
	vel.push_back(299793);
	vel.push_back(299792);
	vel.push_back(299782);
	vector<double> err_vel;
	err_vel.push_back(2);
	err_vel.push_back(4.5);
	err_vel.push_back(25);
	unsigned int L = vel.size();

	double mean = 0;
	double inv_err_sq_sum = 0;
	double Delta = 0;

	// Calculate weighted mean and inner variance.
	for (unsigned int i=0; i<L; i++) {
		Delta = err_vel[i]*err_vel[i];
		mean += vel[i]/Delta;
		inv_err_sq_sum += 1/Delta;
	}
	mean = mean/inv_err_sq_sum;

	// Convert given error on velocities to stddev and variance.
	// 95% confindence means err = 2*stddev, so devide by 2
	vector<double> stddev;
	stddev.reserve(L);
	vector<double> var;
	var.reserve(L);

	for (unsigned int i=0; i<L; i++) {
		stddev[i] = err_vel[i]/2;
		var[i] = stddev[i]*stddev[i];
	}

	// Calculate xi sq
	double xisq = 0;
	for (unsigned int i=0; i<L; i++) {
		double x = (vel[i] - mean)/err_vel[i];
		xisq += x*x;
	}

	double in_var = 0;
	for (unsigned int i=0; i<L; i++) {
		in_var += 1.0/var[i];
	}
	in_var = 1/in_var;

	double ext_var = xisq*in_var;

	cout 	<< scientific
			<< "mean = " << mean << endl
			<< "chisq = " << xisq << endl
			<< "innere varianz = " << in_var << endl
			<< "externe varianz = " << ext_var << endl;

	return 0;
}
