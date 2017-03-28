//============================================================================
// Name        : LineareReg.cpp
// Author      : Thomas Rometsch
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

vector<double> fit_linear_function(vector<double> x, vector<double> y, vector<double> y_err);
vector<double> fit_proportionality(vector<double> x, vector<double> y, vector<double> y_err);


int main() {
	// Spannungs als X-Werte
	vector<double> x = {0.5, 1,	1.5, 2,	2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6};
	// Strom als Y-Werte
	vector<double> y = {0.065, 0.206, 0.405, 0.492, 0.606, 0.782, 0.865, 1.018, 1.199, 1.327, 1.408, 1.627};
	vector<double> err = {0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.04};

	{
		vector<double> fit_parameters = fit_linear_function(x,y,err);

		double a = fit_parameters[0];
		double var_a = fit_parameters[2];
		double b = fit_parameters[1];
		double var_b = fit_parameters[3];
		double chisq = fit_parameters[4];

		double R = 1/b;
		double stddev_R = sqrt(var_b)/b/b;
		double I0 = a;
		double stddev_I0 = sqrt(var_a);

		cout 	<< "====================================" << endl
				<< "====================================" << endl
				<< "Fit: y = a + bx" << endl
				<< scientific
				<< "a = " << a << endl
				<< "Var(a) = " << var_a << endl
				<< "b = " << b << endl
				<< "Var(b) = " << var_b << endl
				<< "chisq/(N-2) = " << chisq << endl
				<< "====================================" << endl
				<< "R = " << R << endl
				<< "stddev_R = " << stddev_R << endl
				<< "I0 = " << I0 << endl
				<< "stddev_I0 = " << stddev_I0 << endl
				<< "====================================" << endl
				<< "====================================" << endl;
	}
	cout << endl;

	{
		vector<double> fit_parameters = fit_proportionality(x,y,err);

		double b = fit_parameters[0];
		double var_b = fit_parameters[1];
		double chisq = fit_parameters[2];

		double R = 1/b;
		double stddev_R = sqrt(var_b)/b/b;

		cout 	<< "====================================" << endl
				<< "====================================" << endl
				<< "Fit: y = bx" << endl
				<< scientific
				<< "b = " << b << endl
				<< "Var(b) = " << var_b << endl
				<< "chisq/(N-2) = " << chisq << endl
				<< "====================================" << endl
				<< "R = " << R << endl
				<< "stddev_R = " << stddev_R << endl
				<< "====================================" << endl
				<< "====================================" << endl;
	}
}

vector<double> fit_linear_function(vector<double> x, vector<double> y, vector<double> y_err) {
	/* Lineare Regression der Funktion y = f(x) = a*x + b
	 * Zurueckgegeben wird ein Vector mit den Eintraegen (a,b,var_a,var_b,chisq)
	 */
	unsigned int N = x.size();

	// Hilfsgrößen berechnen
	double S = 0;
	double Sx = 0;
	double Sy = 0;
	double Sxx = 0;
	double Sxy = 0;

	for (unsigned int i=0; i<N; i++) {
		double invvar = 1.0/(y_err[i]*y_err[i]);
		S += invvar;
		Sx += x[i]*invvar;
		Sy += y[i]*invvar;
		Sxx += x[i]*x[i]*invvar;
		Sxy += x[i]*y[i]*invvar;
	}
	double Det = S*Sxx-Sx*Sx;

	double a = (Sxx*Sy - Sx*Sxy)/Det;
	double b = (S*Sxy-Sx*Sy)/Det;
	double var_a = Sxx/Det;
	double var_b = S/Det;

	double chisq = 0;
	for (unsigned int i=0; i<N; i++) {
		double X = (y[i] - a - b*x[i])/y_err[i];
		chisq += X*X;
	}
	chisq = chisq/(N-2);

	vector<double> res;
	res.push_back(a);
	res.push_back(b);
	res.push_back(var_a);
	res.push_back(var_b);
	res.push_back(chisq);

	return res;
}

vector<double> fit_proportionality(vector<double> x, vector<double> y, vector<double> y_err) {
	/* Lineare Regression der Funktion y = f(x) = a*x
	 * Zurueckgegeben wird ein Vector mit den Eintraegen (a,var_a,chisq)
	 */
	unsigned int N = x.size();

	// Hilfsgrößen berechnen
	double S = 0;
	double Sx = 0;
	double Sy = 0;
	double Sxx = 0;
	double Sxy = 0;

	for (unsigned int i=0; i<N; i++) {
		double invvar = 1.0/(y_err[i]*y_err[i]);
		S += invvar;
		Sx += x[i]*invvar;
		Sy += y[i]*invvar;
		Sxx += x[i]*x[i]*invvar;
		Sxy += x[i]*y[i]*invvar;
	}

	double a = Sxy/Sxx;
	double var_a = 1.0/Sxx;

	double chisq = 0;
	for (unsigned int i=0; i<N; i++) {
		double X = (y[i] - a*x[i])/y_err[i];
		chisq += X*X;
	}
	chisq = chisq/(N-2);

	vector<double> res;
	res.push_back(a);
	res.push_back(var_a);
	res.push_back(chisq);

	return res;
}
