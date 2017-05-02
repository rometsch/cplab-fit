// Beispiel f�r Aufruf von PIKAIA mit C++-Hauptprogramm
// und C++-Routinen f�r zu optimierende Funktionen

// �bersetzung der Programmteile:

// f77 -c fortran_part.f oder
// gfortran -c fortran_part.f  (je nach FORTRAN-Compiler)
// g++ fortran_part.o cpp_part.cpp -lm -lf2c -o programm  bzw.
// g++ fortran_part.o cpp_part.cpp -lm -lgfortran -o programm
// Definition der FORTRAN-Parameter: siehe PIKAIA-Manual

#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <vector>

using namespace std;

extern "C"  // Definition der FORTRAN-Routinen pikaia und rninit
{
    void pikaia_(float(int*, float*), int&, float*, float*, float*, int*);
    void rninit_(int&);
};

// Implement a class to handle ChiSq calculation which can be used with PIKAIA
int Npts;
vector<float> X,Y,Sigma;
float scale_max [5] = {10, 100, 100, 50, 1};
float scale_min [5] = {0,0,0,1,0};

float scale_parameter(float pikaia_param, int ind) {
	return scale_min[ind] + pikaia_param*(scale_max[ind] - scale_min[ind]);
}

float func(int Nparams, float *params, float x) {
	// fit function taking parameters from PIKAIA
	// first two parameters in p are a and b in : f(x) = a*x + b
	const float PI = acos(0.);
	//Scaling
	int Nsins = (Nparams-2)/3;
	vector<float> amp;
	vector<float> freq;
	vector<float> phase;
	float a = scale_parameter(params[0],0);
	float b = scale_parameter(params[1],1);
	for (int i=0; i<Nsins; i++) {
		amp.push_back(scale_parameter(params[2+3*i],2));
		freq.push_back(scale_parameter(params[3+3*i],3));
		phase.push_back(scale_parameter(params[4+3*i],4));
	}
	float sum = a*x+b;
	for (int i=0; i<Nsins; i++) {
		sum += amp[i]*sin( 2*PI*(x/freq[i] + phase[i]) );
	}
	return sum;
}

float calc_ChiSq(int Nparams, float *params) {
	// n = number of parameters
	// x = array of parameters
	float ChiSq = 0;
	for (int i=0; i<Npts; i++) {
		ChiSq += pow( (Y[i] - func(Nparams,params,X[i]) )/Sigma[i] , 2);
	}
	return ChiSq;
}

static float fitness(int *n, float *params) {
	// Calc fitness as 1/ChiSq
	int Nparams = *n;
	return 1./calc_ChiSq(Nparams, params);
}

void fit_data_lin_with_sins(int Nsins, bool verbose) {
	// Initialisierung
	float f; // Fitness
	std::srand(time(NULL));
	int status; // status
	int Nparams = 2+3*Nsins;   // Zahl der Parameter
	int rn=rand();
	float ctrl[]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0};  // Steuerungsfeld
	rninit_(rn);
	float *params = new float[Nparams]; // 2 Parameter
	// Aufruf von PIKAIA
	if(verbose) {
		cout 	<< "Fitting data with a linear and sine functions" << endl
			<< "Number of sine functions = " << Nsins << endl;
	}
	pikaia_(fitness,Nparams,ctrl,params,&f,&status);
	// Rescale params
	vector<float> amp;
	vector<float> freq;
	vector<float> phase;
	float a = scale_parameter(params[0],0);
	float b = scale_parameter(params[1],1);
	for (int i=0; i<Nsins; i++) {
		amp.push_back(scale_parameter(params[2+3*i],2));
		freq.push_back(scale_parameter(params[3+3*i],3));
		phase.push_back(scale_parameter(params[4+3*i],4));
	}
	// Auswerten
	if (verbose) {
		cout 	<< "status = " << status << endl
				<< "ChiSq/(Npts - Nparam) = " << calc_ChiSq(Nparams,params) << endl
				<< "a = " << a << endl
				<< "b = " << b << endl;

		for (int i=0; i<Nsins; i++) {
			cout
				<< "sine numer " << i << endl
				<< "A = " << amp[i] << endl
				<< "P = " << freq[i] << endl
				<< "Phi = " << phase[i] << endl;
		}
	}
	// function to use with gnuplot
	cout << "chisq_reduced_" << Nsins << " = " << calc_ChiSq(Nparams,params) << endl;
	cout 	<< "sin" << Nsins << "(x) = "
			<< a << "*x + "
			<< b;
	for (int i=0; i<Nsins; i++) {
		cout
			<< " + "
			<< amp[i] << "*"
			<< "sin( 2*pi*( x/" << freq[i] << " + "
			<< phase[i] << ") ) ";
	}
	cout << endl;
}

float bsp_func(int *n, float *x){  // Beispielfunktion mit 2 Parametern
    const float PI = acos(0.);     // wird aus PIKAIA heraus aufgerufen
    return x[0]*x[1];
}

float  P1(int *n, float *x) {
	const float PI = acos(0.);
	float rsq = (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5);
	int m = 9;
	float sigsq = 0.15;
	float c = cos(PI*sqrt(rsq)*m);
	return c*c * exp(-rsq/sigsq);
}

float P2(int *n, float *x) {
	float r1sq = (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5);
	float r2sq = (x[0]-0.6)*(x[0]-0.6) + (x[1]-0.1)*(x[1]-0.1);
	float sig1sq = 0.3*0.3;
	float sig2sq = 0.03*0.03;
	return 0.8*exp(-r1sq/sig1sq) + 0.879008*exp(-r2sq/sig2sq);
}

void maximize_function(float func(int*,float*), string comment) {
	// Initialisierung
	float f; // Fitness
	std::srand(time(NULL));
	int status; // status
	int n = 2;   // Zahl der Parameter
	int rn=rand();
	float ctrl[]={-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,0};  // Steuerungsfeld
	rninit_(rn);
	float *x = new float[n]; // 2 Parameter
	// Aufruf von PIKAIA
	cout << "#" << comment << endl;
	pikaia_(func,n,ctrl,x,&f,&status);
	// Auswerten
	cout 	<< "#status = " << status << endl
			<< "#x = " << x[0] << endl
			<< "#y = " << x[1] << endl;
}


int main(){  // aufrufendes Hauptprogramm
	maximize_function(P1,"Fitting P1");
	cout << endl;
	maximize_function(P2,"Fitting P2");
	cout << endl;


	// Initialize data
	ifstream infile("lichtkurve.dat");
	float xin, yin;
	while (infile >> xin >> yin) {
		X.push_back(xin);
		Y.push_back(yin);
		Sigma.push_back(5.0);
	}
	Npts = X.size();


	fit_data_lin_with_sins(0,false);
	cout << endl;
	fit_data_lin_with_sins(1,false);
	cout << endl;
	fit_data_lin_with_sins(3,false);
	cout << endl;
	fit_data_lin_with_sins(5,false);
	cout << endl;



}
