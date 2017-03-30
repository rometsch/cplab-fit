#!/usr/local/bin/python3

import numpy as np;

data = np.loadtxt("Agdecay.dat");

delta_t = 5;

Ns = data[:,1];
times = data[:,0]*delta_t;
errsq = Ns;

l1 = 0.02844;
l2 = 0.003326;

SN1 = np.sum( Ns*np.exp(-l1*times)/errsq );
SN2 = np.sum( Ns*np.exp(-l2*times)/errsq );
S11 = np.sum( np.exp(-2*l1*times)/errsq );
S22 = np.sum( np.exp(-2*l2*times)/errsq );
S12 = np.sum( np.exp(-(l1+l2)*times)/errsq );
S1 = np.sum( np.exp(-l1*times)/errsq );
S2 = np.sum( np.exp(-l2*times)/errsq );
SN = np.sum( Ns/errsq );
S = np.sum( 1.0/errsq );

b = np.array([SN1,SN2,SN]);
M = np.array([  [S11,S12,S1],
                [S12,S22,S2],
                [S1,S2,S]   ]);

coeff_vec = np.linalg.solve(M,b);

N1 = coeff_vec[0];
N2 = coeff_vec[1];
N0 = coeff_vec[2];

print("==============================");
print("Fit data with fixed lambdas")
print("==============================");
print("    N1    =    {}".format(N1));
print("    N2    =    {}".format(N2));
print("    N0    =    {}".format(N0));
print("==============================");
print(" given lambdas")
print("  lambda1 =    {:.3e}".format(l1));
print("  lambda2 =    {:.3e}".format(l2));
print("==============================");
print("==============================");
