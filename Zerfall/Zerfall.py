#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 13:30:22 2017

@author: thomas
"""

import numpy as np;
import scipy.optimize as optimize;

data = np.loadtxt("Agdecay.dat");

delta_t = 5;

gammas = data[:,1];
times = data[:,0]*delta_t;

# Calculate errors of counts as sqrt(N)
errsq = gammas;

def print_results(params, method):
    print("==============================");
    print(" Fitting data with {}".format(method));
    print("==============================");
    print(" lambda 1 = {:.3f}".format(params[0]));
    print(" lambda 2 = {:.3f}".format(params[1]));
    print("     N1   = {:.3f}".format(params[2]));
    print("     N2   = {:.3f}".format(params[3]));
    print("  gamma 0 = {:.3f}".format(params[4]));
    print("==============================");
    print("==============================");

def chisq(params):
    """ Defines the ChiSq for the model and the given data """
    l1 = params[0];
    l2 = params[1];
    N1 = params[2];
    N2 = params[3];
    gamma0 = params[4];
    f_vals = l1*np.exp(-l1*times) + l2*np.exp(-l2*times) + gamma0;
    return np.sum( (gammas - f_vals)**2/errsq )
    # Division by gamma comes from the fact, that err = sqrt(gamma)

def chisq_grad(params):
    """ Defines the ChiSq for the model and the given data """
    l1 = params[0];
    l2 = params[1];
    N1 = params[2];
    N2 = params[3];
    gamma0 = params[4];
    f_vals = l1*np.exp(-l1*times) + l2*np.exp(-l2*times) + gamma0;
    df_dl1 = (1-l1**2)*np.exp(-l1*times)*N1;
    df_dl2 = (1-l1**2)*np.exp(-l2*times)*N2;
    df_dN1 = l1*np.exp(-l1*times);
    df_dN2 = l2*np.exp(-l2*times);
    df_dgamma0 = 1;
    gradient = -2*np.array(
        [   np.sum( df_dl1*f_vals/errsq ),
            np.sum( df_dl2*f_vals/errsq ),
            np.sum( df_dN1*f_vals/errsq ),
            np.sum( df_dN2*f_vals/errsq ),
            np.sum( df_dgamma0*f_vals/errsq )
        ]);
    return gradient;
    # Division by gamma comes from the fact, that err = sqrt(gamma)


initial_guess = np.array([1e-5,1e-6,10000,1000,10]);

res_cg = optimize.minimize(chisq,initial_guess,jac=chisq_grad,method='CG');
# print(params_cg)
res_simplex = optimize.minimize(chisq,initial_guess,method='Nelder-Mead');

print_results(res_simplex.x,"Simplex")
print(res_cg)
