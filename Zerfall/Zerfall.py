#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 13:30:22 2017

@author: thomas
"""

import numpy as np;
import scipy.optimize as optimize;
import matplotlib.pyplot as plt;

data = np.loadtxt("Agdecay.dat");

delta_t = 5;

gammas = data[:,1];
times = data[:,0]*delta_t;

# Calculate errors of counts as sqrt(N)
errsq = gammas;
#errsq = np.ones(gammas.size)

def print_results(params, method):
    print("==============================");
    print(" Fitting data with {}".format(method));
    print("==============================");
    print(" lambda 1 = {:.3e}".format(params[0]));
    print(" lambda 2 = {:.3e}".format(params[1]));
    print("     N1   = {:.3e}".format(params[2]));
    print("     N2   = {:.3e}".format(params[3]));
    print("  gamma 0 = {:.3e}".format(params[4]));
    print("==============================");
    print("  Half-life 1 = {:.3e} sec".format(np.log(2)/params[0]));
    print("  Half-life 2 = {:.3e} sec".format(np.log(2)/params[1]));
    print("==============================");

def func(params,t):
    l1 = params[0];
    l2 = params[1];
    N1 = params[2];
    N2 = params[3];
    gamma0 = params[4];    
    return l1*N1*np.exp(-l1*t) + l2*N2*np.exp(-l2*t) + gamma0;

def chisq(params):
    """ Defines the ChiSq for the model and the given data """
    f_vals = func(params,times);
    return np.sum( (gammas - f_vals)**2/errsq )
    # Division by gamma comes from the fact, that err = sqrt(gamma)

def chisq_grad(params):
    """ Defines the ChiSq for the model and the given data """
    l1 = params[0];
    l2 = params[1];
    N1 = params[2];
    N2 = params[3];
    f_vals = func(params,times);
    df_dl1 = (1-l1**2)*N1*np.exp(-l1*times);
    df_dl2 = (1-l2**2)*N2*np.exp(-l2*times);
    df_dN1 = l1*np.exp(-l1*times);
    df_dN2 = l2*np.exp(-l2*times);
    df_dgamma0 = 1;
    gradient = -2*np.array(
        [   np.sum( df_dl1*(gammas-f_vals)/errsq ),
            np.sum( df_dl2*(gammas-f_vals)/errsq ),
            np.sum( df_dN1*(gammas-f_vals)/errsq ),
            np.sum( df_dN2*(gammas-f_vals)/errsq ),
            np.sum( df_dgamma0*(gammas-f_vals)/errsq )
        ]);
    return gradient;
    # Division by gamma comes from the fact, that err = sqrt(gamma)


def plot_data_fit(params,times):
    """ Plot the data along with the fitted function """
    plt.plot(times, gammas, "b.");
    plt.plot(times, func(params,times), "g");

initial_guess = np.array([0.03,3e-3,28000,34000,5]);

res_simplex = optimize.minimize(chisq,initial_guess,method='Nelder-Mead');
res_cg = optimize.minimize(chisq,initial_guess,jac=chisq_grad,method='CG');
# print(params_cg)

print(res_simplex)
print_results(res_simplex.x,"Simplex");
print(res_cg)
#print_results(res_cg.x,"CG");

plot_data_fit(res_simplex.x,times)

