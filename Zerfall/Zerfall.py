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

Ns = data[:,1];
times = data[:,0]*delta_t;

# Calculate errors of counts as sqrt(N)
errsq = Ns;
# errsq = np.ones(Ns.size)

def print_results(params, method):
    print("==============================");
    print(" Fitting data with {}".format(method));
    print("==============================");
    print(" lambda 1 = {:.3e}".format(params[0]));
    print(" lambda 2 = {:.3e}".format(params[1]));
    print("     N1   = {:.3e}".format(params[2]));
    print("     N2   = {:.3e}".format(params[3]));
    print("     N0   = {:.3e}".format(params[4]));
    print("==============================");
    print("  Half-life 1 = {:.3e} sec".format(np.log(2)/params[0]));
    print("  Half-life 2 = {:.3e} sec".format(np.log(2)/params[1]));
    print("==============================");

def func(params,t):
    l1 = params[0];
    l2 = params[1];
    N1 = params[2];
    N2 = params[3];
    N0 = params[4];
    return N1*np.exp(-l1*t) + N2*np.exp(-l2*t) + N0;

def chisq(params):
    """ Defines the ChiSq for the model and the given data """
    f_vals = func(params,times);
    return np.sum( (Ns - f_vals)**2/errsq )

def chisq_grad(params):
    """ Defines the ChiSq for the model and the given data """
    l1 = params[0];
    l2 = params[1];
    N1 = params[2];
    N2 = params[3];
    f_vals = func(params,times);
    df_dl1 = (-times*N1)*np.exp(-l1*times);
    df_dl2 = (-times*N2)*np.exp(-l2*times);
    df_dN1 = np.exp(-l1*times);
    df_dN2 = np.exp(-l2*times);
    df_dN0 = 1;
    gradient = -2*np.array(
        [   np.sum( df_dl1*(Ns-f_vals)/errsq ),
            np.sum( df_dl2*(Ns-f_vals)/errsq ),
            np.sum( df_dN1*(Ns-f_vals)/errsq ),
            np.sum( df_dN2*(Ns-f_vals)/errsq ),
            np.sum( df_dN0*(Ns-f_vals)/errsq )
        ]);
    return gradient;
    # Division by gamma comes from the fact, that err = sqrt(gamma)


def plot_data_fit(params,times):
    """ Plot the data along with the fitted function """
    plt.plot(times, Ns, "b.");
    plt.plot(times, func(params,times), "g");
    plt.show()

# initial_guess = np.array([0.03,3e-3,655,104,5]);
initial_guess = np.array([1e-1,1e-2,1000,100,10]);

res_simplex = optimize.minimize(chisq,initial_guess,method='Nelder-Mead');
#res_cg = optimize.minimize(chisq,initial_guess,method='CG',jac=chisq_grad);
res_cg = optimize.fmin_cg(chisq, initial_guess, fprime=chisq_grad);
# print(params_cg)

print(res_simplex)
print_results(res_simplex.x,"Simplex");
print(res_cg)
print_results(res_cg,"CG");

# plot_data_fit(res_simplex.x,times)
