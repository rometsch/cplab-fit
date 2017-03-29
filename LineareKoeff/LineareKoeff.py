# -*- coding: utf-8 -*-

import numpy as np;
import matplotlib.pyplot as plt;

def polynomials(X,Y,Err,N,plot=False):
    
    ErrSq = Err*Err;        
    # Construct coefficient vector
    b = np.zeros(N);
    for i in range(N):
        f = lambda x: x**i;
        b[i] = np.sum(f(X)*Y/ErrSq);

    # Construct matrix
    M = np.zeros([N,N])
    for i in range(N):
        for j in range(N):
            f1 = lambda x: x**i;
            f2 = lambda x: x**j;
            M[i,j] = np.sum(f1(X)*f2(X)/ErrSq);
    
    # Solve the linear equation system
    coeff_vec = np.linalg.solve(M,b);
    
    func = np.polynomial.polynomial.Polynomial(coeff_vec);

    ChiSq = np.sum( (Y-func(X))**2/ErrSq )/(X.size-N);

    if plot:
        title = "Polynomial degree {}, ChiSq = {}".format(N-1, ChiSq)    
        plot_data_along_func(X,Y,Err,func,title);
    
    return ChiSq

def legendre(X,Y,Err,N,plot=False):
    ErrSq = Err*Err;    
    
    # Construct coefficient vector
    b = np.zeros(N);
    for i in range(N):
        # First get nth Legendre Polynomial
        mask = np.zeros(N);
        mask[i] = 1;
        f = np.polynomial.legendre.Legendre(mask);
        b[i] = np.sum(f(X)*Y/ErrSq);

    # Construct matrix
    M = np.zeros([N,N])
    for i in range(N):
        for j in range(N):
            mask1 = np.zeros(N);
            mask2 = np.zeros(N);
            mask1[i] = 1;
            mask2[j] = 1;
            f1 = np.polynomial.legendre.Legendre(mask1);
            f2 = np.polynomial.legendre.Legendre(mask2);
            M[i,j] = np.sum(f1(X)*f2(X)/ErrSq);
    
    # Solve the linear equation system
    coeff_vec = np.linalg.solve(M,b);
    
    func = np.polynomial.legendre.Legendre(coeff_vec);

    ChiSq = np.sum( (Y-func(X))**2/ErrSq )/(X.size-N);

    
    if plot:
        title = "Legendre Polynomial degree {}, ChiSq = {}".format(N-1, ChiSq)    
        plot_data_along_func(X,Y,Err,func,title);
    
    return ChiSq


def plot_data_along_func(X,Y,Err,func,title):
    x_small_spaced = np.linspace(X[0],X[-1],500)
    plt.errorbar(X,Y,yerr=Err,fmt="b.");
    plt.plot(x_small_spaced,func(x_small_spaced))
    plt.grid();
    plt.xlabel("x");
    plt.ylabel("y");
    plt.title(title);
    plt.show();
    
    
def main():
    X = np.array([-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9]);
    Y = np.array([81, 50, 35, 27, 26, 60, 106, 189, 318, 520]);
    Err = np.sqrt(Y);
    
    N_max = 8;  # Max order of polynomial
    
    degree = np.arange(N_max);
    poly_chisq = np.zeros(N_max);
    legendre_chisq = np.zeros(N_max);
    for n in range(0,N_max):
        poly_chisq[n] = polynomials(X,Y,Err,n+1, plot=True);
        legendre_chisq[n] = legendre(X,Y,Err,n+1);

    output = ""
    for n in degree:
        output += "\t&\t{}".format(n);
    output += "\\\\\n";
    for chisq in poly_chisq:
        output += "\t&\t{:.2f}".format(chisq);
        
    print(output)
    
    
    
if __name__=="__main__":
    main();