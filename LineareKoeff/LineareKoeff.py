# -*- coding: utf-8 -*-

import numpy as np;
import matplotlib.pyplot as plt;

def polynomials(X,Y,Err,N):    
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

    title = "Polynomial degree {}, ChiSq = {}".format(N-1, ChiSq)    
    
    plot_data_along_func(X,Y,Err,func,title);
    

def plot_data_along_func(X,Y,Err,func,title):
    x_small_spaced = np.linspace(X[0],X[-1],500)
    plt.plot(X,Y,"b.",x_small_spaced,func(x_small_spaced));
    plt.grid();
    plt.xlabel("x");
    plt.ylabel("y");
    plt.title(title);
    plt.show();
    
    
def main():
    X = np.array([-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9]);
    Y = np.array([81, 50, 35, 27, 26, 60, 106, 189, 318, 520]);
    Err = np.sqrt(Y);
    
    N = 5;  # Order of polynomial
    
    polynomials(X,Y,Err,N+1);

if __name__=="__main__":
    main();