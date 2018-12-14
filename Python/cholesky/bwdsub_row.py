# bwdsub_row.py
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

def bwdsub_row(U, y):
    """
    Usage: x = bwdsub_row(U, y) 

    Function to perform row-oriented backwards substitution
 
    Inputs: 
       U is an upper-triangular matrix (n x n numpy matrix)
       y is a vector (n numpy array)

    Outputs:
       x is a vector (n numpy array)
    """

    # imports
    import numpy as np

    # get problem dimensions
    m, n = np.shape(U)

    # check that U and y are compatible
    if (m != n):
        raise ValueError("bwdsub_row error: U is not square")
    if (n != y.size):
        raise ValueError("bwdsub_row error: U and y are incompatible")

    # check that U is nonsingular
    for i in range(n):
        if (U[i,i] == 0.0):
            raise ValueError("bwdsub_row error: U is singular")

    # copy y into solution vector
    x = y

    # loop over matrix rows, performing solve
    for i in range(n-1,-1,-1):
        for j in range(i+1, n):  # update this rhs
            x[i] -= U[i,j]*x[j]
        x[i] /= U[i,i];          # solve row

    return x
