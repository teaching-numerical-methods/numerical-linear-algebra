# fwdsub_row.py
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

def fwdsub_row(L, b):
    """
    Usage: y = fwdsub_row(L, b) 

    Function to perform row-oriented forwards substitution
 
    Inputs: 
       L is a lower-triangular matrix (n x n numpy matrix)
       b is a vector (n numpy array)

    Outputs:
       y is a vector (n numpy array)
    """

    # imports
    import numpy as np

    # get problem dimensions
    m, n = np.shape(L)

    # check that L and b are compatible
    if (m != n):
        raise ValueError("fwdsub_row error: L is not square")
    if (n != b.size):
        raise ValueError("fwdsub_row error: L and b are incompatible")

    # check that L is nonsingular
    for i in range(n):
        if (L[i,i] == 0.0):
            raise ValueError("fwdsub_row error: L is singular")

    # copy b into solution vector
    y = b

    # loop over matrix rows, performing solve
    for i in range(n):
        for j in range(i):   # update this rhs
            y[i] -= L[i,j]*y[j]
        y[i] /= L[i,i]       # solve row

    return y
