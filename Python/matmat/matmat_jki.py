# matmat_jki.py
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

def matmat_jki(A,X):
    """
    Usage: B = matmat_jki(A,X)

    Function to perform matrix-matrix multiplication with "jki" loop ordering.

    Inputs: A is a matrix (m x n numpy matrix)
            X is a matrix (n x p numpy matrix)
    Output: B is a matrix (m x p numpy matrix)
    """

    # imports
    import numpy as np

    # get problem dimensions
    m, n = np.shape(A)
    n2, p = np.shape(X)

    # check that A and X are compatible
    if (n != n2):
        raise ValueError("matmat_jki error: A and X are incompatible")

    # initialize output
    B = np.zeros([m,p], dtype=float)

    # perform product
    for j in range(p):
        for k in range(n):
            B[:,j] += A[:,k] * X[k,j]

    return B
