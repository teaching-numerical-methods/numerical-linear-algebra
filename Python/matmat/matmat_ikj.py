# matmat_ikj.py
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

def matmat_ikj(A,X):
    """
    Usage: B = matmat_ikj(A,X)

    Function to perform matrix-matrix multiplication with "ikj" loop ordering.

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
        raise ValueError("matmat_ikj error: A and X are incompatible")

    # initialize output
    B = np.zeros([m,p], dtype=float)

    # perform product
    for i in range(m):
        for k in range(n):
            B[i,:] += A[i,k] * X[k,:]

    return B
