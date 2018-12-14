# cholesky_op.py
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

def cholesky_op(A):
    """
    Usage: iret = cholesky_op(A) 

    Function to perform an outer-product-oriented Cholesky factorization
 
    Inputs: 
       A is a symmetric matrix (n x n numpy matrix)

    Outputs:
        A is updated in-place to store the Cholesky-factored matrix, A = R'*R, 
         with R stored in the upper-triangular portion of A
       iret is a success/failure flag (0=success, 1=failure)
    """

    # imports
    import numpy as np

    # get problem dimensions
    m, n = np.shape(A)

    # check that A is square
    if (m != n):
        print("cholesky_ip error: A is not square")
        return 1

    # perform factorization in-place
    for i in range(n):            # loop over rows of result
        if (A[i,i] <= 0):         # check positive definite
            print("cholesky_ip error: A is not positive definite")
            return 1
        A[i,i] = np.sqrt(A[i,i])  # set diagonal entry for row
        A[i,i+1:] /= A[i,i]       # update row
        for k in range(i+1,n):    # update remainder of matrix
            A[k,i+1:] -= A[i,k]*A[i,i+1:]

    return 0
