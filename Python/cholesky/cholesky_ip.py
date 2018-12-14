# cholesky_ip.py
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

def cholesky_ip(A):
    """
    Usage: iret = cholesky_ip(A) 

    Function to perform an inner-product-oriented Cholesky factorization
 
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
    for i in range(n):             # loop over rows of result
        A[i,i] -= sum(A[:i,i]**2)  # update diagonal
        if (A[i,i] <= 0):          # check positive definite
            print("cholesky_ip error: A is not positive definite")
            return 1
        A[i,i] = np.sqrt(A[i,i])   # set diagonal entry for row
        for j in range(i+1,n):     # loop over remainder of row
            A[i,j] -= sum(A[:i,i]*A[:i,j])    # update row entry
            A[i,j] /= A[i,i]        # set row entry

    return 0

