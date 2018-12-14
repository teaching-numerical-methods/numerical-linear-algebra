# upper_hess.py
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

def upper_hess(A):
    """
    Usage: H,Q = upper_hess(A) 

    Function to convert a matrix to upper-Hessenberg form via the unitary
    similarity transformation
           H = Q'*A*Q
    where A is a general square matrix, Q is a unitary transformation matrix,
    and H is the upper Hessenberg result.
   
    Note: since this routine returns the full transformation matrix Q (and not
    just its parts), it is less efficient than it *could* be.
   
    Input:    A - square matrix
    Outputs:  Q - unitary matrix
              H - upper Hessenberg matrix
    """

    # imports
    import numpy as np

    # ensure that A is square
    m,n = np.shape(A)
    if (m != n):
        raise ValueError("upper_hess error: matrix must be square")

    # initialize results
    Q = np.eye(n)
    H = A.copy()

    # iterate over columns
    for j in range(n-2):

        # construct reflector to zero out column below first subdiagonal
        b = H[j+1:,j]                            # column to transform
        tau = np.linalg.norm(b)*np.sign(b[0]+np.finfo(float).eps)   # norm of column
        u = b/(b[0]+tau)                         # construct u = (b-y)/(b1+tau), where y = [-tau,0,...,0]
        u[0] = 1.0
        gamma = (tau+b[0])/tau
        Qj = np.eye(n-j-1) - gamma*np.outer(u,u)   # create Householder matrix

        # multiply remaining submatrices of H and Q by Qj
        H[j+1:,:] = Qj @ H[j+1:,:]
        H[:,j+1:] = H[:,j+1:] @ Qj
        Q[:,j+1:] = Q[:,j+1:] @ Qj
   
    return [H, Q]
