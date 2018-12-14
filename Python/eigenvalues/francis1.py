# francis1.py
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

def francis1(Ain,maxit,tol,stype,diags):
    """
    Usage: A,Q,its = francis1(Ain,maxit,tol,stype,diags)

    Function to perform Francis' algorithm of degree one to compute the
    eigen-decomposition of a given matrix A, through iteratively computing the
    similarity transformation 
           Anew = Q'*Aold*Q
    where Aold is an upper-Hessenberg matrix, Q is a unitary transformation matrix,
    and Anew is an upper-Hessenberg matrix with decreasing subdiagonal.  At
    convergence, Anew should be [essentially] upper-triangular.  Specifically,
    it will satisfy
          |Anew(j,j-1)| <= tol*(|Anew(j,j)|+|Anew(j-1,j-1)|),  j=2:n
   
    Note: since this routine returns the full transformation matrix Q (as
    opposed to just the triangular matrix result, it is less efficient than it
    *could* be if only the eigenvalues were desired.
   
    Input:    A - upper-Hessenberg square matrix (Aold)
          maxit - maximum allowed number of iterations
            tol - relative eigenvalue tolerance
          stype - desired shift approach:
                        0 => Rayleigh quotient shift
                     else => Wilkinson shift
          diags - flag to turn on/off diagnostic messages:
                        0 => off
                     else => on
    Outputs:  Q - unitary matrix
              A - upper-triangular matrix (Anew)
            its - number of iterations taken
    """

    # imports
    import numpy as np
    from numpy.linalg import norm
    from francis_step import francis_step

    # copy input matrix into new working matrix
    A = Ain.copy()
    
    # ensure that A is square
    m,n = np.shape(A)
    if (m != n):
        raise ValueError("francis1 error: matrix must be square")

    # ensure that A is essentially upper-Hessenberg
    tmp = A - np.triu(A,-1)  # portion of A that should be zero
    if (norm(tmp,np.inf)/norm(A,np.inf) > 100*np.finfo(float).eps):
        raise ValueError("francis1: matrix must be upper-Hessenberg")

    # initialize results, counter toward completion
    Q = np.eye(n)
    m = n-1

    # perform iteration
    for its in range(1,maxit+1):

        # determine last proper Hessenberg row
        for k in range(m,0,-1):
            if (np.abs(A[k,k-1]) > tol*(np.abs(A[k,k])+np.abs(A[k-1,k-1]))):
                break
        m = k
            
        # check for completion
        if ((m == 1) and (np.abs(A[1,0]) <= tol*(np.abs(A[1,1])+np.abs(A[0,0])))):
            break
   
        # set shift
        if (stype == 0):       # Rayleigh quotient shift
            rho = A[m,m]
        else:                  # Wilkinson shift

            # set terms from quadratic equation, a*lambda^2 + b*lambda + c = 0,
            # for eigenvalues of submatrix A[m-1:m,m-1:m]
            a = 1.0
            b = -A[m-1,m-1] - A[m,m]
            c = A[m-1,m-1]*A[m,m] - A[m-1,m]*A[m,m-1]

            # compute eigenvalues of submatrix
            disc = b**2-4.0*a*c
            if (disc < 0):
                sqdisc = np.sqrt(complex(disc))
            else:
                sqdisc = np.sqrt(disc)
            s1 = (-b + sqdisc)/(2.0*a)
            s2 = (-b - sqdisc)/(2.0*a)

            # set shift as closest eigenvalue to A[m,m] component
            if (np.abs(s1-A[m,m]) < np.abs(s2-A[m,m])):
                rho = s1
            else:
                rho = s2

        # output some diagnostic information to the screen
        if (diags):
            print("   iter ",its,":  submatrix 0 :",m,",  shift = ",rho)
   
        # perform one iteration of algorithm on remaining submatrix
        Asub, Qt = francis_step(A[:m+1,:m+1],rho)
        if isinstance(Asub[0,0], complex):
            A = A.astype(complex)
            Q = Q.astype(complex)
        A[:m+1,:m+1] = Asub
        
        # update transformation matrix
        Q[:,:m+1] = Q[:,:m+1] @ Qt

    return [A, Q, its]
