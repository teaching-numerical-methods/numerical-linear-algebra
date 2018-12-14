# francis_step.py
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019


#----------------------------------------
# utility routines

def rot(x1,x2):
    """
    Usage: c,s = rot(x1,x2)
  
    Function to compute cosine and sine components of a rotation matrix,
    given the two arguments:
  
    Input:   x1 - value to use in elimination
             x2 - value to be eliminated
    Outputs:  c - cosine term
              s - sine term
    """

    # imports
    import numpy as np

    beta = max(np.abs(x1),np.abs(x2))
    if (beta == 0):  # compute c, s values
       c = 1.0
       s = 0.0
    else:
       x1 = x1/beta
       x2 = x2/beta
       v = np.sqrt(x1**2+x2**2)
       c = x1/v
       s = x2/v
    return [c,s]


def apply_rot_left(Mi,Mj,c,s):
    """
    Usage: [QMi,QMj] = apply_rot_left(Mi,Mj,c,s)
   
    Applies a rotation matrix multiply on the left,
              QM = Q'*M
    where Q is a rotation matrix with block [c -s; s c]
  
    Input:    Mi - first row from matrix to be multiplied
              Mj - second row from matrix to be multiplied
               c - cosine term in rotator
               s - sine term in rotator
    Outputs: QMi - first row of resulting multiply
             QMj - second row of resulting multiply
    """
    QMi = c*Mi + s*Mj
    QMj = c*Mj - s*Mi

    return [QMi, QMj]

  
def apply_rot_right(Mi,Mj,c,s):
    """
    Usage: [MQi,MQj] = apply_rot_right(Mi,Mj,c,s)
  
    Applies a rotation matrix multiply on the right
             MQ = M*Q
    where Q is a rotation matrix with block [c -s; s c]
  
    Input:    Mi - first column from matrix to be multiplied
              Mj - second column from matrix to be multiplied
               c - cosine term in rotator
               s - sine term in rotator
    Outputs: MQi - first column of resulting multiply
             MQj - second column of resulting multiply
    """
    MQi = c*Mi + s*Mj
    MQj = c*Mj - s*Mi
    return [MQi, MQj]



#----------------------------------------
# primary routine

def francis_step(Ain,rho):
    """
    Usage: Anew, Q = francis_step(Aold,rho)
   
    Function to perform one iteration of Francis' algorithm of degree one.
    Computes the unitary similarity transformation 
           Anew = Q'*Aold*Q
    where Aold is an upper-Hessenberg matrix, Q is a unitary transformation matrix,
    and Anew is an upper-Hessenberg result, through introducing a bulge in Aold (via
    a rotator), and chasing it back out (via a sequence of rotators).
   
    Note: since this routine returns the full transformation matrix Q (and not
    just its parts), it is less efficient than it *could* be.
   
    Input:    A - upper-Hessenberg square matrix (Aold)
              rho - desired shift
    Outputs:  Q - unitary matrix
              A - upper-Hessenberg matrix result (Anew)
    """

    # imports
    import numpy as np

    # get matrix size
    n = np.size(Ain,1)   # assumed square

    # shift the input matrix by rho
    A = Ain - rho*np.eye(n)

    # introduce bulge
    c,s = rot(A[0,0], A[1,0])     # compute rotation matrix components

    # apply rotation to A on left: A = Qij'*A
    v1, v2 = apply_rot_left(A[0,:],A[1,:],c,s)
    A[0,:] = v1
    A[1,:] = v2

    # apply rotation to A on right: A = A*Qij
    v1, v2 = apply_rot_right(A[:,0],A[:,1],c,s)
    A[:,0] = v1
    A[:,1] = v2

    # create Q
    if isinstance(rho, complex):
        Q = np.eye(n, dtype=complex)
    else:
        Q = np.eye(n)
    Q[0,0] = c
    Q[0,1] = -s
    Q[1,0] = s
    Q[1,1] = c

  
    # chase the bulge
    for j in range(n-2):

        # construct rotator to zero non-Hessenberg part of column
        c, s = rot(A[j+1,j], A[j+2,j])

        # apply rotation to A on left: A = Qij'*A
        v1, v2 = apply_rot_left(A[j+1,:],A[j+2,:],c,s)
        A[j+1,:] = v1
        A[j+2,:] = v2

        # apply rotation to A on right: A = A*Qij
        v1, v2 = apply_rot_right(A[:,j+1],A[:,j+2],c,s)
        A[:,j+1] = v1
        A[:,j+2] = v2

        # update Q
        v1, v2 = apply_rot_right(Q[:,j+1],Q[:,j+2],c,s)
        Q[:,j+1] = v1
        Q[:,j+2] = v2

    # undo shift on result
    A = A + rho*np.eye(n)

    return [A, Q]
