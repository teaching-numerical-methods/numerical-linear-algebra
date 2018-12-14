# matvec_col.py
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

def matvec_col(A, x):
    """
    Usage: b = matvec_col(A, x) 

    Function to perform column-oriented matrix multiplication
 
    Inputs: 
       A is a matrix (m x n numpy matrix)
       x is a vector (n numpy array)

    Outputs:
       b is a vector (m numpy array)
    """

    # imports
    import numpy as np

    # get problem dimensions
    m, n = np.shape(A)

    # check that A and x are compatible
    if (n != x.size):
        raise ValueError("matvec_col error: A and x are incompatible")
        
    # initialize output
    b = np.zeros(m, dtype=float)
    
    # perform product
    for j in range(n):
        b += A[:,j]*x[j]

    return b


      
