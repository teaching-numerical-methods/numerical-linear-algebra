# matvec_row.py
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

def matvec_row(A, x):
    """
    Usage: b = matvec_row(A, x) 

    Function to perform row-oriented matrix multiplication
 
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
        raise ValueError("matvec_row error: A and x are incompatible")
        
    # initialize output
    b = np.zeros(m, dtype=float)
    
    # perform product
    for i in range(m):
        b[i] += sum(A[i,:]*x)

    return b


      
