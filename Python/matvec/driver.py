#!/usr/bin/env python3
#
# Script to test matrix-vector products
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

# imports
import time
import numpy as np
from numpy.linalg import norm
from matvec_row import *
from matvec_col import *

# set testing values
mvals = [1000, 2000, 4000, 8000]
nvals = [2000, 4000, 8000, 16000]
nsizes = 4

# run tests
for k in range(nsizes):

    # set the problem size
    m = mvals[k]
    n = nvals[k]

    # display current problem size
    print("\nTesting matrix-vector products with a ", m, " x ", n, " matrix:")

    # allocate the matrix & vectors of this size
    A = np.zeros([m,n], dtype=float)
    x = np.zeros(n, dtype=float)

    # fill A and x with values
    for i in range(m):
        for j in range(n):
            A[i,j] = 1.0*(1 + i - j)/(n+m)
    for j in range(n):
        x[j] = 1.0*j/n

    # perform row-based product
    stime = time.time()
    b = matvec_row(A, x)
    runtime = time.time()-stime
    b_err = norm(b-A@x)
    print("   matvec_row:  time = ", runtime, ", error = ", b_err)

    # perform column-based product
    stime = time.time()
    b = matvec_col(A, x)
    runtime = time.time()-stime
    b_err = norm(b-A@x)
    print("   matvec_col:  time = ", runtime, ", error = ", b_err)

    # call 'dot' for product
    stime = time.time()
    b = np.dot(A,x)
    runtime = time.time()-stime
    b_err = norm(b-A@x)
    print("   dot:         time = ", runtime, ", error = ", b_err)

    # call 'matmul' for product
    stime = time.time()
    b = np.matmul(A,x)
    runtime = time.time()-stime
    b_err = norm(b-A@x)
    print("   matmul:      time = ", runtime, ", error = ", b_err)

    # call '@' operator for product
    stime = time.time()
    b = A @ x
    runtime = time.time()-stime
    b_err = norm(b-A@x)
    print("   @ operator:  time = ", runtime, ", error = ", b_err)
