#!/usr/bin/env python3
#
# Script to test matrix-matrix products
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

# imports
import time
import numpy as np
from matmat_ijk import *
from matmat_ikj import *
from matmat_jik import *
from matmat_jki import *
from matmat_kij import *
from matmat_kji import *

# set testing values
mvals = [50, 100, 200, 400]
nvals = [100, 200, 400, 800]
pvals = [75, 150, 300, 600]
nsizes = 4

# run tests
for k in range(nsizes):

    # set the problem size
    m = mvals[k]
    n = nvals[k]
    p = pvals[k]

    # display current problem size
    print("\nTesting matrix-matrix products: m = ",m," n = ",n," p = ",p)

    # allocate the matrix & vectors of this size
    A = np.zeros([m,n], dtype=float)
    X = np.zeros([n,p], dtype=float)

    # fill A and X with values
    for i in range(m):
        for j in range(n):
            A[i,j] = 1.0*(1 + i - j)/(n+m)
    for i in range(n):
        for j in range(p):
            X[i,j] = 1.0*(1 - i + j)/(n+p)

    # perform product 1
    stime = time.time()
    B = matmat_ijk(A, X)
    runtime = time.time()-stime
    B_err = np.max(np.max(np.abs(B-A@X)))
    print("   matmat_ijk:  time = ", runtime,", error = ", B_err)

    # perform product 2
    stime = time.time()
    B = matmat_ikj(A, X)
    runtime = time.time()-stime
    B_err = np.max(np.max(np.abs(B-A@X)))
    print("   matmat_ikj:  time = ", runtime,", error = ", B_err)

    # perform product 3
    stime = time.time()
    B = matmat_jik(A, X)
    runtime = time.time()-stime
    B_err = np.max(np.max(np.abs(B-A@X)))
    print("   matmat_jik:  time = ", runtime,", error = ", B_err)

    # perform product 4
    stime = time.time()
    B = matmat_jki(A, X)
    runtime = time.time()-stime
    B_err = np.max(np.max(np.abs(B-A@X)))
    print("   matmat_jki:  time = ", runtime,", error = ", B_err)

    # perform product 5
    stime = time.time()
    B = matmat_kij(A, X)
    runtime = time.time()-stime
    B_err = np.max(np.max(np.abs(B-A@X)))
    print("   matmat_kij:  time = ", runtime,", error = ", B_err)

    # perform product 6
    stime = time.time()
    B = matmat_kji(A, X)
    runtime = time.time()-stime
    B_err = np.max(np.max(np.abs(B-A@X)))
    print("   matmat_kji:  time = ", runtime,", error = ", B_err)

    # call 'dot' for product 
    stime = time.time()
    B = np.dot(A, X)
    runtime = time.time()-stime
    B_err = np.max(np.max(np.abs(B-A@X)))
    print("   dot:         time = ", runtime,", error = ", B_err)

    # call 'matmul' for product 
    stime = time.time()
    B = np.matmul(A, X)
    runtime = time.time()-stime
    B_err = np.max(np.max(np.abs(B-A@X)))
    print("   matmul:      time = ", runtime,", error = ", B_err)

    # call '@' operator for product 
    stime = time.time()
    B = A @ X
    runtime = time.time()-stime
    B_err = np.max(np.max(np.abs(B-A@X)))
    print("   @ operator:  time = ", runtime,", error = ", B_err)
