#!/usr/bin/env python3
#
# Script to test Cholesky factorizations
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

# imports
import time
import numpy as np
from numpy.linalg import norm
from cholesky_ip import *
from cholesky_op import *
from fwdsub_row import *
from bwdsub_row import *

# set testing values
nvals = [500, 700, 900, 1100]

# run tests
for n in nvals:

    # allocate the matrices & vectors of this size
    A = np.zeros([n,n], dtype=float)
    R = np.zeros([n,n], dtype=float)
    Rt = np.zeros([n,n], dtype=float)
    xtrue = np.zeros(n, dtype=float)
    x = np.zeros(n, dtype=float)
    y = np.zeros(n, dtype=float)
    b = np.zeros(n, dtype=float)

    # fill A and xtrue with values
    for i in range(n):
        for j in range(n):
            A[i,j] = 1.0/(1.0 + 5.0*abs(i - j))
    for i in range(n):
        xtrue[i] = 1.0*(1 - i)/n

    # compute b from A and xtrue
    b = A@xtrue

    # start first test
    print("\nTesting Cholesky (IP) factorizations: n = ",n)

    # compute Cholesky decomposition, placing result in R
    R = A.copy()
    stime = time.time()
    if (cholesky_ip(R) != 0):
        print("cholesky_ip failed")
    runtime = time.time()-stime
    print("   cholesky time = ", runtime)

    # fill Rt as the transpose of R
    Rt = np.transpose(R)

    # solve linear system
    stime = time.time()
    y = fwdsub_row(Rt, b)
    x = bwdsub_row(R, y)
    runtime = time.time()-stime
    print("   solve time = ", runtime)

    # check error
    err_norm = norm(x-xtrue)
    print("   solution error = ", err_norm)


    # reset A, xtrue and b to original values (since Python is sometimes
    # 'call-by-value' and other times 'call-by-reference')
    for i in range(n):
        for j in range(n):
            A[i,j] = 1.0/(1.0 + 5.0*abs(i - j))
    for i in range(n):
        xtrue[i] = 1.0*(1 - i)/n
    b = A@xtrue

    # start second test
    print("\nTesting Cholesky (OP) factorizations: n = ",n)

    # compute Cholesky decomposition, storing result in R
    R = A.copy()
    stime = time.time()
    if (cholesky_op(R) != 0):
        print("cholesky_op failed")
    runtime = time.time()-stime
    print("   cholesky time = ", runtime)

    # fill Rt as the transpose of R
    Rt = np.transpose(R)

    # solve linear system
    stime = time.time()
    y = fwdsub_row(Rt, b)
    x = bwdsub_row(R, y)
    runtime = time.time()-stime
    print("   solve time = ", runtime)

    # check error
    err_norm = norm(x-xtrue)
    print("   solution error = ", err_norm)

