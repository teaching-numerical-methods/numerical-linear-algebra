#!/usr/bin/env python3
#
# Script to test out approaches for converting a matrix to upper Hessenberg form.
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

# imports
import numpy as np
from numpy.linalg import norm
from upper_hess import upper_hess

# adjust output precision
np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)

# non-symmetric matrix
n = 10
A = np.random.rand(n,n)
print("first matrix (non-symmetric):")
print(A)
input("Press Enter to continue...")

# convert to upper-Hessenberg form
H, Q = upper_hess(A)
print("upper-Hessenberg H = Q^T A Q:")
print(H)
print("checks:")
print("  norm of lower portion of H = ", norm(H - np.triu(H,-1)))
print("  ||Q^H Q - I|| = ", norm(Q.T.conj() @ Q - np.eye(n)))
print("  ||Q Q^H - I|| = ", norm(Q @ Q.T.conj() - np.eye(n)))
print("  ||Q^H A Q - H|| = ", norm(Q.T.conj() @ A @ Q -H), "\n")


# symmetric matrix
n = 10
A = np.random.rand(n,n)
A = A + A.T
print("second matrix (symmetric):")
print(A)
input("Press Enter to continue...")

# convert to upper-Hessenberg form
H, Q = upper_hess(A)
print("upper-Hessenberg H = Q^T A Q:")
print(H)
print("checks:");
print("  norm of lower portion of H = ", norm(H - np.triu(H,-1)))
print("  ||Q^T Q - I|| = ", norm(Q.T @ Q - np.eye(n)))
print("  ||Q Q^T - I|| = ", norm(Q @ Q.T - np.eye(n)))
print("  ||Q^T A Q - H|| = ", norm(Q.T @ A @ Q - H), "\n")
