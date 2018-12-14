#!/usr/bin/env python3
#
# Script to show what happens with successive Francis iterations on a matrix
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

# imports
import numpy as np
from numpy.linalg import norm
from upper_hess import upper_hess
from francis1 import francis1
from francis_step import francis_step

# adjust output precision
np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)

# symmetric matrix
n = 10
A = np.random.rand(n,n)
A = A + A.T
print("first matrix A (symmetric):")
print(A)
input("Press Enter to continue...")

# convert to upper-Hessenberg form
H, Q = upper_hess(A)
print("upper-Hessenberg H = Q^T A Q:")
print(H)
input("Press Enter to continue...")
print("checks:")
print("  norm of lower portion of H = ", norm(H - np.triu(H,-1)))
print("  ||Q^T Q - I|| = ", norm(Q.T @ Q - np.eye(n)))
print("  ||Q Q^T - I|| = ", norm(Q @ Q.T - np.eye(n)))
print("  ||Q^T A Q - H|| = ", norm(Q.T @ A @ Q - H))
print("  subdiagonal of H:")
print(np.diag(H,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(H,-1)),"\n")
input("Press Enter to continue...")


# perform a few Francis steps (Rayleigh quotient shift) to show what happens
Hold = H.copy()
rho = Hold[n-1,n-1]
Hnew, U = francis_step(Hold,rho)
print("first Francis iteration of degree 1 (rho = ",rho,"), H = U^T H U:")
print(Hnew)
input("Press Enter to continue...")
print("checks:")
print("  ||U^T U - I|| = ", norm(U.T @ U - np.eye(n)))
print("  ||U U^T - I|| = ", norm(U @ U.T - np.eye(n)))
print("  ||U^T H U - Hnew|| = ", norm(U.T @ Hold @ U - Hnew))
print("  subdiagonal of H:")
print(np.diag(Hnew,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(Hnew,-1)),"\n")
input("Press Enter to continue...")

Hold = Hnew.copy()
rho = Hold[n-1,n-1]
Hnew, U = francis_step(Hold,rho)
print("second Francis iteration of degree 1 (rho = ",rho,"), H = U^T H U:")
print(Hnew)
input("Press Enter to continue...")
print("checks:")
print("  ||U^T U - I|| = ", norm(U.T @ U - np.eye(n)))
print("  ||U U^T - I|| = ", norm(U @ U.T - np.eye(n)))
print("  ||U^T H U - Hnew|| = ", norm(U.T @ Hold @ U - Hnew))
print("  subdiagonal of H:")
print(np.diag(Hnew,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(Hnew,-1)),"\n")
input("Press Enter to continue...")

Hold = Hnew.copy()
rho = Hold[n-1,n-1]
Hnew, U = francis_step(Hold,rho)
print("third Francis iteration of degree 1 (rho = ",rho,"), H = U^T H U:")
print(Hnew)
input("Press Enter to continue...")
print("checks:")
print("  ||U^T U - I|| = ", norm(U.T @ U - np.eye(n)))
print("  ||U U^T - I|| = ", norm(U @ U.T - np.eye(n)))
print("  ||U^T H U - Hnew|| = ", norm(U.T @ Hold @ U - Hnew))
print("  subdiagonal of H:")
print(np.diag(Hnew,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(Hnew,-1)),"\n")
input("Press Enter to continue...")

print("Note that convergence is most rapid at bottom, and the last row is essentially triangular.")
print("Let's proceed to a non-symmetric matrix.")
input("Press Enter to continue...")


# non-symmetric matrix
n = 6   # must be even
v = np.random.rand(n,n//2) + np.complex(0,1)*np.random.rand(n,n//2)
V = np.hstack([v, np.conj(v)])
d = np.random.rand(n//2) + np.complex(0,1)*np.random.rand(n//2)
D = np.diag(np.hstack([d, np.conj(d)] ))
A = np.real(V @ D @ np.linalg.inv(V))
print('second matrix A (non-symmetric):')
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
print("  ||Q^H A Q - H|| = ", norm(Q.T.conj() @ A @ Q - H))
print("  subdiagonal of H:")
print(np.diag(H,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(H,-1)), "\n")
input("Press Enter to continue...")


# perform a few Francis steps (Rayleigh quotient shift) to show what happens
Hold = H.copy()
rho = Hold[n-1,n-1]
Hnew, U = francis_step(Hold,rho)
print("first Francis iteration of degree 1 (rho = ",rho,"), H = U^H H U:")
print(Hnew)
input("Press Enter to continue...")
print("checks:")
print("  ||U^H U - I|| = ", norm(U.T.conj() @ U - np.eye(n)))
print("  ||U U^H - I|| = ", norm(U @ U.T.conj() - np.eye(n)))
print("  ||U^H H U - Hnew|| = ", norm(U.T.conj() @ Hold @ U - Hnew))
print("  subdiagonal of H:")
print(np.diag(Hnew,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(Hnew,-1)),"\n")
input("Press Enter to continue...")

Hold = Hnew.copy()
rho = Hold[n-1,n-1]
Hnew, U = francis_step(Hold,rho)
print("second Francis iteration of degree 1 (rho = ",rho,"), H = U^H H U:")
print(Hnew)
input("Press Enter to continue...")
print("checks:")
print("  ||U^H U - I|| = ", norm(U.T.conj() @ U - np.eye(n)))
print("  ||U U^H - I|| = ", norm(U @ U.T.conj() - np.eye(n)))
print("  ||U^H H U - Hnew|| = ", norm(U.T.conj() @ Hold @ U - Hnew))
print("  subdiagonal of H:")
print(np.diag(Hnew,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(Hnew,-1)),"\n")
input("Press Enter to continue...")

Hold = Hnew.copy()
rho = Hold[n-1,n-1]
Hnew, U = francis_step(Hold,rho)
print("third Francis iteration of degree 1 (rho = ",rho,"), H = U^H H U:")
print(Hnew)
input("Press Enter to continue...")
print("checks:")
print("  ||U^H U - I|| = ", norm(U.T.conj() @ U - np.eye(n)))
print("  ||U U^H - I|| = ", norm(U @ U.T.conj() - np.eye(n)))
print("  ||U^H H U - Hnew|| = ", norm(U.T.conj() @ Hold @ U - Hnew))
print("  subdiagonal of H:")
print(np.diag(Hnew,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(Hnew,-1)),"\n")
input("Press Enter to continue...")

print("Note that convergence stagnates (due to complex eigenvalues but real arithmetic).")
print("This used the Rayleigh quotient shift; let's try Wilkinson.")
input("Press Enter to continue...")




# non-symmetric matrix, Wilkinson shift
print("recall our upper-Hessenberg H = Q^H A Q:")
print(H)
print("  subdiagonal of H:")
print(np.diag(H,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(H,-1)), "\n")
input("Press Enter to continue...")


# perform a few Francis steps (Wilkinson shift) to show what happens
Hold = H.copy()
s, vecs = np.linalg.eig(Hold[n-2:,n-2:])
print("s = ",s)
if (np.abs(s[0]-Hold[n-1,n-1]) < np.abs(s[1]-Hold[n-1,n-1])):
   rho = s[0]
else:
   rho = s[1]
Hnew, U = francis_step(Hold,rho)
print("first Francis iteration of degree 1:  shift =", rho)
print("  H = U^H H U:")
print(Hnew)
input("Press Enter to continue...")
print("checks:")
print("  ||U^H U - I|| = ", norm(U.T.conj() @ U - np.eye(n)))
print("  ||U U^H - I|| = ", norm(U @ U.T.conj() - np.eye(n)))
print("  ||U^H H U - Hnew|| = ", norm(U.T.conj() @ Hold @ U - Hnew))
print("  subdiagonal of H:")
print(np.diag(Hnew,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(Hnew,-1)),"\n")
input("Press Enter to continue...")

Hold = Hnew.copy()
s, vecs = np.linalg.eig(Hold[n-2:,n-2:])
if (np.abs(s[0]-Hold[n-1,n-1]) < np.abs(s[1]-Hold[n-1,n-1])):
   rho = s[0]
else:
   rho = s[1]
Hnew, U = francis_step(Hold,rho)
print("second Francis iteration of degree 1:  shift =", rho)
print("  H = U^H H U:")
print(Hnew)
input("Press Enter to continue...")
print("checks:")
print("  ||U^H U - I|| = ", norm(U.T.conj() @ U - np.eye(n)))
print("  ||U U^H - I|| = ", norm(U @ U.T.conj() - np.eye(n)))
print("  ||U^H H U - Hnew|| = ", norm(U.T.conj() @ Hold @ U - Hnew))
print("  subdiagonal of H:")
print(np.diag(Hnew,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(Hnew,-1)),"\n")
input("Press Enter to continue...")

Hold = Hnew.copy()
s, vecs = np.linalg.eig(Hold[n-2:,n-2:])
if (np.abs(s[0]-Hold[n-1,n-1]) < np.abs(s[1]-Hold[n-1,n-1])):
   rho = s[0]
else:
   rho = s[1]
Hnew, U = francis_step(Hold,rho)
print("third Francis iteration of degree 1:  shift =", rho)
print("  H = U^H H U:")
print(Hnew)
input("Press Enter to continue...")
print("checks:")
print("  ||U^H U - I|| = ", norm(U.T.conj() @ U - np.eye(n)))
print("  ||U U^H - I|| = ", norm(U @ U.T.conj() - np.eye(n)))
print("  ||U^H H U - Hnew|| = ", norm(U.T.conj() @ Hold @ U - Hnew))
print("  subdiagonal of H:")
print(np.diag(Hnew,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(Hnew,-1)),"\n")
input("Press Enter to continue...")

Hold = Hnew.copy()
s, vecs = np.linalg.eig(Hold[n-2:,n-2:])
if (np.abs(s[0]-Hold[n-1,n-1]) < np.abs(s[1]-Hold[n-1,n-1])):
   rho = s[0]
else:
   rho = s[1]
Hnew, U = francis_step(Hold,rho)
print("fourth Francis iteration of degree 1:  shift =", rho)
print("  H = U^H H U:")
print(Hnew)
input("Press Enter to continue...")
print("checks:")
print("  ||U^H U - I|| = ", norm(U.T.conj() @ U - np.eye(n)))
print("  ||U U^H - I|| = ", norm(U @ U.T.conj() - np.eye(n)))
print("  ||U^H H U - Hnew|| = ", norm(U.T.conj() @ Hold @ U - Hnew))
print("  subdiagonal of H:")
print(np.diag(Hnew,-1))
print("  norm of subdiagonal of H = ", norm(np.diag(Hnew,-1)),"\n")
input("Press Enter to continue...")

print("Convergence is restored.")
