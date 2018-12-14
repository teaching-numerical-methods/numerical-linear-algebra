#!/usr/bin/env python3
#
# Script to demonstrate sparse LU factorizations (and reorderings).
#
# Daniel R. Reynolds
# SMU Mathematics
# Math 5316
# Spring 2019

# imports
import time
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import linalg as la
import matplotlib.pyplot as plt


##################
# utility routines

def makeplots(D):
    """
    Usage: makeplots(D)

    This routine creates 3 plots:
       1. D and its L,U factors
       2. D using the approximate minimum degree reordering, and its
          L,U factors)
       3. D using the reverse Cuthill-McKee reordering, and its L,U
       factors)
    """

    # get size of D
    m,n = D.shape
    
    # figure 1: original structure and L,U factors
    lu = la.splu(D, permc_spec='NATURAL')
    LU = lu.L+lu.U
    fig, axarr = plt.subplots(1,2)
    axarr[0].spy(D)
    axarr[0].set_title('Original matrix (nnz = ' + str(D.getnnz()) + ')')
    axarr[1].spy(LU.A)
    axarr[1].set_title('Original: L+U (nnz = ' + str(LU.getnnz()) + ')')


    # figure 2: approximate minimum degree on D.T@D
    lu = la.splu(D, permc_spec='MMD_ATA')
    LU = lu.L+lu.U
    Pr = csc_matrix((n,n))
    Pr[lu.perm_r, np.arange(n)] = 1
    Pc = csc_matrix((n,n))
    Pc[np.arange(n), lu.perm_c] = 1
    PD = Pr * D * Pc
    fig, axarr = plt.subplots(1,2)
    axarr[0].spy(PD.A)
    axarr[0].set_title('Approximate minimum degree on A.T @ A (nnz = ' + str(PD.getnnz()) + ')')
    axarr[1].spy(LU.A)
    axarr[1].set_title('AMD A.T @ A: L+U (nnz = ' + str(LU.getnnz()) + ')')


    # figure 3: approximate minimum degree on A.T+A
    lu = la.splu(D, permc_spec='MMD_AT_PLUS_A')
    LU = lu.L+lu.U
    Pr = csc_matrix((n,n))
    Pr[lu.perm_r, np.arange(n)] = 1
    Pc = csc_matrix((n,n))
    Pc[np.arange(n), lu.perm_c] = 1
    PD = Pr * D * Pc
    fig, axarr = plt.subplots(1,2)
    axarr[0].spy(PD.A)
    axarr[0].set_title('Approximate minimum degree on A.T + A (nnz = ' + str(PD.getnnz()) + ')')
    axarr[1].spy(LU.A)
    axarr[1].set_title('AMD A.T + A: L+U (nnz = ' + str(LU.getnnz()) + ')')


    # figure 4: approximate minimum degree on column ordering
    fig, axarr = plt.subplots(1,2)
    lu = la.splu(D, permc_spec='COLAMD')
    LU = lu.L+lu.U
    Pr = csc_matrix((n,n))
    Pr[lu.perm_r, np.arange(n)] = 1
    Pc = csc_matrix((n,n))
    Pc[np.arange(n), lu.perm_c] = 1
    PD = (Pr * D * Pc)
    axarr[0].spy(PD.A)
    axarr[0].set_title('Approximate minimum degree column order (nnz = ' + str(PD.getnnz()) + ')')
    axarr[1].spy(LU.A)
    axarr[1].set_title('COLAMD: L+U (nnz = ' + str(LU.getnnz()) + ')')

    plt.show()
    


def diff_2D(Nx,Ny):
    """
    Usage: D = diff_2D(Nx,Ny)
    
    This routine creates the diffusion matrix resulting from the equation
    \[
         u - \Delta u,
    \]
    where $u \in \Real$ is defined on the square domain [0,1] x [0,1], which
    is discretized using Nx points in the x-direction, and Ny points in the 
    y-direction, and the Laplace operator is discretized using the standard 
    2nd-order 5 point stencil.  Homogeneous Dirichlet boundary conditions are
    assumed just outside the domain.
    
    inputs:
        Nx       # spatial points in the x-direction of the domain
        Ny       # spatial points in the y-direction of the domain
    
    outputs:
        D        REAL (Nx*Ny) x (Nx*Ny) sparse (CSC) matrix
    """
    
    # set indexing function from 2D physical space to 1D index space
    def ij(i,j):
        return(j*Nx + i)
    
    # initialize the output matrix 
    D = np.zeros((Nx*Ny,Nx*Ny), dtype=float)

    # set differencing constants
    dx = 1.0/(Nx-1)
    dy = 1.0/(Ny-1)
    Dx2i = 1.0/dx/dx
    Dy2i = 1.0/dy/dy
    Diag = 1.0 + 2.0*(Dx2i + Dy2i)

    # iterate over the domain
    for iy in range(Ny):
        for ix in range(Nx):

            # set the matrix entries for this row of D
            D[ ij(ix,iy), ij(ix,iy) ] = Diag
            if (ix > 0):
                D[ ij(ix,iy), ij(ix-1,iy) ] = -Dx2i
            if (ix < Nx-1):
                D[ ij(ix,iy), ij(ix+1,iy) ] = -Dx2i
            if (iy > 0):
                D[ ij(ix,iy), ij(ix,iy-1) ] = -Dy2i
            if (iy < Ny-1):
                D[ ij(ix,iy), ij(ix,iy+1) ] = -Dy2i

    Dcsc = csc_matrix(D)
    return Dcsc


def diff_3D(Nx,Ny,Nz):
    """
    Usage: D = diff_3D(Nx,Ny,Nz)
    
    This routine creates the diffusion matrix resulting from the equation
    \[
         u - \Delta u,
    \]
    where $u \in \Real$ is defined on the cube domain [0,1] x [0,1] x [0,1], 
    which is discretized using Nx points in the x-direction, Ny points in the 
    y-direction, Nz points in the z-direction, and the Laplace operator is 
    discretized using the standard 2nd-order 7 point stencil.  Homogeneous 
    Dirichlet boundary conditions are assumed just outside the domain.
    
    inputs:
        Nx       # spatial points in the x-direction of the domain
        Ny       # spatial points in the y-direction of the domain
        Nz       # spatial points in the z-direction of the domain
    
    outputs:
        D        REAL (Nx*Ny*Nz) x (Nx*Ny*Nz) sparse (CSC) matrix
    """

    # set indexing function from 3D physical space to 1D index space
    def ijk(i,j,k):
        return (k*Nx*Ny + j*Nx + i)

    # initialize the output matrix
    D = np.zeros((Nx*Ny*Nz,Nx*Ny*Nz), dtype=float)

    # set differencing constants
    dx = 1.0/(Nx-1)
    dy = 1.0/(Ny-1)
    dz = 1.0/(Nz-1)
    Dx2i = 1.0/dx/dx
    Dy2i = 1.0/dy/dy
    Dz2i = 1.0/dz/dz
    Diag = 1.0 + 2.0*(Dx2i + Dy2i + Dz2i)

    # iterate over the domain
    for iz in range(Nz):
        for iy in range(Ny):
            for ix in range(Nx):
     
                # set the matrix entries for this row of D
                D[ ijk(ix,iy,iz), ijk(ix,iy,iz) ] = Diag
                if (ix > 0):
                    D[ ijk(ix,iy,iz), ijk(ix-1,iy,iz) ] = -Dx2i
                if (ix < Nx-1):
                    D[ ijk(ix,iy,iz), ijk(ix+1,iy,iz) ] = -Dx2i
                if (iy > 0):
                    D[ ijk(ix,iy,iz), ijk(ix,iy-1,iz) ] = -Dy2i
                if (iy < Ny-1):
                    D[ ijk(ix,iy,iz), ijk(ix,iy+1,iz) ] = -Dy2i
                if (iz > 0):
                    D[ ijk(ix,iy,iz), ijk(ix,iy,iz-1) ] = -Dz2i
                if (iz < Nz-1):
                    D[ ijk(ix,iy,iz), ijk(ix,iy,iz+1) ] = -Dz2i

    Dcsc = csc_matrix(D)
    return Dcsc



##################
# script

if __name__ == "__main__":

    # problem 1: small 2D diffusion matrix
    print("problem 1: small 2D diffusion matrix")
    D = diff_2D(5,10)
    makeplots(D)
    input("Press enter to continue")

    # problem 2: larger 2D diffusion matrix
    print("problem 2: larger 2D diffusion matrix");
    D = diff_2D(50,100)
    makeplots(D)
    input("Press enter to continue")

    # problem 3: small 3D diffusion matrix
    print("problem 3: small 3D diffusion matrix");
    D = diff_3D(5,8,10)
    makeplots(D)
    input("Press enter to continue")

    # problem 4: larger 3D diffusion matrix
    print("problem 4: larger 3D diffusion matrix");
    D = diff_3D(20,25,30)
    makeplots(D)

