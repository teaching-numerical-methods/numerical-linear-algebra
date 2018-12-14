function D = diff_2D(Nx,Ny)
% Usage: D = diff_2D(Nx,Ny)
%
% This routine creates the diffusion matrix resulting from the equation
% \[
%      u - \Delta u,
% \]
% where $u \in \Real$ is defined on the square domain [0,1] x [0,1], which
% is discretized using Nx points in the x-direction, and Ny points in the 
% y-direction, and the Laplace operator is discretized using the standard 
% 2nd-order 5 point stencil.  Homogeneous Dirichlet boundary conditions are
% assumed just outside the domain.
%
% inputs:
%     Nx       # spatial points in the x-direction of the domain
%     Ny       # spatial points in the y-direction of the domain
%
% outputs:
%     D        REAL (Nx*Ny) x (Nx*Ny) sparse matrix
%
% Daniel R. Reynolds
% Math5316 @ SMU
% Spring 2019

% set indexing function from 2D physical space to 1D index space
ij = @(i,j) (j-1)*Nx + i;

% initialize the output matrix 
D = sparse(Nx*Ny,Nx*Ny);

% set differencing constants
dx = 1/(Nx-1);
dy = 1/(Ny-1);
Dx2i = 1/dx/dx;
Dy2i = 1/dy/dy;
Diag = 1 + 2*(Dx2i + Dy2i);

% iterate over the domain
for iy=1:Ny
   for ix=1:Nx
     
      % set the matrix entries for this row of D
      D( ij(ix,iy), ij(ix,iy) ) = Diag;
      if (ix > 1)
	 D( ij(ix,iy), ij(ix-1,iy) ) = -Dx2i;
      end
      if (ix < Nx)
	 D( ij(ix,iy), ij(ix+1,iy) ) = -Dx2i;
      end
      if (iy > 1)
	 D( ij(ix,iy), ij(ix,iy-1) ) = -Dy2i;
      end
      if (iy < Ny)
	 D( ij(ix,iy), ij(ix,iy+1) ) = -Dy2i;
      end
      
   end
end


% end of function