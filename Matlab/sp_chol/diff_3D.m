function D = diff_3D(Nx,Ny,Nz)
% Usage: D = diff_3D(Nx,Ny,Nz)
%
% This routine creates the diffusion matrix resulting from the equation
% \[
%      u - \Delta u,
% \]
% where $u \in \Real$ is defined on the cube domain [0,1] x [0,1] x [0,1], 
% which is discretized using Nx points in the x-direction, Ny points in the 
% y-direction, Nz points in the z-direction, and the Laplace operator is 
% discretized using the standard 2nd-order 7 point stencil.  Homogeneous 
% Dirichlet boundary conditions are assumed just outside the domain.
%
% inputs:
%     Nx       # spatial points in the x-direction of the domain
%     Ny       # spatial points in the y-direction of the domain
%     Nz       # spatial points in the z-direction of the domain
%
% outputs:
%     D        REAL (Nx*Ny*Nz) x (Nx*Ny*Nz) sparse matrix
%
% Daniel R. Reynolds
% Math5316 @ SMU
% Spring 2019


% set indexing function from 3D physical space to 1D index space
ijk = @(i,j,k) (k-1)*Nx*Ny + (j-1)*Nx + i;

% initialize the output matrix
D = sparse(Nx*Ny*Nz,Nx*Ny*Nz);

% set differencing constants
dx = 1/(Nx-1);
dy = 1/(Ny-1);
dz = 1/(Nz-1);
Dx2i = 1/dx/dx;
Dy2i = 1/dy/dy;
Dz2i = 1/dz/dz;
Diag = 1 + 2*(Dx2i + Dy2i + Dz2i);

% iterate over the domain
for iz=1:Nz
   for iy=1:Ny
      for ix=1:Nx
     
	 % set the matrix entries for this row of D
	 D( ijk(ix,iy,iz), ijk(ix,iy,iz) ) = Diag;
	 if (ix > 1)
	    D( ijk(ix,iy,iz), ijk(ix-1,iy,iz) ) = -Dx2i;
	 end
	 if (ix < Nx)	 
	    D( ijk(ix,iy,iz), ijk(ix+1,iy,iz) ) = -Dx2i;
	 end
	 if (iy > 1)
	    D( ijk(ix,iy,iz), ijk(ix,iy-1,iz) ) = -Dy2i;
	 end
	 if (iy < Ny)
	    D( ijk(ix,iy,iz), ijk(ix,iy+1,iz) ) = -Dy2i;
	 end
	 if (iz > 1)
	    D( ijk(ix,iy,iz), ijk(ix,iy,iz-1) ) = -Dz2i;
	 end
	 if (iz < Nz)
	    D( ijk(ix,iy,iz), ijk(ix,iy,iz+1) ) = -Dz2i;
	 end
	 
      end
   end
end


% end of function