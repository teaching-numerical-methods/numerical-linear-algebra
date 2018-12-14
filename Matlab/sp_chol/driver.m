% Script to demonstrate sparse Cholesky factorizations (and reorderings).
%
% Daniel R. Reynolds
% Math5316 @ SMU
% Spring 2019

clear

% problem 1: small 2D diffusion matrix
fprintf('problem 1: small 2D diffusion matrix\n');
Nx = 5;
Ny = 10;
D = diff_2D(Nx,Ny);
makeplots(D)
pause

% problem 2: larger 2D diffusion matrix
fprintf('problem 2: larger 2D diffusion matrix\n');
Nx = 50;
Ny = 100;
D = diff_2D(Nx,Ny);
makeplots(D)
pause

% problem 3: small 3D diffusion matrix
fprintf('problem 3: small 3D diffusion matrix\n');
Nx = 5;
Ny = 8;
Nz = 10;
D = diff_3D(Nx,Ny,Nz);
makeplots(D)
pause

% problem 4: larger 3D diffusion matrix
fprintf('problem 4: larger 3D diffusion matrix\n');
Nx = 20;
Ny = 25;
Nz = 30;
D = diff_3D(Nx,Ny,Nz);
makeplots(D)


% end of script