% Script to demonstrate sparse LU factorizations (and reorderings).
%
% Daniel R. Reynolds
% Math5316 @ SMU
% Spring 2019

clear

% problem 1: small 2D diffusion matrix
fprintf('problem 1: small 2D diffusion matrix\n');
D = diff_2D(5,10);
makeplots(D)
pause

% problem 2: larger 2D diffusion matrix
fprintf('problem 2: larger 2D diffusion matrix\n');
D = diff_2D(50,100);
makeplots(D)
pause

% problem 3: small 3D diffusion matrix
fprintf('problem 3: small 3D diffusion matrix\n');
D = diff_3D(5,8,10);
makeplots(D)
pause

% problem 4: larger 3D diffusion matrix
fprintf('problem 4: larger 3D diffusion matrix\n');
D = diff_3D(20,25,30);
makeplots(D)


% end of script
