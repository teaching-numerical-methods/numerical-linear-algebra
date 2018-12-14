function makeplots(D)
% Usage: makeplots(D)
%
% This routine creates 3 plots:
%    1. D and its Cholesky factor
%    2. D using the approximate minimum degree reordering, and its
%       Cholesky factor)
%    3. D using the reverse Cuthill-McKee reordering, and its Cholesky
%       factor)
%
% Daniel R. Reynolds
% Math5316 @ SMU
% Spring 2019


%   figure 1: original structure and Cholesky factor
figure(1)
subplot(1,2,1)
spy(D)
title(sprintf('Original matrix (nnz = %i)',nnz(D)), 'FontSize', 12)
R = chol(D);
subplot(1,2,2)
spy(R+R')
title(sprintf('Original Cholesky (nnz = %i)',nnz(R+R')), 'FontSize', 12)


%   figure 2: symmetric approximate minimum degree
figure(2)
p = symamd(D);
Damd = D(p,p);
subplot(1,2,1)
spy(Damd)
title(sprintf('Approximate minimum degree (nnz = %i)',nnz(Damd)), 'FontSize', 12)
R = chol(Damd);
subplot(1,2,2)
spy(R+R')
title(sprintf('AMD Cholesky (nnz = %i)',nnz(R+R')), 'FontSize', 12)


%   figure 3: symmetric reverse Cuthill-McKee 
figure(3)
p = symrcm(D);
Drcm = D(p,p);
subplot(1,2,1)
spy(Drcm)
title(sprintf('Reverse Cuthill-McKee (nnz = %i)',nnz(Drcm)), 'FontSize', 12)
R = chol(Drcm);
subplot(1,2,2)
spy(R+R')
title(sprintf('RCM Cholesky (nnz = %i)',nnz(R+R')), 'FontSize', 12)


% end of script