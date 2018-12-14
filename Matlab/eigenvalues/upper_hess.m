function [H,Q] = upper_hess(A)
% Usage: [H,Q] = upper_hess(A)
%
% Function to convert a matrix to upper-Hessenberg form via the unitary
% similarity transformation
%        H = Q'*A*Q
% where A is a general square matrix, Q is a unitary transformation matrix,
% and H is the upper Hessenberg result.
%
% Note: since this routine returns the full transformation matrix Q (and not
% just its parts), it is less efficient than it *could* be.
%
% Input:    A - square matrix
% Outputs:  Q - unitary matrix
%           H - upper Hessenberg matrix
%
% Daniel R. Reynolds
% Math 5316 @ SMU
% Spring 2019

% ensure that A is square
[m,n] = size(A);
if (m ~= n)
   error('upper_hess: matrix must be square')
end

% initialize results
Q = eye(n);
H = A;

% iterate over columns
for j=1:n-2

   % construct reflector to zero out column below first subdiagonal
   b = H(j+1:n,j);               % column to transform
   tau = norm(b)*sign(b(1)+eps); % norm of column
   u = b/(b(1)+tau);             % construct u = (b-y)/(b1+tau), where y = [-tau,0,...,0]'
   u(1) = 1;
   gamma = (tau+b(1))/tau;
   Qj = eye(n-j) - gamma*u*u';  % create Householder matrix

   % multiply remaining submatrices of H and Q by Qj
   H(j+1:n,:) = Qj*H(j+1:n,:);
   H(:,j+1:n) = H(:,j+1:n)*Qj;
   Q(:,j+1:n) = Q(:,j+1:n)*Qj;
   
end


% end function