function A = cholesky_op(A)
% Usage: A = cholesky_op(A) 
%
% Function to perform an outer-product-oriented Cholesky factorization
% 
% Inputs: 
%    A is a symmetric matrix (n x n)
%
% Outputs:
%    A is the Cholesky-factored matrix, A = R'*R, with R embedded
%    in the upper-triangular portion of A
%
% Daniel R. Reynolds
% SMU Mathematics
% Math 5316
% Spring 2019

% get problem dimensions
[m,n] = size(A);

% check that A is square
if (m ~= n)
   error('cholesky_op error: A is not square')
end

% perform factorization in-place
for i=1:n                        % loop over rows of result
   if (A(i,i) <= 0)              % check positive definite
      error('cholesky_op error: A is not positive definite')
   end
   A(i,i) = sqrt(A(i,i));        % set diagonal entry for row
   for j = i+1:n                 % update row
      A(i,j) = A(i,j) / A(i,i);
   end
   for k = i+1:n                 % update remainder of matrix
      for j = i+1:n
         A(k,j) = A(k,j) - A(i,k)*A(i,j);   
      end
   end
end
