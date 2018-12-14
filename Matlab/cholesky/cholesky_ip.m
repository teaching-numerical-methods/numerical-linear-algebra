function A = cholesky_ip(A)
% Usage: A = cholesky_ip(A) 
%
% Function to perform an inner-product-oriented Cholesky factorization
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
   error('cholesky_ip error: A is not square')
end

% perform factorization in-place
for i=1:n                        % loop over rows of result
   for k=1:i-1                   % update diagonal
      A(i,i) = A(i,i) - A(k,i)*A(k,i);
   end
   if (A(i,i) <= 0)              % check positive definite
      error('cholesky_ip error: A is not positive definite')
   end
   A(i,i) = sqrt(A(i,i));        % set diagonal entry for row
   for j=i+1:n                   % loop over remainder of row
      for k=1:i-1                % update row entry
         A(i,j) = A(i,j) - A(k,i)*A(k,j);   
      end
      A(i,j) = A(i,j) / A(i,i);  % set row entry
   end
end

