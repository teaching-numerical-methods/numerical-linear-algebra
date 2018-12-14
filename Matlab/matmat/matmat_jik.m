function B = matmat_jik(A,X)
% Usage: B = matmat_jik(A,X)
%
% Function to perform matrix-matrix multiplication with "jik" loop ordering.
%
% Inputs: A is a matrix (m x n numpy matrix)
%         X is a matrix (n x p numpy matrix)
% Output: B is a matrix (m x p numpy matrix)
%
% Daniel R. Reynolds
% SMU Mathematics
% Math 5316
% Spring 2019

% get problem dimensions
[m,n]  = size(A);
[n2,p] = size(X);

% check that A and X are compatible
if (n ~= n2)
   error('matmat_jik error: A and X are incompatible')
end

% initialize output
B = zeros(m,p);

% perform product
for j=1:p
   for i=1:m
      for k=1:n
         B(i,j) = B(i,j) + A(i,k)*X(k,j);
      end
   end
end
