function b = matvec_col(A, x) 
% Usage: b = matvec_col(A, x) 
%
% Function to perform column-oriented matrix multiplication
% 
% Inputs: 
%    A is a matrix (m x n)
%    x is a vector (n)
%
% Outputs:
%    b is a vector (m)
%
% Daniel R. Reynolds
% SMU Mathematics
% Math 5316
% Spring 2019

% get problem dimensions
[m,n] = size(A);

% check that A and x are compatible
if (n ~= length(x))
   error('matvec_col error: A and x are incompatible')
end

% initialize output
b = zeros(m,1);

% perform product
for j = 1:n
   for i=1:m
      b(i) = b(i) + A(i,j)*x(j);
   end
end

