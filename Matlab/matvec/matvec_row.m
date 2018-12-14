function b = matvec_row(A, x) 
% Usage: b = matvec_row(A, x) 
%
% Function to perform row-oriented matrix multiplication
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
   error('matvec_row error: A and x are incompatible')
end

% initialize output
b = zeros(m,1);

% perform product
for i=1:m
   for j=1:n
      b(i) = b(i) + A(i,j)*x(j);
   end
end
