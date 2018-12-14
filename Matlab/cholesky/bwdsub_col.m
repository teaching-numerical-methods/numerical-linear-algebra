function x = bwdsub_col(U, y) 
% Usage: x = bwdsub_col(U, y) 
%
% Function to perform column-oriented backwards substitution
% 
% Inputs: 
%    U is an upper-triangular matrix (n x n)
%    y is a vector (n)
%
% Outputs:
%    x is a vector (n)
%
% Daniel R. Reynolds
% SMU Mathematics
% Math 5316
% Spring 2019

% get problem dimensions
[m,n] = size(U);

% check that U and y are compatible
if (m ~= n)
   error('bwdsub_col error: U is not square')
end
if (n ~= size(y))
   error('bwdsub_col error: U and y are incompatible')
end

% check that U is nonsingular
for i=1:n
   if (U(i,i) == 0)
      error('bwdsub_col error: U is singular');
   end
end

% copy y into solution vector
x = y;

% loop over matrix columns, performing solve
for j=n:-1:1
   x(j) = x(j)/U(j,j);    % solve row
   for i=1:j-1            % update remaining rhs
      x(i) = x(i) - U(i,j)*x(j);
   end
end
