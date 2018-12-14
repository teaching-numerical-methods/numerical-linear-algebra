function y = fwdsub_col(L, b) 
% Usage: y = fwdsub_col(L, b) 
%
% Function to perform column-oriented forwards substitution
% 
% Inputs: 
%    L is a lower-triangular matrix (n x n)
%    b is a vector (n)
%
% Outputs:
%    y is a vector (n)
%
% Daniel R. Reynolds
% SMU Mathematics
% Math 5316
% Spring 2019

% get problem dimensions
[m,n] = size(L);

% check that L and b are compatible
if (m ~= n)
   error('fwdsub_col error: L is not square')
end
if (n ~= size(b))
   error('fwdsub_col error: L and b are incompatible')
end

% check that L is nonsingular
for i=1:n
   if (L(i,i) == 0)
      error('fwdsub_col error: L is singular');
   end
end

% copy b into solution vector
y = b;

% loop over matrix columns, performing solve
for j=1:n
   y(j) = y(j)/L(j,j);    % solve row
   for i=j+1:n            % update all remaining rhs
      y(i) = y(i) - L(i,j)*y(j);
   end
end
