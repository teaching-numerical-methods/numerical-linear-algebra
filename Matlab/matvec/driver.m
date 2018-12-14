% Script to test matrix-vector products
%
% Daniel R. Reynolds
% SMU Mathematics
% Math 5316
% Spring 2019
clear

% set testing values
mvals = [1000, 2000, 4000, 8000];
nvals = [2000, 4000, 8000, 16000];
nsizes = 4;

% run tests
for k=1:nsizes

   % set the problem size
   m = mvals(k);
   n = nvals(k);

   % display current problem size
   fprintf('\nTesting matrix-vector products with a %i x %i matrix:\n', m, n);

   % allocate the matrix & vector of this size
   A = zeros(m,n);
   x = zeros(n,1);

   % fill A and x with values
   for i = 1:m
      for j = 1:n
         A(i,j) = (1 + i - j)/(n+m);
      end
   end
   for j = 1:n
      x(j) = j/n;
   end
   
   % perform row-based product
   stime = tic;
   b = matvec_row(A, x);
   runtime = toc(stime);
   b_err = norm(b-A*x);
   fprintf('   matvec_row:  time = %.4f, error = %.4e\n', runtime, b_err);
   
   % perform column-based product
   stime = tic;
   b = matvec_col(A, x);
   runtime = toc(stime);
   b_err = norm(b-A*x);
   fprintf('   matvec_col:  time = %.4f, error = %.4e\n', runtime, b_err);
   
   % call '*' operator for product
   stime = tic;
   b = A*x;
   runtime = toc(stime);
   b_err = norm(b-A*x);
   fprintf('   * operator:  time = %.4f, error = %.4e\n', runtime, b_err);

end
