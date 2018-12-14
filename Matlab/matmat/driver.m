% Script to test matrix-matrix products
%
% Daniel R. Reynolds
% SMU Mathematics
% Math 5316
% Spring 2019
clear

% set testing values
mvals = [50, 100, 200, 400];
nvals = [100, 200, 400, 800];
pvals = [75, 150, 300, 600];
nsizes = 4;

% run tests
for k=1:nsizes

   % set the problem size
   m = mvals(k);
   n = nvals(k);
   p = pvals(k);

   % display current problem size
   fprintf('\nTesting matrix-matrix products: m = %i n = %i p = %i\n',m,n,p);

   % allocate the matrix & vectors of this size
   A = zeros(m,n);
   X = zeros(n,p);

   % fill A and X with values
   for i = 1:m
      for j = 1:n
         A(i,j) = (1 + i - j)/(n+m);
      end
   end
   for i = 1:n
      for j = 1:p
         X(i,j) = (1 - i + j)/(n+p);
      end
   end
   
   % perform product 1
   stime = tic;
   B = matmat_ijk(A, X);
   runtime = toc(stime);
   B_err = max(max(abs(B-A*X)));
   fprintf('   matmat_ijk:  time = %.4f,  error = %.5e\n', runtime, B_err);
   
   % perform product 2
   stime = tic;
   B = matmat_ikj(A, X);
   runtime = toc(stime);
   B_err = max(max(abs(B-A*X)));
   fprintf('   matmat_ikj:  time = %.4f,  error = %.5e\n', runtime, B_err);

   % perform product 3
   stime = tic;
   B = matmat_jik(A, X);
   runtime = toc(stime);
   B_err = max(max(abs(B-A*X)));
   fprintf('   matmat_jik:  time = %.4f,  error = %.5e\n', runtime, B_err);

   % perform product 4
   stime = tic;
   B = matmat_jki(A, X);
   runtime = toc(stime);
   B_err = max(max(abs(B-A*X)));
   fprintf('   matmat_jki:  time = %.4f,  error = %.5e\n', runtime, B_err);


   % perform product 5
   stime = tic;
   B = matmat_kij(A, X);
   runtime = toc(stime);
   B_err = max(max(abs(B-A*X)));
   fprintf('   matmat_kij:  time = %.4f,  error = %.5e\n', runtime, B_err);


   % perform product 6
   stime = tic;
   B = matmat_kji(A, X);
   runtime = toc(stime);
   B_err = max(max(abs(B-A*X)));
   fprintf('   matmat_kji:  time = %.4f,  error = %.5e\n', runtime, B_err);


   % call '*' operator for product
   stime = tic;
   B = A * X;
   runtime = toc(stime);
   B_err = max(max(abs(B-A*X)));
   fprintf('   * operator:  time = %.4f,  error = %.5e\n', runtime, B_err);

end
