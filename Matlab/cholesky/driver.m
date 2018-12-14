% Script to test Cholesty factorizations
%
% Daniel R. Reynolds
% SMU Mathematics
% Math 5316
% Spring 2019
clear

% set testing values
nvals = [500, 700, 900, 1100];

% run tests
for n = nvals

   % allocate the matrices & vectors of this size
   A = zeros(n,n);
   Rt = zeros(n,n);
   xtrue = zeros(n,1);
   x = zeros(n,1);
   y = zeros(n,1);
   b = zeros(n,1);

   % fill A and xtrue with values
   for i = 1:n
      for j = 1:n
         A(i,j) = 1/(1 + 5*abs(i - j));
      end
   end
   for i=1:n
      xtrue(i) = (1 - i)/n;
   end

   % start first test
   fprintf('\nTesting Cholesky (IP) factorizations: n = %i\n',n);

   % compute b from A and xtrue
   b = A*xtrue;

   % compute Cholesky decomposition, placing result in R
   stime = tic;
   R = cholesky_ip(A);
   runtime = toc(stime);
   fprintf('   cholesky time = %.4f\n', runtime);

   % fill Rt as the transpose of R
   Rt = R';

   % solve linear system
   stime = tic;
   y = fwdsub_col(Rt, b);
   x = bwdsub_col(R, y);
   runtime = toc(stime);
   fprintf('   solve time = %.4f\n', runtime);

   % check error
   err_norm = norm(x-xtrue);
   fprintf('   solution error = %.4f\n', err_norm);

   % start second test
   fprintf('\nTesting Cholesky (OP) factorizations: n = %i\n',n);

   % compute Cholesky decomposition, storing result in R
   stime = tic;
   R = cholesky_op(A);
   runtime = toc(stime);
   fprintf('   cholesky time = %.4f\n', runtime);

   % fill Rt as the transpose of R
   Rt = R';

   % solve linear system
   stime = tic;
   y = fwdsub_col(Rt, b);
   x = bwdsub_col(R, y);
   runtime = toc(stime);
   fprintf('   solve time = %.4f\n', runtime);

   % check error
   err_norm = norm(x-xtrue);
   fprintf('   solution error = %.4f\n', err_norm);

end
