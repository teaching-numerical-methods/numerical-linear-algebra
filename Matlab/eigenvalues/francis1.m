function [A,Q,its] = francis1(A,maxit,tol,stype,diags)
% Usage: [A,Q,its] = francis1(A,maxit,tol,stype,diags)
%
% Function to perform Francis' algorithm of degree one to compute the
% eigen-decomposition of a given matrix A, through iteratively computing the
% similarity transformation 
%        Anew = Q'*Aold*Q
% where Aold is an upper-Hessenberg matrix, Q is a unitary transformation matrix,
% and Anew is an upper-Hessenberg matrix with decreasing subdiagonal.  At
% convergence, Anew should be [essentially] upper-triangular.  Specifically,
% it will satisfy
%       |Anew(j,j-1)| <= tol*(|Anew(j,j)|+|Anew(j-1,j-1)|),  j=2:n
%
% Note: since this routine returns the full transformation matrix Q (as
% opposed to just the triangular matrix result, it is less efficient than it
% *could* be if only the eigenvalues were desired.
%
% Input:    A - upper-Hessenberg square matrix (Aold)
%       maxit - maximum allowed number of iterations
%         tol - relative eigenvalue tolerance
%       stype - desired shift approach:
%                     0 => Rayleigh quotient shift
%                  else => Wilkinson shift
%       diags - flag to turn on/off diagnostic messages:
%                     0 => off
%                  else => on
% Outputs:  Q - unitary matrix
%           A - upper-triangular matrix (Anew)
%         its - number of iterations taken
%
% Daniel R. Reynolds
% Math 5316 @ SMU
% Spring 2019

% ensure that A is square
[m,n] = size(A);
if (m ~= n)
   error('francis1: matrix must be square')
end

% ensure that A is essentially upper-Hessenberg
tmp = A - triu(A,-1);  % portion of A that should be zero
if (norm(tmp,inf)/norm(A,inf) > 100*eps)
   error('francis1: matrix must be upper-Hessenberg')
end

% initialize results, counter toward completion
Q = eye(n);
m = n;

% perform iteration
for its = 1:maxit

   % determine last proper Hessenberg row
   for m=m:-1:2
      if (abs(A(m,m-1)) > tol*(abs(A(m,m))+abs(A(m-1,m-1))))
	 break
      end
   end
   
   % check for completion
   if ((m == 2) && (abs(A(2,1)) <= tol*(abs(A(2,2))+abs(A(1,1)))))
      break;
   end
   
   % set shift
   if (stype == 0)       % Rayleigh quotient shift
      rho = A(m,m);
   else                  % Wilkinson shift

      % set terms from quadratic equation, a*lambda^2 + b*lambda + c = 0,
      % for eigenvalues of submatrix A(m-1:m,m-1:m)
      a = 1;
      b = -A(m-1,m-1) - A(m,m);
      c = A(m-1,m-1)*A(m,m) - A(m-1,m)*A(m,m-1);

      % compute eigenvalues of submatrix
      s1 = (-b + sqrt(b^2-4*a*c))/(2*a);
      s2 = (-b - sqrt(b^2-4*a*c))/(2*a);

      % set shift as closest eigenvalue to A(m,m) component
      if (abs(s1-A(m,m)) < abs(s2-A(m,m)))
	 rho = s1;
      else
	 rho = s2;
      end
   end

   % output some diagnostic information to the screen
   if (diags)
      if (isreal(rho))
	 fprintf('   iter %4i:  submatrix 1:%i,  shift = %g\n', its, m, rho);
      else
	 fprintf('   iter %4i:  submatrix 1:%i,  shift = %g + %gi\n', its, m, real(rho), imag(rho));
      end
   end   
   
   % perform one iteration of algorithm on remaining submatrix
   [A(1:m,1:m),Qt] = francis_step(A(1:m,1:m),rho);

   % update transformation matrix
   Q(:,1:m) = Q(:,1:m)*Qt;

end

% end function
