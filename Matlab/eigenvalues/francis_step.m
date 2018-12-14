function [A,Q] = francis_step(A,rho)
% Usage: [A,Q] = francis_step(A,rho)
%
% Function to perform one iteration of Francis' algorithm of degree one.
% Computes the unitary similarity transformation 
%        Anew = Q'*Aold*Q
% where Aold is an upper-Hessenberg matrix, Q is a unitary transformation matrix,
% and Anew is an upper-Hessenberg result, through introducing a bulge in Aold (via
% a rotator), and chasing it back out (via a sequence of rotators).
%
% Note: since this routine returns the full transformation matrix Q (and not
% just its parts), it is less efficient than it *could* be.
%
% Input:    A - upper-Hessenberg square matrix (Aold)
%           rho - desired shift
% Outputs:  Q - unitary matrix
%           A - upper-Hessenberg matrix result (Anew)
%
% Daniel R. Reynolds
% Math 5316 @ SMU
% Spring 2019

% get matrix size
n = size(A,1);   % assumed square

% shift the input matrix by rho
A = A - rho*eye(n);

%----------
% introduce bulge 
[c,s] = rot(A(1,1), A(2,1));     % compute rotation matrix components

% apply rotation to A on left: A = Qij'*A
[A(1,:),A(2,:)] = apply_rot_left(A(1,:),A(2,:),c,s);

% apply rotation to A on right: A = A*Qij
[A(:,1),A(:,2)] = apply_rot_right(A(:,1),A(:,2),c,s);

% create Q
Q = eye(n);
Q([1:2],[1:2]) = [c -s; s c];

  
%----------
% chase the bulge
for j=1:n-2

   % construct rotator to zero non-Hessenberg part of column
   [c,s] = rot(A(j+1,j), A(j+2,j));

   % apply rotation to A on left: A = Qij'*A
   [A(j+1,:),A(j+2,:)] = apply_rot_left(A(j+1,:),A(j+2,:),c,s);

   % apply rotation to A on right: A = A*Qij
   [A(:,j+1),A(:,j+2)] = apply_rot_right(A(:,j+1),A(:,j+2),c,s);

   % update Q
   [Q(:,j+1),Q(:,j+2)] = apply_rot_right(Q(:,j+1),Q(:,j+2),c,s);

end

% undo shift on result
A = A + rho*eye(n);

return

% end of function




%----------------------------------------
% nested functions

  function [c,s] = rot(x1,x2)
  % Usage: [c,s] = rot(x1,x2)
  %
  % Function to compute cosine and sine components of a rotation matrix,
  % given the two arguments:
  %
  % Input:   x1 - value to use in elimination
  %          x2 - value to be eliminated
  % Outputs:  c - cosine term
  %           s - sine term
    beta = max(abs(x1),abs(x2));
    if (beta == 0)  % compute c, s values
       c = 1;
       s = 0;
    else
       x1 = x1/beta;
       x2 = x2/beta;
       v = sqrt(x1^2+x2^2);
       c = x1/v;
       s = x2/v;
    end
  end  % end of rot function

  
  function [QMi,QMj] = apply_rot_left(Mi,Mj,c,s)
  % Usage: [QMi,QMj] = apply_rot_left(Mi,Mj,c,s)
  %
  % Applies a rotation matrix multiply on the left,
  %           QM = Q'*M
  % where Q is a rotation matrix with block [c -s; s c]
  %
  % Input:    Mi - first row from matrix to be multiplied
  %           Mj - second row from matrix to be multiplied
  %            c - cosine term in rotator
  %            s - sine term in rotator
  % Outputs: QMi - first row of resulting multiply
  %          QMj - second row of resulting multiply
    QMi = c*Mi + s*Mj;
    QMj = c*Mj - s*Mi;

  end  % end of apply_rot_left function

  
  function [MQi,MQj] = apply_rot_right(Mi,Mj,c,s)
  % Usage: [MQi,MQj] = apply_rot_right(Mi,Mj,c,s)
  %
  % Applies a rotation matrix multiply on the right
  %          MQ = M*Q
  % where Q is a rotation matrix with block [c -s; s c]
  %
  % Input:    Mi - first column from matrix to be multiplied
  %           Mj - second column from matrix to be multiplied
  %            c - cosine term in rotator
  %            s - sine term in rotator
  % Outputs: MQi - first column of resulting multiply
  %          MQj - second column of resulting multiply
    MQi = c*Mi + s*Mj;
    MQj = c*Mj - s*Mi;

  end  % end of apply_rot_right function


end  % end of francis_step function
