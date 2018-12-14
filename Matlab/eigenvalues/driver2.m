% Script to show what happens with successive Francis iterations on a matrix.
%
% Daniel R. Reynolds
% Math 5316 @ SMU
% Spring 2019

clear


% symmetric matrix
n = 10;
A = rand(n,n);
A = A+A';
fprintf('first matrix A (symmetric):\n')
disp(A);
pause

% convert to upper-Hessenberg form
[H,Q] = upper_hess(A);
fprintf('upper-Hessenberg H = Q''AQ:\n');
disp(H);
pause
fprintf('checks:\n');
fprintf('  norm of lower portion of H = %g\n', norm(H - triu(H,-1)));
fprintf('  ||Q''Q - I|| = %g\n', norm(Q'*Q-eye(n)));
fprintf('  ||QQ'' - I|| = %g\n', norm(Q*Q'-eye(n)));
fprintf('  ||Q''AQ - H|| = %g\n', norm(Q'*A*Q-H));
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(H,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(H,-1)));
pause


% perform a few Francis steps (Rayleigh quotient shift) to show what happens
Hold = H;
rho = Hold(n,n);
[Hnew,U] = francis_step(Hold,rho);
fprintf('first Francis iteration of degree 1 (rho = %g): H = U''HU:\n',rho);
disp(Hnew);
pause
fprintf('checks:\n');
fprintf('  ||U''U - I|| = %g\n', norm(U'*U-eye(n)));
fprintf('  ||UU'' - I|| = %g\n', norm(U*U'-eye(n)));
fprintf('  ||U''HU - Hnew|| = %g\n', norm(U'*Hold*U-Hnew));
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(Hnew,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(Hnew,-1)));
pause

Hold = Hnew;
rho = Hold(n,n);
[Hnew,U] = francis_step(Hold,rho);
fprintf('second Francis iteration of degree 1 (rho = %g): H = U''HU:\n',rho);
disp(Hnew);
pause
fprintf('checks:\n');
fprintf('  ||U''U - I|| = %g\n', norm(U'*U-eye(n)));
fprintf('  ||UU'' - I|| = %g\n', norm(U*U'-eye(n)));
fprintf('  ||U''HU - Hnew|| = %g\n', norm(U'*Hold*U-Hnew));
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(Hnew,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(Hnew,-1)));
pause

Hold = Hnew;
rho = Hold(n,n);
[Hnew,U] = francis_step(Hold,rho);
fprintf('third Francis iteration of degree 1 (rho = %g): H = U''HU:\n',rho);
disp(Hnew);
pause
fprintf('checks:\n');
fprintf('  ||U''U - I|| = %g\n', norm(U'*U-eye(n)));
fprintf('  ||UU'' - I|| = %g\n', norm(U*U'-eye(n)));
fprintf('  ||U''HU - Hnew|| = %g\n', norm(U'*Hold*U-Hnew));
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(Hnew,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(Hnew,-1)));
pause

fprintf('Note that convergence is most rapid at bottom, and the last row is essentially triangular.\n\n');
fprintf('Let''s proceed to a non-symmetric matrix.\n\n');
pause


% non-symmetric matrix
n = 6;   % must be even
v = rand(n,n/2) + i*rand(n,n/2);
V = [v, conj(v)];
d = rand(n/2,1) + i*rand(n/2,1);
D = diag([d; conj(d)]);
A = real(V*D*inv(V));
fprintf('second matrix A (non-symmetric):\n')
disp(A);
pause

% convert to upper-Hessenberg form
[H,Q] = upper_hess(A);
fprintf('upper-Hessenberg H = Q''AQ:\n');
disp(H);
fprintf('checks:\n');
fprintf('  norm of lower portion of H = %g\n', norm(H - triu(H,-1)));
fprintf('  ||Q''Q - I|| = %g\n', norm(Q'*Q-eye(n)));
fprintf('  ||QQ'' - I|| = %g\n', norm(Q*Q'-eye(n)));
fprintf('  ||Q''AQ - H|| = %g\n', norm(Q'*A*Q-H));
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(H,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(H,-1)));
pause


% perform a few Francis steps (Rayleigh quotient shift) to show what happens
Hold = H;
rho = Hold(n,n);
[Hnew,U] = francis_step(Hold,rho);
fprintf('first Francis iteration of degree 1 (rho = %g): H = U''HU:\n',rho);
disp(Hnew);
pause
fprintf('checks:\n');
fprintf('  ||U''U - I|| = %g\n', norm(U'*U-eye(n)));
fprintf('  ||UU'' - I|| = %g\n', norm(U*U'-eye(n)));
fprintf('  ||U''HU - Hnew|| = %g\n', norm(U'*Hold*U-Hnew));
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(Hnew,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(Hnew,-1)));
pause

Hold = Hnew;
rho = Hold(n,n);
[Hnew,U] = francis_step(Hold,rho);
fprintf('second Francis iteration of degree 1 (rho = %g): H = U''HU:\n',rho);
disp(Hnew);
pause
fprintf('checks:\n');
fprintf('  ||U''U - I|| = %g\n', norm(U'*U-eye(n)));
fprintf('  ||UU'' - I|| = %g\n', norm(U*U'-eye(n)));
fprintf('  ||U''HU - Hnew|| = %g\n', norm(U'*Hold*U-Hnew));
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(Hnew,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(Hnew,-1)));
pause

Hold = Hnew;
rho = Hold(n,n);
[Hnew,U] = francis_step(Hold,rho);
fprintf('third Francis iteration of degree 1 (rho = %g): H = U''HU:\n',rho);
disp(Hnew);
pause
fprintf('checks:\n');
fprintf('  ||U''U - I|| = %g\n', norm(U'*U-eye(n)));
fprintf('  ||UU'' - I|| = %g\n', norm(U*U'-eye(n)));
fprintf('  ||U''HU - Hnew|| = %g\n', norm(U'*Hold*U-Hnew));
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(Hnew,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(Hnew,-1)));
pause

fprintf('Note that convergence stagnates (due to complex eigenvalues but real arithmetic).\n\n');
fprintf('This used the Rayleigh quotient shift; let''s try Wilkinson.\n\n');
pause




% non-symmetric matrix, Wilkinson shift
fprintf('recall our upper-Hessenberg H = Q''AQ:\n');
disp(H);
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(H,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(H,-1)));
pause


% perform a few Francis steps (Wilkinson shift) to show what happens
Hold = H;
s = eig(Hold(n-1:n,n-1:n));
if (abs(s(1)-Hold(n,n)) < abs(s(2)-Hold(n,n)))
   rho = s(1);
else
   rho = s(2);
end
[Hnew,U] = francis_step(Hold,rho);
fprintf('first Francis iteration of degree 1:\n  shift =');
disp(rho)
fprintf('  H = U''HU:\n');
disp(Hnew);
pause
fprintf('checks:\n');
fprintf('  ||U''U - I|| = %g\n', norm(U'*U-eye(n)));
fprintf('  ||UU'' - I|| = %g\n', norm(U*U'-eye(n)));
fprintf('  ||U''HU - Hnew|| = %g\n', norm(U'*Hold*U-Hnew));
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(Hnew,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(Hnew,-1)));
pause

Hold = Hnew;
s = eig(Hold(n-1:n,n-1:n));
if (abs(s(1)-Hold(n,n)) < abs(s(2)-Hold(n,n)))
   rho = s(1);
else
   rho = s(2);
end
[Hnew,U] = francis_step(Hold,rho);
fprintf('second Francis iteration of degree 1:\n  shift =');
disp(rho)
fprintf('  H = U''HU:\n');
disp(Hnew);
pause
fprintf('checks:\n');
fprintf('  ||U''U - I|| = %g\n', norm(U'*U-eye(n)));
fprintf('  ||UU'' - I|| = %g\n', norm(U*U'-eye(n)));
fprintf('  ||U''HU - Hnew|| = %g\n', norm(U'*Hold*U-Hnew));
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(Hnew,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(Hnew,-1)));
pause

Hold = Hnew;
s = eig(Hold(n-1:n,n-1:n));
if (abs(s(1)-Hold(n,n)) < abs(s(2)-Hold(n,n)))
   rho = s(1);
else
   rho = s(2);
end
[Hnew,U] = francis_step(Hold,rho);
fprintf('third Francis iteration of degree 1:\n  shift =');
disp(rho)
fprintf('  H = U''HU:\n');
disp(Hnew);
pause
fprintf('checks:\n');
fprintf('  ||U''U - I|| = %g\n', norm(U'*U-eye(n)));
fprintf('  ||UU'' - I|| = %g\n', norm(U*U'-eye(n)));
fprintf('  ||U''HU - Hnew|| = %g\n', norm(U'*Hold*U-Hnew));
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(Hnew,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(Hnew,-1)));
pause

Hold = Hnew;
s = eig(Hold(n-1:n,n-1:n));
if (abs(s(1)-Hold(n,n)) < abs(s(2)-Hold(n,n)))
   rho = s(1);
else
   rho = s(2);
end
[Hnew,U] = francis_step(Hold,rho);
fprintf('fourth Francis iteration of degree 1:\n  shift =');
disp(rho)
fprintf('  H = U''HU:\n');
disp(Hnew);
pause
fprintf('checks:\n');
fprintf('  ||U''U - I|| = %g\n', norm(U'*U-eye(n)));
fprintf('  ||UU'' - I|| = %g\n', norm(U*U'-eye(n)));
fprintf('  ||U''HU - Hnew|| = %g\n', norm(U'*Hold*U-Hnew));
fprintf('  subdiagonal of H:\n');
format short e
disp(diag(Hnew,-1)')
format short 
fprintf('  norm of subdiagonal of H = %g\n\n', norm(diag(Hnew,-1)));
pause

fprintf('Convergence is restored.\n\n');
