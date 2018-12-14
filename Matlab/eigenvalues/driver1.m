% Script to test out approaches for converting a matrix to upper Hessenberg form.
%
% Daniel R. Reynolds
% Math 5316 @ SMU
% Spring 2019

clear

% non-symmetric matrix
n = 10;
A = rand(n,n);
fprintf('first matrix (non-symmetric):\n')
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
fprintf('  ||Q''AQ - H|| = %g\n\n', norm(Q'*A*Q-H));


% symmetric matrix
n = 10;
A = rand(n,n);
A = A+A';
fprintf('second matrix (symmetric):\n')
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
fprintf('  ||Q''AQ - H|| = %g\n\n', norm(Q'*A*Q-H));


% end of script