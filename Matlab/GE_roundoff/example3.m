% Example: step-by-step solution a 3x3 linear system using 4-digit
% precision without row-swapping, first with small pivot, second with
% rescaled row to "fix" small pivot
%
% Daniel R. Reynolds
% Math5316 @ SMU
% Spring 2019

clear
d = 4;

% Set up matrix and right-hand side as augmented linear system
A = [0.002 1.231 2.471 3.704;
     1.196 3.165 2.543 6.904;
     1.475 4.271 2.142 7.888];
format short e
disp('original system (in augmented matrix form):')
disp(A)
pause

% First, perform elimination without row-swapping
%    first stage
for i=2:3
   m = roundsd(A(i,1)/A(1,1),d);
   A(i,2:4) = roundsd(A(i,2:4) - roundsd(m*A(1,2:4),d),d);
   A(i,1) = 0;
end
disp('after one stage of elimination:')
disp(A)
pause
%    second stage
m = roundsd(A(3,2)/A(2,2),d);
A(3,3:4) = roundsd(A(3,3:4) - roundsd(m*A(2,3:4),d),d);
A(3,2) = 0;
disp('after second stage of elimination:')
disp(A)
pause

% compute final solution
x = zeros(3,1);
b = A(:,4);
for i=3:-1:1
   for j=i+1:3
      b(i) = roundsd(b(i) - roundsd(A(i,j)*x(j),d),d);
   end
   x(i) = roundsd(b(i)/A(i,i),d);
end
x1 = x;
disp('final solution:')
disp(x1)


% reset augmented matrix system, this time with rescaled first row
A = [2     1231  2471  3704;
     1.196 3.165 2.543 6.904;
     1.475 4.271 2.142 7.888];

% Second, do the same problem using partial pivoting
disp('Repeating, but now with rescaled first equation:')
disp('original system (in augmented matrix form):')
disp(A)
pause

% First, perform elimination without row-swapping
%    first stage
for i=2:3
   m = roundsd(A(i,1)/A(1,1),d);
   A(i,2:4) = roundsd(A(i,2:4) - roundsd(m*A(1,2:4),d),d);
   A(i,1) = 0;
end
disp('after one stage of elimination:')
disp(A)
pause
%    second stage
m = roundsd(A(3,2)/A(2,2),d);
A(3,3:4) = roundsd(A(3,3:4) - roundsd(m*A(2,3:4),d),d);
A(3,2) = 0;
disp('after second stage of elimination:')
disp(A)
pause

% compute final solution
x = zeros(3,1);
b = A(:,4);
for i=3:-1:1
   for j=i+1:3
      b(i) = roundsd(b(i) - roundsd(A(i,j)*x(j),d),d);
   end
   x(i) = roundsd(b(i)/A(i,i),d);
end
x2 = x;
disp('final solution:')
disp(x2)

% compare results
disp('We now compare the results of these two calculations:')
A1 = [0.002 1.231 2.471;
      1.196 3.165 2.543;
      1.475 4.271 2.142];
b1 = [3.704; 6.904; 7.888];
A2 = [2     1231  2471;
      1.196 3.165 2.543;
      1.475 4.271 2.142];
b2 = [3704; 6.904; 7.888];
disp('first computed solution:')
disp(x1)
disp('second computed solution:')
disp(x2)
pause

disp('residual in first solution:')
disp(b1 - A1*x1)
disp('residual in second solution:')
disp(b2 - A2*x2)
