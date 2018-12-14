% Example: step-by-step LU factorization of the 4x4 Hilbert matrix, rounding
% each entry to 3 significant digits, versus using full double precision
%
% Daniel R. Reynolds
% Math5316 @ SMU
% Spring 2019

clear

% set number of digits to use in second example
d = 3;

% First, do the 4x4 Hilbert matrix using full precision
%    set up matrix, display
format long
H = hilb(4);
disp('original matrix:')
disp(H)
pause
%    first stage
L = eye(4,4);
for i=2:4
   L(i,1) = H(i,1)/H(1,1);
   H(i,2:4) = H(i,2:4) - L(i,1)*H(1,2:4);
   H(i,1) = 0;
end
disp('after one stage of elimination:')
disp('  L =')
disp(L)
disp('  H =')
disp(H)
pause
%    second stage
for i=3:4
   L(i,2) = H(i,2)/H(2,2);
   H(i,3:4) = H(i,3:4) - L(i,2)*H(2,3:4);
   H(i,2) = 0;
end
disp('after second stage of elimination:')
disp('  L =')
disp(L)
disp('  H =')
disp(H)
pause
%    swap rows 3 and 4
tmp = H(3,:);
H(3,:) = H(4,:);
H(4,:) = tmp;
tmp = L(3,1:2);
L(3,1:2) = L(4,1:2);
L(4,1:2) = tmp;
disp('after swapping rows 3 and 4:')
disp('  L =')
disp(L)
disp('  H =')
disp(H)
pause
%    third stage
L(4,3) = H(4,3)/H(3,3);
H(4,4) = H(4,4) - L(4,3)*H(3,4);
H(4,3) = 0;
disp('after last stage of elimination:')
disp('  L =')
disp(L)
disp('  H =')
disp(H)
pause

% store results
L16 = L;
H16 = H;


% Second, do the 4x4 Hilbert matrix using 3-digit precision
%    set up matrix, display
format long
H = hilb(4);
H = roundsd(H,d);
disp('Repeating with 3-digit precision; original matrix:')
disp(H)
pause
%    first stage
L = eye(4);
for i=2:4
   L(i,1) = roundsd(H(i,1)/H(1,1),d);
   H(i,2:4) = roundsd(H(i,2:4) - roundsd(L(i,1)*H(1,2:4),d),d);
   H(i,1) = 0;
end
disp('after one stage of elimination:')
disp('  L =')
disp(L)
disp('  H =')
disp(H)
pause
%    second stage
for i=3:4
   L(i,2) = roundsd(H(i,2)/H(2,2),d);
   H(i,3:4) = roundsd(H(i,3:4) - roundsd(L(i,2)*H(2,3:4),d),d);
   H(i,2) = 0;
end
disp('after second stage of elimination:')
disp('  L =')
disp(L)
disp('  H =')
disp(H)
pause
%    swap rows 3 and 4
tmp = H(3,:);
H(3,:) = H(4,:);
H(4,:) = tmp;
tmp = L(3,1:2);
L(3,1:2) = L(4,1:2);
L(4,1:2) = tmp;
disp('after swapping rows 3 and 4:')
disp('  L =')
disp(L)
disp('  H =')
disp(H)
pause
%    third stage
L(4,3) = roundsd(H(4,3)/H(3,3),d);
H(4,4) = roundsd(H(4,4) - roundsd(L(4,3)*H(3,4),d),d);
H(4,3) = 0;
disp('after last stage of elimination:')
disp('  L =')
disp(L)
disp('  H =')
disp(H)
pause

H3 = H;
L3 = L;


% now compare results
disp('We now compare the results of these two calculations:')
format
disp('L matrices:')
disp('  16-digit:')
disp(L16)
disp('  3-digit:')
disp(L3)
disp('  relative error:')
err = zeros(4,4);
for i=1:4
   for j=1:i
      err(i,j) = abs(L16(i,j)-L3(i,j))/abs(L16(i,j));
   end
end
disp(err)
pause

disp('U matrices:')
disp('  16-digit:')
disp(H16)
disp('  3-digit:')
disp(H3)
disp('  relative error:')
err = zeros(4,4);
for i=1:4
   for j=i:4
      err(i,j) = abs(H16(i,j)-H3(i,j))/abs(H16(i,j));
   end
end
disp(err)
