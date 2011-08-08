function dy = fbruss1(t, y)
% usage: dy = fbruss1(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% model parameters
% test 1 -- initially very fast, then slower evolution
a = 1.2; 
b = 2.5; 
ep = 1e-5;

% form the ODE RHS
dy = [a - (y(3)+1)*y(1) + y(1)*y(1)*y(2);
      y(3)*y(1) - y(1)*y(1)*y(2);
      (b-y(3))/ep - y(3)*y(1)];

% end function
