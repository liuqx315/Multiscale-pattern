function dy = fbruss3(t, y)
% usage: dy = fbruss3(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% model parameters
% test3 -- initially very rapid, then equilibrium or steady growth
a = 0.5; 
b = 3; 
ep = 5e-4;

% form the ODE RHS
dy = [a - (y(3)+1)*y(1) + y(1)*y(1)*y(2);
      y(3)*y(1) - y(1)*y(1)*y(2);
      (b-y(3))/ep - y(3)*y(1)];

% end function
