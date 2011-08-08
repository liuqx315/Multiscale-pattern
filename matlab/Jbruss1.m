function J = Jbruss1(t, y)
% usage: J = Jbruss1(t, y)
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

% form the ODE Jacobian
J = [-(y(3)+1) + 2*y(1)*y(2),   y(1)*y(1),     -y(1);
        y(3) - 2*y(1)*y(2),    -y(1)*y(1),      y(1);
             -y(3),               0,       -1/ep - y(1)];

% end function
