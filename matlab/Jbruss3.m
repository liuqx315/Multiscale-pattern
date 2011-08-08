function J = Jbruss3(t, y)
% usage: J = Jbruss3(t, y)
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

% form the ODE Jacobian
J = [-(y(3)+1) + 2*y(1)*y(2),   y(1)*y(1),     -y(1);
        y(3) - 2*y(1)*y(2),    -y(1)*y(1),      y(1);
             -y(3),               0,       -1/ep - y(1)];

% end function
