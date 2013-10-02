function J = J_bruss(t, y)
% usage: J = J_bruss(t, y)
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% model parameters
global Pdata;
a = Pdata.a; 
b = Pdata.b; 
ep = Pdata.ep;

% form the ODE Jacobian
J = [-(y(3)+1) + 2*y(1)*y(2),   y(1)*y(1),     -y(1);
        y(3) - 2*y(1)*y(2),    -y(1)*y(1),      y(1);
             -y(3),               0,       -1/ep - y(1)];

% end function
