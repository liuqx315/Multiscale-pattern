function J = Ji_VanderPol(t, y)
% usage: J = Ji_VanderPol(t, y)
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
ep = Pdata.ep;
ep2 = Pdata.ep2;

% extract variables
u = y(1);
v = y(2);

% form the ODE Jacobian
Juu = 0;
Juv = 0;
Jvu = -ep2/ep;
Jvv = (1-u^2)/ep;
J = [Juu, Juv; Jvu, Jvv];

% end function
