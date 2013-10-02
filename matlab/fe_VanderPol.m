function dy = fe_VanderPol(t, y)
% usage: dy = fe_VanderPol(t, y)
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

% form the ODE RHS
du = v;
dv = 0;

dy = [du; dv];

% end function
