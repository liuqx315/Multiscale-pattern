function dy = f_analytic(t, y)
% usage: dy = f_analytic(t, y)
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
lam = -1/ep;

% form the ODE RHS
dy = lam*y + 1/(1+t^2) - lam*atan(t);

% end function
