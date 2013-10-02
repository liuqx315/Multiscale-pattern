function J = J_analytic(t, y)
% usage: J = J_analytic(t, y)
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

% form the ODE Jacobian
J = lam;

% end function
