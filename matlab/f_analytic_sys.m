function dy = f_analytic_sys(t, y)
% usage: dy = f_analytic_sys(t, y)
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
lam = Pdata.lam;

% form the ODE RHS
V = [1 -1 1; -1 2 1; 0 -1 2];
Vinv = 0.25*[5 1 -3; 2 2 -2; 1 1 1];
D = [-0.5 0 0; 0 -0.1 0; 0 0 lam];
dy = V*D*Vinv*y;

% end function
