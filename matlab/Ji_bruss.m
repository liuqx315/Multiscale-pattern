function J = Ji_bruss(t, y)
% usage: J = Ji_bruss(t, y)
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

% form the ODE Jacobian
J = [0, 0, 0; 0, 0, 0; 0, 0, -1/ep];

% end function
